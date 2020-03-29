#include "trajec.h"

// private function prototypes
double factorial (int n);
bool doublesEqual (double a, double b);
double indelTrajectoryLikelihood (const vector<int>& zoneLengths, const IndelParams& params, double time);
bool anyIdenticalNeighbors (const vector<int>& list);
double exitRateForZoneLength (int zoneLength, const IndelParams& params);
vector<double> exitRatesForZoneLengths (const vector<int>& zoneLengths, const IndelParams& params);
double insertionRate (int k, const IndelParams& params);
double deletionRate (int k, const IndelParams& params);
double transitionRateForZoneLengthChange (int srcZoneLength, int destZoneLength, const IndelParams& params);
vector<double> transitionRatesForZoneLengths (const vector<int>& zoneLengths, const IndelParams& params);
bool trajectoryIsValid (const vector<int>& zoneLengths, int nInserted, int nDeleted);
int indelTrajectoryDegeneracy (const vector<int>& zoneLengths);
int fastIndelTrajectoryDegeneracy (const vector<int>& zoneLengths);
vector<Machine> makeEventMachines (const vector<int>& zoneLengths);
Machine makeTrajectoryMachine (const vector<int>& zoneLengths);
Machine makeInsertionMachine (int k);
Machine makeDeletionMachine (int k);
int countTotalInsertions (const vector<int>& zoneLengths);
int countTotalDeletions (const vector<int>& zoneLengths);
void logTrajectoryLikelihood (const vector<int>& zoneLengths, double pTraj, const IndelParams& params, double time, const ChopZoneConfig& config, bool isValid);
  
// function definitions
bool doublesEqual (double a, double b) {
  const double diff = abs (a - b);
  return diff < std::numeric_limits<double>::epsilon();
}

// Algorithm 1 of (Miklos, Lunter & Holmes, 2004)
// exitRates = chi
// transitionRates = r
double trajectoryLikelihood (const vector<double>& exitRates,
			     const vector<double>& transitionRates,
			     double time) {
  if (time <= 0)
    throw std::runtime_error ("Trajectory must have finite duration");
  if (exitRates.size() == 0)
    throw std::runtime_error ("There must be at least one state in the trajectory");
  if (exitRates.size() != transitionRates.size() + 1)
    throw std::runtime_error ("Mismatch between numbers of states & transition rates");
  for (int i = 0; i < transitionRates.size(); ++i) {
    if (transitionRates[i] < 0)
      throw std::runtime_error ("Transition rates must all be nonnegative");
    if (transitionRates[i] > exitRates[i])
      throw std::runtime_error ("Exit rates must be at least as great as transition rates");
  }

  // switch to the notation from the paper
  const vector<double>& zeta = exitRates;
  const vector<double>& r = transitionRates;
  const double T = time;
  vector<double> chi (1, exitRates[0]);  // unique exit rates
  vector<double> d (1, 0.);  // exit rate degeneracy
  vector<vector<double> > c (1, vector<double> (1, 1.));
  for (int i = 1; i < zeta.size(); ++i) {
    int j = -1;
    for (int j_test = 0; j_test < chi.size(); ++j_test)
      if (doublesEqual (chi[j_test], zeta[i])) {
	j = j_test;
	break;
      }
    if (j < 0) {
      chi.push_back (zeta[i]);
      d.push_back (0);
    } else
      ++d[j];
    vector<vector<double> > u;  // u[n][k] = new value of c[n][k]
    const int M = chi.size() - 1;
    for (int n = 0; n < chi.size(); ++n) {
      vector<double> u_n;
      for (int k = 0; k <= d[n]; ++k) {
        double u_nk = 0;
        if (!doublesEqual (chi[n], zeta[i]))
          for (int j = k; j <= d[n]; ++j)
            u_nk -= c[n][j] * factorial(j) / (factorial(k) * pow (chi[n] - zeta[i], j - k + 1));
        else if (k == 0) {
          for (int m = 0; m <= M; ++m)
            if (m != n)
              for (int j = k; j <= d[m]; ++j)
                u_nk += c[m][j] * factorial(j) / pow (chi[m] - zeta[i], j+1);
        } else
          u_nk = c[n][k-1] / k;
        u_n.push_back (u_nk);
      }
      u.push_back (u_n);
    }
    swap (c, u);
  }
  double P = 0;
  const int M = chi.size() - 1;
  for (int n = 0; n <= M; ++n) {
    double T_poly = 0;
    for (int k = 0; k <= d[n]; ++k)
      T_poly += c[n][k] * pow (T, k);
    P += exp (-chi[n]*T) * T_poly;
  }
  for (double rate: r)
    P *= rate;
  return P;
}

vguard<double> precomputedFactorial (1, 1);
double factorial (int n) {
  if (n < 0)
    throw std::runtime_error ("Factorial function defined for nonnegative integers only");
  for (int k = precomputedFactorial.size(); k <= n; ++k)
    precomputedFactorial.push_back (precomputedFactorial[k-1] * k);
  return precomputedFactorial[n];
}

// The trajectory in Figure 1 of MLH2004 is
//  AAAM -> AABBBBAM -> BBAM -> BBM
// The zone lengths for this trajectory are
//  3 -> 7 -> 3 -> 2
// Note that there are other valid trajectories with these zone lengths, e.g.
//  AAAM -> AAABBBBM -> ABBM -> BBM
// There are also invalid trajectories with the same zone lengths, e.g.
//  AAAM -> AAABBBBM -> AAAM -> AAM  (not allowed; we can't have any A's left at the end)
// We count the number of valid trajectories using a finite-state machine approach.
double indelTrajectoryLikelihood (const vector<int>& zoneLengths, const IndelParams& params, double time) {
  const int N = zoneLengths.size();
  if (time < 0)
    throw std::runtime_error ("Time must be nonnegative");
  if (N == 0)
    throw std::runtime_error ("There must be at least one state in the trajectory");
  for (int n = 0; n < zoneLengths.size(); ++n) {
    if (zoneLengths[n] < 1)
      throw std::runtime_error ("All zone lengths must be positive integers");
  }
  if (anyIdenticalNeighbors (zoneLengths))
    throw std::runtime_error ("No two adjacent states in the trajectory can be identical");
  const vector<double> exitRates = exitRatesForZoneLengths (zoneLengths, params);
  const vector<double> transitionRates = transitionRatesForZoneLengths (zoneLengths, params);
  return trajectoryLikelihood (exitRates, transitionRates, time) * fastIndelTrajectoryDegeneracy (zoneLengths);
}

vector<double> exitRatesForZoneLengths (const vector<int>& zoneLengths, const IndelParams& params) {
  vector<double> exitRates;
  for (int n = 0; n < zoneLengths.size(); ++n)
    exitRates.push_back (exitRateForZoneLength (zoneLengths[n], params));
  return exitRates;
}

vector<double> transitionRatesForZoneLengths (const vector<int>& zoneLengths, const IndelParams& params) {
  vector<double> transitionRates;
  for (int n = 1; n < zoneLengths.size(); ++n)
    transitionRates.push_back (transitionRateForZoneLengthChange (zoneLengths[n-1], zoneLengths[n], params));
  return transitionRates;
}

bool anyIdenticalNeighbors (const vector<int>& list) {
  for (int n = 0; n + 1 < list.size(); ++n)
    if (list[n] == list[n+1])
      return true;
  return false;
}

double exitRateForZoneLength (int zoneLength, const IndelParams& params) {
  const double lambdaSum = params.gamma * params.mu * (1-params.r) * (1-params.r) / (1 - params.gamma * params.r);
  const double muSum = params.mu * (1 - params.r);
  return zoneLength * (lambdaSum + muSum);
}

double insertionRate (int k, const IndelParams& params) {
  return params.gamma * params.mu * (1 - params.r) * (1 - params.r) * pow (params.gamma * params.r, k-1);
}

double deletionRate (int k, const IndelParams& params) {
  return params.mu * (1-params.r) * (1-params.r) * pow (params.r, k-1);
}

double transitionRateForZoneLengthChange (int srcZoneLength, int destZoneLength, const IndelParams& params) {
  return (destZoneLength > srcZoneLength
          ? insertionRate (destZoneLength - srcZoneLength, params)
          : deletionRate (srcZoneLength - destZoneLength, params));
}

const vector<string> wildChars = { string("A"), string("B") };
const Machine wildEcho = Machine::wildEcho (wildChars);
int indelTrajectoryDegeneracy (const vector<int>& zoneLengths) {
  if (zoneLengths.size() == 1)  // catch the special case where there is no machine to create
    return 1;
  const Machine machine = makeTrajectoryMachine (zoneLengths);
  const EvaluatedMachine eval (machine);
  SeqPair seqPair;
  seqPair.input.seq = vector<string> (zoneLengths[0] - 1, string("A"));
  seqPair.output.seq = vector<string> (zoneLengths[zoneLengths.size() - 1] - 1, string("B"));
  const RollingOutputForwardMatrix fwd (eval, seqPair);
  return (int) (.5 + exp (fwd.logLike()));
}

Machine makeTrajectoryMachine (const vector<int>& zoneLengths) {
  Machine machine = wildEcho;
  const vector<Machine> eventMachines = makeEventMachines (zoneLengths);
  for (const auto& eventMachine: eventMachines)
    machine = Machine::compose (machine, eventMachine);
  return machine;
}

vector<Machine> makeEventMachines (const vector<int>& zoneLengths) {
  vector<Machine> machines;
  for (int n = 0; n + 1 < zoneLengths.size(); ++n) {
    const int deltaLength = zoneLengths[n+1] - zoneLengths[n];
    const Machine deltaMachine = (deltaLength > 0
				  ? makeInsertionMachine (deltaLength)
				  : makeDeletionMachine (-deltaLength));
    machines.push_back (deltaMachine);
  }
  return machines;
}

vector<Machine> insertionMachine (1, Machine::null());
Machine makeInsertionMachine (int k) {
  while (insertionMachine.size() <= k)
    insertionMachine.push_back (Machine::concatenate
				(Machine::concatenate (wildEcho,
						       Machine::generator (vector<string> (insertionMachine.size(),
											   string("B")))),
				 wildEcho));
  return insertionMachine[k];
}

vector<Machine> deletionMachine (1, Machine::null());
const Machine singleDelete = Machine::wildSingleRecognizer (wildChars);
Machine makeDeletionMachine (int k) {
  while (deletionMachine.size() <= k)
    deletionMachine.push_back (Machine::concatenate
			       (Machine::concatenate (wildEcho,
						      Machine::repeat (singleDelete, deletionMachine.size())),
				wildEcho));
  return deletionMachine[k];
}

// Calculate chop zone probabilities (internal zones i.e. not at the ends of the sequence)
double chopZoneLikelihood (int nDeleted, int nInserted, const IndelParams& params, double time, const ChopZoneConfig& config) {
  double prob = 0;
  int minEvents = 0;
  if (nDeleted)
    ++minEvents;
  if (nInserted)
    ++minEvents;
  for (int events = minEvents; events <= config.maxEvents; ++events) {
    if (events == 0) {  // only true if nDeleted == nInserted == 0
      const vector<int> zoneLengths (1, 1);
      const double pTraj = indelTrajectoryLikelihood (zoneLengths, params, time);
      prob += pTraj;
      logTrajectoryLikelihood (zoneLengths, pTraj, params, time, config, true);
    } else {  // events > 0
      vector<int> zoneLengths (events + 1, 1);
      zoneLengths[0] = nDeleted + 1;
      zoneLengths[events] = nInserted + 1;
      bool finished = false;
      while (!finished) {
	const bool isValid = trajectoryIsValid (zoneLengths, nInserted, nDeleted);
        const double pTraj = (isValid
			      ? indelTrajectoryLikelihood (zoneLengths, params, time)
			      : 0);
	prob += pTraj;
	logTrajectoryLikelihood (zoneLengths, pTraj, params, time, config, isValid);
        if (events == 1)
          finished = true;
        else
          for (int i = 1; true; ++i)
            if (++zoneLengths[i] > config.maxLen) {
              if (i == events - 1) {
                finished = true;
                break;
              } else
                zoneLengths[i] = 1;
            } else
              break;
      }
    }
  }
  return prob;
}

void logTrajectoryLikelihood (const vector<int>& zoneLengths, double pTraj, const IndelParams& params, double time, const ChopZoneConfig& config, bool isValid) {
  if (config.verbose && ((isValid && pTraj > 0) || config.verbose > 1)) {
    cerr << "Zone lengths: [" << to_string_join(zoneLengths) << "]  ";
    if (isValid) {
      if (config.verbose > 2) {
	const auto exitRates = exitRatesForZoneLengths (zoneLengths, params);
	cerr << "Exit rates: [" << to_string_join(exitRates) << "]  ";
	if (config.verbose > 3) {
	  const auto transitionRates = transitionRatesForZoneLengths (zoneLengths, params);
	  cerr << "Transition rates: [" << to_string_join(transitionRates) << "]  ";
	  if (config.verbose > 4) {
	    const auto trajLike = trajectoryLikelihood (exitRates, transitionRates, time);
	    const auto degen = indelTrajectoryDegeneracy (zoneLengths);
	    cerr << "P(traj) " << trajLike << " #traj=" << degen << "  ";
	  }
	}
      }
      cerr << "Probability: " << pTraj << endl;
      if (isValid && pTraj > 0 && config.verbose > 5) {
	if (config.verbose > 6)
	  for (const auto& e: makeEventMachines (zoneLengths))
	    e.writeJson (cerr);
	makeTrajectoryMachine (zoneLengths).writeJson (cerr);
      }
    } else
      cerr << "Invalid" << endl;
  }
}

bool trajectoryIsValid (const vector<int>& zoneLengths, int nInserted, int nDeleted) {
  const int trajectoryInsertions = countTotalInsertions (zoneLengths);
  const int trajectoryDeletions = countTotalDeletions (zoneLengths);
  return !anyIdenticalNeighbors (zoneLengths) && trajectoryInsertions >= nInserted && trajectoryDeletions >= nDeleted;
}

int countTotalInsertions (const vector<int>& zoneLengths) {
  int count = 0;
  for (int i = 1; i < zoneLengths.size(); ++i)
    count += max (zoneLengths[i] - zoneLengths[i-1], 0);
  return count;
}

int countTotalDeletions (const vector<int>& zoneLengths) {
  int count = 0;
  for (int i = 1; i < zoneLengths.size(); ++i)
    count -= min (zoneLengths[i] - zoneLengths[i-1], 0);
  return count;
}

vector<vector<double> > chopZoneLikelihoods (const IndelParams& params, double time, const ChopZoneConfig& config) {
  vector<vector<double> > probs;
  for (int nDeleted = 0; nDeleted <= config.maxLen; ++nDeleted) {
    vector<double> pd;
    for (int nInserted = 0; nInserted <= config.maxLen; ++nInserted)
      pd.push_back (chopZoneLikelihood (nDeleted, nInserted, params, time, config));
    probs.push_back (pd);
  }
  return probs;
}

int fastIndelTrajectoryDegeneracy (const vector<int>& zoneLengths) {
  int result = 0;
  bool verified = false;
  switch (zoneLengths.size()) {
  case 1:
  case 2:
    result = 1;
    break;
  case 3:
    {
      const int z0 = zoneLengths[0],
	z1 = zoneLengths[1],
	z2 = zoneLengths[2];
      const bool ins1 = z1 > z0, ins2 = z2 > z1;
      const int e1 = abs(z1 - z0), e2 = abs(z2 - z1);
      if (ins1 && ins2)  // ins,ins
	result = z0 * z1;
      else if (ins1 && !ins2)  // ins,del
	result = z2 > 1 ? (z0 > 1 ? 2 : (z1 - e2)) : z0;
      else if (!ins1 && ins2)  // del,ins
	result = 1;
      else  // del,del
	result = z0 - e1;
    }
    break;
  case 4:
    {
      const int z0 = zoneLengths[0],
	z1 = zoneLengths[1],
	z2 = zoneLengths[2],
	z3 = zoneLengths[3];
      const bool ins1 = z1 > z0, ins2 = z2 > z1, ins3 = z3 > z2;
      const int e1 = abs(z1 - z0), e2 = abs(z2 - z1), e3 = abs(z3 - z2);
      if (ins1 && ins2 && ins3) {  // ins,ins,ins
	result = z0 * z1 * z2;
      } else if (ins1 && ins2 && !ins3) {  // ins,ins,del
	if (z3 > 1) {
	  if (z0 > 1) {
	    result = 2 * ((e1 + 1) + (z2 - e3));
	  } else
	    result = (e1 + 1) * (z2 - e3);
	} else
	  result = (z0 * z1) * (z2 - e3);
      } else if (ins1 && !ins2 && ins3) {  // ins,del,ins
	result = (z2 > 1 ? (z0 > 1 ? 2 : (z1 - e2)) : z0) * z2;
      } else if (ins1 && !ins2 && !ins3) {  // ins,del,del
	if (z3 > 1) {
	  if (z0 > 1) {
	    if (e2 <= e1)
	      result = 2 * ((e1 + 1 - e2) + (z2 - e3));
	    else
	      result = 2 * (z1 - e2) * (z2 - e3);
	  } else
	    result = (z1 - e2) * (z2 - e3);
	} else
	  result = z0 * (z1 - e2) * (z2 - e3);
      } else if (!ins1 && ins2 && ins3) {  // del,ins,ins
	result = z1 * z2;
      } else if (!ins1 && ins2 && !ins3) {  // del,ins,del
	result = (z0 - e1) * (z3 > 1 ? (z1 > 1 ? 2 : (z2 - e3)) : z1);
      } else if (!ins1 && !ins2 && ins3) {  // del,del,ins
	result = z0 - e1;
      } else {  // del,del,del
	result = (z0 - e1) * (z1 - e2);
      }
    }
    break;
  default:
    result = indelTrajectoryDegeneracy (zoneLengths);
    verified = true;
    break;
  }
  if (!verified) {
    const int expected = indelTrajectoryDegeneracy (zoneLengths);
    if (result != expected) {
      cerr << "Zone lengths: " << to_string_join (zoneLengths) << endl;
      cerr << "Expected degeneracy: " << expected << endl;
      cerr << "Calculated degeneracy: " << result << endl;
      throw std::runtime_error ("Fast degeneracy calculation failed");
    }
  }
  return result;
}
