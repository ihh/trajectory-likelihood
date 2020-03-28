#include "trajec.h"

// private function prototypes
double factorial (int n);
bool doublesEqual (double a, double b);
double indelTrajectoryLikelihood (const vector<int>& zoneLengths, const IndelParams& params, double time);
bool anyIdenticalNeighbors (const vector<int>& list);
double exitRateForZoneLength (int zoneLength, const IndelParams& params);
double insertionRate (int k, const IndelParams& params);
double deletionRate (int k, const IndelParams& params);
double transitionRateForZoneLengthChange (int srcZoneLength, int destZoneLength, const IndelParams& params);
int indelTrajectoryDegeneracy (const vector<int>& zoneLengths);
Machine makeInsertionMachine (int k);
Machine makeDeletionMachine (int k);
int countTotalInsertions (const vector<int>& zoneLengths);
int countTotalDeletions (const vector<int>& zoneLengths);

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
    if (zoneLengths[n] < 0)
      throw std::runtime_error ("All zone lengths must be nonnegative integers");
    if (n < zoneLengths.size() - 1 && zoneLengths[n] == 0)
      throw std::runtime_error ("Only the final zone length in the trajectory can be zero");
  }
  if (anyIdenticalNeighbors (zoneLengths))
    throw std::runtime_error ("No two adjacent states in the trajectory can be identical");
  vector<double> exitRates, transitionRates;
  for (int n = 0; n < zoneLengths.size(); ++n) {
    exitRates.push_back (exitRateForZoneLength (zoneLengths[n], params));
    if (n > 0)
      transitionRates.push_back (transitionRateForZoneLengthChange (zoneLengths[n-1], zoneLengths[n], params));
  }
  return trajectoryLikelihood (exitRates, transitionRates, time) * indelTrajectoryDegeneracy (zoneLengths);
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
  return (zoneLength + 1) * (lambdaSum + muSum);
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
  Machine machine = wildEcho;
  for (int n = 0; n + 1 < zoneLengths.size(); ++n) {
    const int deltaLength = zoneLengths[n+1] - zoneLengths[n];
    const Machine deltaMachine = (deltaLength > 0
				  ? makeInsertionMachine (deltaLength)
				  : makeDeletionMachine (-deltaLength));
    machine = Machine::compose (machine, deltaMachine);
  }
  const EvaluatedMachine eval (machine);
  SeqPair seqPair;
  seqPair.input.seq = vector<string> (zoneLengths[0], string("A"));
  seqPair.output.seq = vector<string> (zoneLengths[zoneLengths.size() - 1], string("B"));
  const RollingOutputForwardMatrix fwd (eval, seqPair);
  return (int) (.5 + exp (fwd.logLike()));
}

vector<Machine> insertionMachine (1, Machine::null());
Machine makeInsertionMachine (int k) {
  while (insertionMachine.size() <= k)
    insertionMachine.push_back (Machine::concatenate
				(Machine::concatenate (wildEcho,
						       Machine::wildGenerator (vector<string> (insertionMachine.size(),
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
      const double pTraj = indelTrajectoryLikelihood (vector<int> (1, 1), params, time);
      prob += pTraj;
      if (config.verbose)
	cerr << "Zone lengths: [ 1 ]  Probability: " << pTraj << endl;
    } else {  // events > 0
      vector<int> zoneLengths (events + 1, 1);
      zoneLengths[0] = nDeleted + 1;
      zoneLengths[events] = nInserted + 1;
      bool finished = false;
      while (!finished) {
        const int trajectoryInsertions = countTotalInsertions (zoneLengths);
        const int trajectoryDeletions = countTotalDeletions (zoneLengths);
        const double pTraj = (!anyIdenticalNeighbors (zoneLengths) && trajectoryInsertions >= nInserted && trajectoryDeletions >= nDeleted
			      ? indelTrajectoryLikelihood (zoneLengths, params, time)
			      : 0);
	prob += pTraj;
	if (config.verbose)
	  cerr << "Zone lengths: [ " << to_string_join (zoneLengths, ", ") << " ]  Probability: " << pTraj << endl;
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
