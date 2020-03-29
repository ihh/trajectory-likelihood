#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "trajec.h"

namespace TrajectoryLikelihood {

#define Assert(assertion,...) do { if (!(assertion)) Abort("Assertion Failed: " __VA_ARGS__); } while (0)

// Private function prototypes
void Abort(const char* error, ...);
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
int indelTrajectoryDegeneracy (const vector<int>& zoneLengths, bool lastResidueConserved = true);  // lastResidueConserved is true for chop zones inside the sequence & at the left end
bool runLengthEncodedSequenceHasNoAncestralResidues (const vector<int>& seq);
bool runLengthEncodedSequenceIsValid (const vector<int>& seq);
int runLengthEncodedSequenceLength (const vector<int>& seq);
void appendToRunLengthEncodedSequence (vector<int>& seq, int chunk);
void mutateRunLengthEncodedSequence (const vector<int>& ancestor, int pos, int delta, vector<int>& descendant, int expectedLen);
int countTotalInsertions (const vector<int>& zoneLengths);
int countTotalDeletions (const vector<int>& zoneLengths);
void logTrajectoryLikelihood (const vector<int>& zoneLengths, double pTraj, const IndelParams& params, double time, const ChopZoneConfig& config, bool isValid);

// Templates
// sgn
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Function definitions
void Abort(const char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  fprintf(stderr,"Abort: ");
  vfprintf(stderr,error,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
  throw runtime_error("Abort");
}

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

vector<double> precomputedFactorial (1, 1);
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
// We count the number of valid trajectories by enumeration, using run-length encoding to compress the sequences.
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
  return trajectoryLikelihood (exitRates, transitionRates, time) * indelTrajectoryDegeneracy (zoneLengths);
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
  if (config.verbose)
    cerr << "Likelihood for (" << nDeleted << " deletions, " << nInserted << " insertions) is " << prob << endl;
  return prob;
}

void logTrajectoryLikelihood (const vector<int>& zoneLengths, double pTraj, const IndelParams& params, double time, const ChopZoneConfig& config, bool isValid) {
  if (config.verbose > 1 && ((isValid && pTraj > 0) || config.verbose > 2)) {
    cerr << "Zone lengths: [" << to_string_join(zoneLengths) << "]  ";
    if (isValid) {
      if (config.verbose > 3) {
	const auto exitRates = exitRatesForZoneLengths (zoneLengths, params);
	cerr << "Exit rates: [" << to_string_join(exitRates) << "]  ";
	if (config.verbose > 4) {
	  const auto transitionRates = transitionRatesForZoneLengths (zoneLengths, params);
	  cerr << "Transition rates: [" << to_string_join(transitionRates) << "]  ";
	  if (config.verbose > 5) {
	    const auto trajLike = trajectoryLikelihood (exitRates, transitionRates, time);
	    const auto degen = indelTrajectoryDegeneracy (zoneLengths);
	    cerr << "P(traj) " << trajLike << " #traj=" << degen << "  ";
	  }
	}
      }
      cerr << "Probability: " << pTraj << endl;
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

int runLengthEncodedSequenceLength (const vector<int>& seq) {
  int len = 0;
  for (auto x: seq)
    len += abs(x);
  return len;
}

bool runLengthEncodedSequenceIsValid (const vector<int>& seq) {
  for (int i = 0; i + 1 < seq.size(); ++i)
    if (sgn(seq[i]) == sgn(seq[i+1]))
      return false;
  return true;
}

void appendToRunLengthEncodedSequence (vector<int>& seq, int chunk) {
  if (seq.size() && sgn(chunk) == sgn(seq.back()))
    seq.back() += chunk;
  else
    seq.push_back (chunk);
}

void mutateRunLengthEncodedSequence (const vector<int>& ancestor, int pos, int delta, vector<int>& descendant, int expectedLen) {
  //  cerr << "Mutating (" << to_string_join(ancestor) << ") at " << pos << " by " << delta << endl;
  Assert (runLengthEncodedSequenceIsValid(ancestor), "invalid sequence!");
  Assert (pos + max(-delta,0) <= runLengthEncodedSequenceLength(ancestor), "overflow!");
  descendant.clear();
  int idx = 0;
  while (pos > 0) {
    const int nextChunkSize = abs (ancestor[idx]);
    if (pos < nextChunkSize)
      break;
    descendant.push_back (ancestor[idx++]);
    pos -= nextChunkSize;
  }
  if (delta > 0) {  // insertion
    if (pos > 0 && idx < ancestor.size())
      descendant.push_back (sgn(ancestor[idx]) * pos);
    appendToRunLengthEncodedSequence (descendant, -delta);
    if (idx < ancestor.size()) {
      appendToRunLengthEncodedSequence (descendant, sgn(ancestor[idx]) * (abs(ancestor[idx]) - pos));
      ++idx;
    }
  } else {  // deletion
    int delSize = -delta, nextChunkSize = abs(ancestor[idx]);
    if (pos > 0) {
      descendant.push_back (sgn(ancestor[idx]) * pos);
      nextChunkSize -= pos;
    }
    while (delSize > 0 && delSize >= nextChunkSize) {
      delSize -= nextChunkSize;
      nextChunkSize = ++idx >= ancestor.size() ? 0 : abs(ancestor[idx]);
    }
    if (nextChunkSize)
      appendToRunLengthEncodedSequence (descendant, sgn(ancestor[idx++]) * (nextChunkSize - delSize));
  }
  while (idx < ancestor.size())
    appendToRunLengthEncodedSequence (descendant, ancestor[idx++]);
  //  cerr << "Mutated (" << to_string_join(ancestor) << ") at " << pos << " by " << delta << " yielding (" << to_string_join(descendant) << ")" << endl;
  Assert (runLengthEncodedSequenceLength(descendant) == expectedLen, "length is wrong");
}

bool runLengthEncodedSequenceHasNoAncestralResidues (const vector<int>& seq) {
  return seq.size() == 0 || (seq.size() == 1 && seq[0] < 0);
}

int indelTrajectoryDegeneracy (const vector<int>& zoneLengths, bool lastResidueConserved) {
  //  cerr << "Calculating degeneracy for zone lengths (" << to_string_join(zoneLengths) << ")" << endl;
  const int conservedResidues = lastResidueConserved ? 1 : 0;
  if (zoneLengths.size() == 1 && zoneLengths[0] == conservedResidues)
    return 1;
  int result = 0;
  const int nEvents = zoneLengths.size() - 1;
  vector<int> delta (nEvents), eventPos (nEvents, 0), maxEventPos (nEvents);
  for (int i = 0; i < nEvents; ++i) {
    delta[i] = zoneLengths[i+1] - zoneLengths[i];
    maxEventPos[i] = zoneLengths[i] + min (delta[i], 0) - 1;
  }
  vector<vector<int> > zoneSeq (zoneLengths.size());
  if (zoneLengths[0] > conservedResidues)
    zoneSeq[0].push_back (zoneLengths[0] - conservedResidues);
  //  cerr << "Event offsets: (" << to_string_join(eventPos) << "), max (" << to_string_join(maxEventPos) << ")" << endl;
  for (int i = 0; i < nEvents; ++i)
    mutateRunLengthEncodedSequence (zoneSeq[i], eventPos[i], delta[i], zoneSeq[i+1], zoneLengths[i+1] - conservedResidues);
  while (true) {
    if (runLengthEncodedSequenceHasNoAncestralResidues (zoneSeq[nEvents]))
      ++result;
    int i = nEvents - 1;
    while (i >= 0)
      if (++eventPos[i] > maxEventPos[i])
	eventPos[i--] = 0;
      else
	break;
    if (i < 0)
      break;
    //    cerr << "Event offsets: (" << to_string_join(eventPos) << "), max (" << to_string_join(maxEventPos) << ")" << endl;
    for (int j = i; j < nEvents; ++j)
      mutateRunLengthEncodedSequence (zoneSeq[j], eventPos[j], delta[j], zoneSeq[j+1], zoneLengths[j+1] - conservedResidues);
  }

  return result;
}

}
