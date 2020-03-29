#include <numeric>
#include <iostream>

#include "simulate.h"

namespace TrajectoryLikelihood {

  // Private function prototypes
  vector<int> runLengthEncodeSequence (const vector<int>& seq);

  // Function definitions

  // runLengthEncodeSequence converts a sequence of 0's (for inserted residues) and +ve numbers (for ancestral residues)
  // to a run-length encoded form as in trajec.cpp
  vector<int> runLengthEncodeSequence (const vector<int>& seq) {
    vector<int> rle;
    int lastSign = 0, count = 0;
    for (int c: seq) {
      const int currentSign = (c == 0 ? -1 : +1);
      if (currentSign == lastSign)
	++count;
      else {
	if (count)
	  rle.push_back (lastSign * count);
	count = 1;
	lastSign = currentSign;
      }
    }
    if (count)
      rle.push_back (lastSign * count);
    return rle;
  }
  
  void simulateIndels (vector<int>& seq, const IndelParams& params, double time, mt19937& rnd, int verbose, int insertedValue) {
    uniform_real_distribution<double> uniform (0.0, 1.0);
    const double insRate = params.totalInsertionRatePerSite();
    const double delRate = params.totalRightwardDeletionRatePerSite();
    if (verbose > 2)
      cerr << "Initial sequence length is " << seq.size() << ", time remaining is " << time << endl;
    while (time > 0) {
      const int seqLen = seq.size();
      const double totalInsRate = insRate * (seqLen + 1), totalDelRate = delRate * seqLen;
      const double totalIndelRate = totalInsRate + totalDelRate;
      const double timeToNextEvent = -log(uniform(rnd)) / totalIndelRate;
      if (verbose > 3)
	cerr << "Sequence length is " << seq.size() << ", insertion rate " << totalInsRate << ", deletion rate " << totalDelRate << ", time to next event is " << timeToNextEvent << endl;
      time -= timeToNextEvent;
      if (time < 0)
	break;
      double r = uniform(rnd) * totalIndelRate;
      int pos = 0, k = 1;
      if (r < totalInsRate) {
	// insertion
	while (r > insRate && pos <= seqLen) {
	  r -= insRate;
	  ++pos;
	}
	while ((r -= params.insertionRate(k)) > 0 && (pos + k) < seqLen)
	  ++k;
	seq.insert (seq.begin() + pos, k, insertedValue);
	if (verbose > 2)
	  cerr << "With time " << time << " remaining, " << k << "-residue insertion at position " << pos << " increased sequence length to " << seq.size() << endl;
      } else {
	// deletion
	r -= totalInsRate;
	while (r > delRate && pos < seqLen) {
	  r -= delRate;
	  ++pos;
	}
	while ((r -= params.rightwardDeletionRate(k)) > 0 && (pos + k) < seqLen)
	  ++k;
	seq.erase (seq.begin() + pos, seq.begin() + pos + k);
	if (verbose > 2)
	  cerr << "With time " << time << " remaining, " << k << "-residue deletion at position " << pos << " reduced sequence length to " << seq.size() << endl;
      }
      if (verbose > 4)
	cerr << "Sequence: (" << to_string_join (runLengthEncodeSequence (seq)) << ")" << endl;
    }
    if (verbose > 3)
      cerr << "Final sequence: (" << to_string_join (runLengthEncodeSequence (seq)) << ")" << endl;
  }

  SimulationCounts chopZoneSimulatedCounts (const IndelParams& params, double time, const SimulationConfig& config, mt19937& rnd) {
    SimulationCounts simCounts (config.maxGapLen);
    for (int trial = 0; trial < config.trials; ++trial) {
      if (config.verbose > 1)
	cerr << "Beginning simulation trial #" << (trial + 1) << endl;
      vector<int> seq (config.initSeqLen);
      iota (seq.begin(), seq.end(), 1);
      simulateIndels (seq, params, time, rnd, config.verbose, 0);
      int pos = 0;
      while (pos < seq.size() && seq[pos] == 0)
	++pos;
      int endPos = ((int) seq.size()) - 1;
      while (endPos >= 0 && seq[endPos] == 0)
	--endPos;
      while (pos < endPos) {
	const int lastAncIndex = seq[pos];
	int nInsertions = 0;
	while (seq[++pos] == 0)
	  ++nInsertions;
	const int nDeletions = seq[pos] - lastAncIndex - 1;
	if (nInsertions <= config.maxGapLen && nDeletions <= config.maxGapLen)
	  ++simCounts.count[nInsertions][nDeletions];
	else
	  ++simCounts.overflowCount;
      }
    }
    return simCounts;
  }

  vector<vector<double> > chopZoneSimulatedProbabilities (const IndelParams& params, double time, const SimulationConfig& config, mt19937& rnd) {
    const SimulationCounts simCounts = chopZoneSimulatedCounts (params, time, config, rnd);
    vector<vector<double> > probs (config.maxGapLen + 1, vector<double> (config.maxGapLen + 1, 0.));
    int total = simCounts.overflowCount;
    for (const auto& v: simCounts.count)
      for (int c: v)
	total += c;
    for (int i = 0; i <= config.maxGapLen; ++i)
      for (int j = 0; j <= config.maxGapLen; ++j)
	probs[i][j] = simCounts.count[i][j] / (double) total;
    return probs;
  }
}
