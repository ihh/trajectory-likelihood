#include <random>

#include "trajec.h"

namespace TrajectoryLikelihood {
  struct SimulationConfig {
    int initSeqLen;
    int maxGapLen;
    int trials;
    int verbose;
    SimulationConfig() : initSeqLen(1000), maxGapLen(100), trials(100), verbose(0) { }
    SimulationConfig (int i, int g, int t, int v) : initSeqLen(i), maxGapLen(g), trials(t), verbose(v) { }
  };

  struct SimulationCounts {
    vector<vector<int> > count;
    int overflowCount;
    SimulationCounts (int len) : count (len + 1, vector<int> (len + 1, 0)), overflowCount(0) { }
  };
  
  SimulationCounts chopZoneSimulatedCounts (const IndelParams& params, double time, const SimulationConfig& config, mt19937& rnd);   // result is indexed [nDeleted][nInserted]
  vector<vector<double> > chopZoneSimulatedProbabilities (const IndelParams& params, double time, const SimulationConfig& config, mt19937& rnd);
  void simulateIndels (vector<int>& seq, const IndelParams& params, double time, mt19937& rnd, int verbose, int insertedValue = 0);
}
