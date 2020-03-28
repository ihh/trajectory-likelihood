const machineboss = require ('machineboss'),
      { initPromise, runWithFilesSync } = machineboss

const extend = function() {
  let a = arguments[0]
  Array.from(arguments).slice(1).forEach ((b) => Object.keys(b).forEach ((k) => a[k] = b[k]))
  return a
}

// Algorithm 1 of (Miklos, Lunter & Holmes, 2004)
// exitRates = chi
// transitionRates = r
const trajectoryLikelihood = (exitRates, transitionRates, time) => {
  if (typeof(time) !== 'number')
    throw new Error ("Time must be numeric")
  if (time <= 0)
    throw new Error ("Trajectory must have finite duration")
  if (exitRates.length === 0)
    throw new Error ("There must be at least one state in the trajectory")
  if (exitRates.filter ((r) => typeof(r) !== 'number').length)
    throw new Error ("Exit rates must be numeric")
  if (exitRates.length !== transitionRates.length + 1)
    throw new Error ("Mismatch between numbers of exit rates & transition rates")
  if (transitionRates.filter ((r) => typeof(r) !== 'number').length)
    throw new Error ("Transition rates must be numeric")
  if (transitionRates.filter ((r) => r < 0).length)
    throw new Error ("Transition rates must all be nonnegative")
  if (transitionRates.filter ((r, i) => r > exitRates[i]).length)
    throw new Error ("Exit rates must be at least as great as transition rates")
  // switch to the notation from the paper
  let zeta = exitRates,
      r = transitionRates,
      T = time,
      chi = [exitRates[0]],  // unique exit rates
      d = [0],  // exit rate degeneracy
      c = [[1]]   // c[n][k] = coefficient of T^k exp(-chi[n]*T) in likelihood
  zeta.slice(1).forEach ((zeta_i) => {
    const j = chi.indexOf (zeta_i)
    if (j < 0) {
      chi.push (zeta_i)
      d.push (0)
    } else
      ++d[j]
    let u = [],  // u[n][k] = new value of c[n][k]
        M = chi.length - 1
    chi.forEach ((chi_n, n) => {
      let u_n = []
      for (let k = 0; k <= d[n]; ++k) {
        let u_nk = 0
        if (chi_n != zeta_i)
          for (let j = k; j <= d[n]; ++j)
            u_nk -= c[n][j] * factorial(j) / (factorial(k) * Math.pow (chi[n] - zeta_i, j - k + 1))
        else if (k == 0) {
          for (let m = 0; m <= M; ++m)
            if (m != n)
              for (let j = k; j <= d[m]; ++j)
                u_nk += c[m][j] * factorial(j) / Math.pow (chi[m] - zeta_i, j+1)
        } else
          u_nk = c[n][k-1] / k
        u_n.push (u_nk)
      }
      u.push (u_n)
    })
    c = u
  })
  let waitProb = 0,
      M = chi.length - 1
  for (let n = 0; n <= M; ++n) {
    let T_poly = 0
    for (let k = 0; k <= d[n]; ++k)
      T_poly += c[n][k] * Math.pow (T, k)
    waitProb += Math.exp (-chi[n]*T) * T_poly
  }
  return r.reduce ((P, rate) => P*rate, waitProb)
}

let precomputedFactorial = [1]
const factorial = (n) => {
  if (n < 0 || Math.floor(n) != n)
    throw new Error ("Factorial function defined for nonnegative integers only")
  for (let k = precomputedFactorial.length; k <= n; ++k)
    precomputedFactorial.push (precomputedFactorial[k-1] * k)
  return precomputedFactorial[n]
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
const indelTrajectoryLikelihood = (zoneLengths, params, time) => {
  const { gamma, mu, r } = params
  const N = zoneLengths.length
  if ([gamma,mu,r].filter ((p) => typeof(p) !== 'number').length)
    throw new Error ("Parameters {gamma,mu,r} must be defined numeric")
  if (typeof(time) !== 'number')
    throw new Error ("Time must be numeric")
  if (time < 0)
    throw new Error ("Time must be nonnegative")
  if (N === 0)
    throw new Error ("There must be at least one state in the trajectory")
  if (zoneLengths.filter ((l) => typeof(l) !== 'number').length)
    throw new Error ("Zone lengths must be numeric")
  if (zoneLengths.filter ((l) => l < 0 || Math.round(l) != l).length)
    throw new Error ("All zone lengths must be nonnegative integers")
  if (zoneLengths.filter ((l, n) => l == 0 && n < N-1).length)
    throw new Error ("Only the final zone length in the trajectory can be zero")
  if (countIdenticalNeighbors (zoneLengths))
    throw new Error ("No two adjacent states in the trajectory can be identical")
  const exitRates = zoneLengths.map ((l) => exitRateForZoneLength (l, params))
  const transitionRates = zoneLengths.slice(1).map ((l, n) => transitionRateForZoneLengthChange (zoneLengths[n], l, params))
  return trajectoryLikelihood (exitRates, transitionRates, time) * indelTrajectoryDegeneracy (zoneLengths)
}

const countIdenticalNeighbors = (list) => {
  return list.filter ((l, n) => n < list.length - 1 && l == list[n+1]).length
}

const exitRateForZoneLength = (zoneLength, params) => {
  const { gamma, mu, r } = params
  const lambdaSum = gamma * mu * (1-r) * (1-r) / (1 - gamma*r)
  const muSum = mu * (1-r)
  return (zoneLength + 1) * (lambdaSum + muSum)
}

const insertionRate = (k, params) => {
  const { gamma, mu, r } = params
  return gamma * mu * (1-r) * (1-r) * Math.pow (gamma*r, k-1)
}

const deletionRate = (k, params) => {
  const { gamma, mu, r } = params
  return mu * (1-r) * (1-r) * Math.pow (r, k-1)
}

const transitionRateForZoneLengthChange = (srcZoneLength, destZoneLength, params) => {
  const { gamma, mu, r } = params
  return (destZoneLength > srcZoneLength
          ? insertionRate (destZoneLength - srcZoneLength, params)
          : deletionRate (srcZoneLength - destZoneLength, params))
}

const indelTrajectoryDegeneracy = (zoneLengths) => {
  if (zoneLengths.length == 1)  // catch the special case where there is no machine to create
    return 1
  let machines = []
  for (let n = 0; n < zoneLengths.length - 1; ++n) {
    const deltaLength = zoneLengths[n+1] - zoneLengths[n]
    machines.push (deltaLength > 0
                   ? makeInsertionMachine (deltaLength)
                   : makeDeletionMachine (-deltaLength))
  }
  const inputSeq = "A".repeat (zoneLengths[0])
  const outputSeq = "B".repeat (zoneLengths[zoneLengths.length - 1])
  const args = machines.join(" --compose ").split(" ")
        .concat (["--input-chars", inputSeq,
                  "--output-chars", outputSeq,
                  "--loglike"])
  console.warn("boss "+args.join(" "))
  const mbResult = runWithFilesSync (args)
  const mbJson = JSON.parse (mbResult.stdout),
        logForwardProb = mbJson[0][2]
  return logForwardProb == '-Infinity' ? 0 : Math.round(Math.exp(logForwardProb))
}

const makeInsertionMachine = (k) => {
  return "--begin --echo-wild AB --concat --generate-chars " + "B".repeat(k) + " --concat --echo-wild AB --end"
}

const makeDeletionMachine = (k) => {
  return "--begin --echo-wild AB --concat --begin --recognize-one AB --repeat " + k + " --end --concat --echo-wild AB --end"
}

// Config for chop zone probability calculations
const defaultChopZoneLikelihoodConfig = { maxEvents: 3,
                                          maxLen: 100 }
const getConfig = (config) => extend ({},
                                      defaultChopZoneLikelihoodConfig,
                                      config || {})

// Calculate chop zone probabilities (internal zones i.e. not at the ends of the sequence)
const chopZoneLikelihood = (nDeleted, nInserted, params, time, config) => {
  const { maxEvents, maxLen } = getConfig (config)
  let prob = 0
  let minEvents = 0
  if (nDeleted)
    ++minEvents
  if (nInserted)
    ++minEvents
  for (let events = minEvents; events <= maxEvents; ++events) {
    if (events == 0)  // only true if nDeleted == nInserted == 0
      prob += indelTrajectoryLikelihood ([0], params, time)
    else {  // events > 0
      let zoneLengths = new Array (events + 1).fill (1)
      zoneLengths[0] = nDeleted + 1
      zoneLengths[events] = nInserted + 1
      let finished = false
      while (!finished) {
        console.warn(zoneLengths)
        if (!countIdenticalNeighbors (zoneLengths))
          prob += indelTrajectoryLikelihood (zoneLengths, params, time)
        if (events == 1)
          finished = true
        else
          for (let i = 1; true; ++i)
            if (++zoneLengths[i] > maxLen) {
              if (i == events - 1) {
                finished = true
                break
              } else
                zoneLengths[i] = 1
            } else
              break
      }
    }
  }
  return prob
}

const chopZoneLikelihoods = (params, time, config) => {
  const { maxEvents, maxLen } = getConfig (config)
  let probs = []
  for (let nDeleted = 0; nDeleted <= maxLen; ++nDeleted) {
    let pd = []
    for (let nInserted = 0; nInserted <= maxLen; ++nInserted)
      pd.push (chopZoneLikelihood (nDeleted, nInserted, params, time, config))
    probs.push (pd)
  }
  return probs
}

const testCZL = (params, time, config) => {
  const probs = chopZoneLikelihoods (params, time, config)
  let total = 0
  probs.forEach ((pd) => {
    pd.forEach ((p) => total += p)
    console.log (JSON.stringify (pd))
  })
  console.log ("Total: " + total)
}

module.exports = { trajectoryLikelihood,
                   factorial,
                   initPromise,
                   indelTrajectoryDegeneracy,
                   indelTrajectoryLikelihood,
                   testCZL,
                   defaultChopZoneLikelihoodConfig,
                   chopZoneLikelihood,
                   chopZoneLikelihoods }
