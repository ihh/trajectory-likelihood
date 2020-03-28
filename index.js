const machineboss = require ('machineboss')

// Algorithm 1 of (Miklos, Lunter & Holmes, 2004)
// exitRates = chi
// transitionRates = r
const trajectoryLikelihood = (exitRates, transitionRates, time) => {
  if (T <= 0)
    throw new Error ("Trajectory must have finite duration")
  if (exitRates.length === 0)
    throw new Error ("There must be at least one state in the trajectory")
  if (exitRates.length + 1 !== transitionRates.length)
    throw new Error ("Mismatch between numbers of exit rates & transition rates")
  if (transitionRates.filter ((r) => r < 0))
    throw new Error ("Transition rates must all be nonnegative")
  if (transitionRates.filter ((r, i) => r <= exitRates[i]))
    throw new Error ("Exit rates must be at least as great as transition rates")
  // switch to the notation from the paper
  let zeta = exitRates,
      r = transitionRates,
      T = time,
      chi = [exitRates[0]],  // unique exit rates
      d = [0],  // exit rate degeneracy
      c = [[1]]   // c[n][k] = coefficient of T^k exp(-chi[n]*T) in likelihood
  zeta.slice(1).forEach ((zeta_i, i) => {
    const j = chi.indexOf (zeta_i)
    if (j < 0) {
      chi.push (zeta_i)
      d.push (0)
    } else
      ++d[j]
    let u = [],  // u[n][k] = new value of c[n][k]
        M = chi.length + 1
    chi.forEach ((chi_n, n) => {
      let u_n = []
      for (let k = 0; k <= d[n]; ++k) {
        let u_nk = 0
        if (chi_n != zeta_i)
          for (let j = k; j <= d[n]; ++j)
            u_nk -= c[n][j] * factorial(j) / (factorial(k) * Math.pow (chi[n] - zeta[i], j - k + 1))
        else if (k == 0) {
          for (let m = 0; m <= M; ++m)
            if (m != n)
              for (let j = k; j <= d[m]; ++j)
                u_nk += c[m][j] * factorial(j) / Math.pow (chi[m] - zeta[i], j+1)
        } else
          u_nk = c[n][k-1] / k
        u_n.push (u_nk)
      }
      u.push (u_n)
    })
    c = u
  })
  let waitProb = 0,
      M = chi.length + 1
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
// The zone lengths are
//  4 -> 8 -> 4 -> 3
const indelTrajectoryLikelihood = async (zoneLengths, params, time) => {
  const { gamma, mu, r } = params
  const N = zoneLengths.length
  if (N === 0)
    throw new Error ("There must be at least one state in the trajectory")
  if (zoneLengths.filter ((l) => l < 0).length)
    throw new Error ("All zone lengths must be nonnegative")
  if (zoneLengths.filter ((l, n) => l == 0 && n < N-1).length)
    throw new Error ("Only the final zone length in the trajectory can be zero")
  if (zoneLengths.filter ((l, n) => n < N-1 && l == zoneLengths[n+1]).length)
    throw new Error ("No two adjacent states in the trajectory can be identical")
  const exitRates = zoneLengths.map ((l) => exitRateForZoneLength(l))
  const transitionRates = zoneLengths.slice(1).map ((l, n) => transitionRateForZoneLengthChange (zoneLengths[n-1], l, params))
  return trajectoryLikelihood (exitRates, transitionRates, time) * indelTrajectoryDegeneracy (zoneLengths)
}

const exitRateForZoneLength = (zoneLength, params) => {
  const { gamma, mu, r } = params
  const lambdaSum = gamma * mu * (1-r) * (1-r) / (1 - gamma*r)
  const muSum = mu * (1-r)
  return zoneLength * (lambdaSum + muSum)
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
          ? insertionRate (destZoneLength - srcZoneLength)
          : deletionRate (destZoneLength - srcZoneLength))
}

const indelTrajectoryDegeneracy = async (zoneLengths) => {
  let eventMachine = []
  for (let n = 0; n < zoneLengths.length - 1; ++n) {
    const deltaLength = zoneLengths[n+1] - zoneLengths[n]
    eventMachine.push (deltaLength > 0
                       ? makeInsertionMachine (deltaLength)
                       : makeDeletionMachine (deltaLength))
  }
  const inputSeq = String.repeat ("A", zoneLengths[0])
  const outputSeq = String.repeat ("X", zoneLengths[zoneLengths.length - 1])
  const mbResult = await machineboss.runWithFiles (eventMachine
                                                   .concat (['--input-chars', inputSeq,
                                                             '--output-chars', outputSeq,
                                                             '--loglike']))
  const mbJson = JSON.parse (mbResult.stdout)
  return Math.exp (mbJson[0][2])
}

const extend = function() {
  let a = arguments[0]
  Array.from(arguments).slice(1).forEach ((b) => Object.keys(b).forEach ((k) => a[k] = b[k]))
  return a
}

const makeTransition = (protoTrans, dest) => {
  return extend ({}, protoTrans, { to: dest })
}

// Ancestral residue = "A"
// Inserted residue = "X"
const makeMachine = (protoTrans, copies) => {
  let state = []
  for (let n = 0; n < copies; ++n)
    state.push (extend ({ id: n,
                          trans: ((n ? [] : [{ to: 0, in: "A", out: "A" },
                                             { to: 0, in: "X", out: "X" }])
                                  .concat (protoTrans.map (pt => makeTransition (pt, n+1)))) }))
  state.push ({ id: copies, trans: [{ to: copies, in: "A", out: "A" },
                                    { to: copies, in: "X", out: "X" }] })
  return { state }
}

const makeInsertionMachine = (k) => {
  return makeMachine ([{ out: "X" }], k)
}

const makeDeletionMachine = (k) => {
  return makeMachine ([{ in: "A" }, { in: "X" }], k)
}

module.exports = { trajectoryLikelihood,
                   indelTrajectoryLikelihood,
                   factorial }
