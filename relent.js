#!/usr/bin/env node
// emacs mode -*-JavaScript-*-

const fs = require('fs'),
      getopt = require('node-getopt')

// parse command-line options
const opt = getopt.create([
  ['n' , 'normalize'        , 'normalize distributions'],
  ['h' , 'help'             , 'display this help message']
])              // create Getopt instance
      .bindHelp()     // bind option 'help' to default action
      .parseSystem() // parse command line

const files = opt.argv
if (files.length < 2)
  throw new Error ("Please specify at least two distribution files")

const normalize = opt.options.normalize
const distrib = files.map ((filename) => {
  let d = fs.readFileSync(filename).toString().split("\n").map ((line) => line.split(" ").map ((x) => parseFloat(x)).filter ((p) => !isNaN(p))).filter ((line) => line.length)
  if (normalize) {
    const norm = d.reduce ((sum, row) => row.reduce ((sum, p) => sum + p, sum), 0)
    d = d.map ((row) => row.map ((p) => p / norm))
  }
  return d
})

const entropy = (dist) => {
  let s = 0
  dist.forEach ((row) => {
    row.forEach ((p) => {
      if (p > 0)
        s -= p * Math.log (p)
    })
  })
  return s / Math.log(2)
}

const relent = (dist1, dist2) => {
  let d = 0, norm = 0
  dist2.forEach ((row, i) => {
    row.forEach ((p, j) => {
      if (i < dist1.length && j < dist1[i].length)
        norm += dist1[i][j]
    })
  })
  dist1.forEach ((row, i) => {
    row.forEach ((p, j) => {
      if (normalize)
        p /= norm
      const inRange = i < dist2.length && j < dist2[i].length
      if (p > 0 && (inRange || !normalize))
        d += p * Math.log (p / (inRange ? dist2[i][j] : 0))
    })
  })
  return d / Math.log(2)
}

const moments = (dist) => {
  let ei = 0, ed = 0, ei2 = 0, ed2 = 0, eid = 0, norm = 0
  dist.forEach ((row, i) => {
    row.forEach ((p, d) => {
      ei += i*p
      ed += d*p
      ei2 += i*i*p
      ed2 += d*d*p
      eid += i*d*p
      norm += p
    })
  })
  ei /= norm
  ed /= norm
  ei2 /= norm
  ed2 /= norm
  eid /= norm
  return { ei, ed, ei2, ed2, eid,
           vi: ei2 - ei*ei,
           vd: ed2 - ed*ed,
           cid: eid - ei*ed }
}

files.forEach ((file, i) => {
  const m = moments (distrib[i])
  console.log ("Ei " + file + " " + m.ei)
  console.log ("Ed " + file + " " + m.ed)
  console.log ("Vi " + file + " " + m.vi)
  console.log ("Vd " + file + " " + m.vd)
  console.log ("Cid " + file + " " + m.cid)
})

files.forEach ((file, i) => {
  console.log ("S " + file + " " + entropy (distrib[i]))
})

files.forEach ((file1, i) => {
  files.forEach ((file2, j) => {
    if (i != j)
      console.log ("D " + file1 + " " + file2 + " " + relent (distrib[i], distrib[j]))
  })
})
