#!/usr/bin/env node
// emacs mode -*-JavaScript-*-

const fs = require('fs'),
      getopt = require('node-getopt')

// parse command-line options
const opt = getopt.create([
  ['h' , 'help'             , 'display this help message']
])              // create Getopt instance
      .bindHelp()     // bind option 'help' to default action
      .parseSystem() // parse command line

const files = opt.argv
if (files.length < 2)
  throw new Error ("Please specify at least two distribution files")

const distrib = files.map ((filename) => {
  return fs.readFileSync(filename).toString().split("\n").map ((line) => line.split(" ").map ((x) => parseFloat(x)))
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
  let d = 0
  dist1.forEach ((row, i) => {
    row.forEach ((p, j) => {
      if (p > 0)
        d += p * Math.log (p / dist2[i][j])
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
}

files.forEach ((file, i) => {
  console.log ("S(" + file + ")= " + entropy (distrib[i]))
})

files.forEach ((file1, i) => {
  files.forEach ((file2, j) => {
    if (i != j)
      console.log ("D(" + file1 + "||" + file2 + ")= " + relent (distrib[i], distrib[j]))
  })
})
