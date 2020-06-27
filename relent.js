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

const distrib = files.map ((filename) => {
  let d = fs.readFileSync(filename).toString().split("\n").map ((line) => line.split(" ").map ((x) => parseFloat(x)).filter ((p) => !isNaN(p))).filter ((line) => line.length)
  if (opt.options.normalize) {
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
  let d = 0
  dist1.forEach ((row, i) => {
    row.forEach ((p, j) => {
      if (p > 0)
        d += p * Math.log (p / dist2[i][j])
    })
  })
  return d / Math.log(2)
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
