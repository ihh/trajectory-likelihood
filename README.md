# trajectory-likelihood

Calculate indel probabilities using the trajectory enumeration method of Miklos, Lunter &amp; Holmes (2004),
[A "Long Indel" model for evolutionary sequence alignment](https://www.ncbi.nlm.nih.gov/pubmed/14694074).

For clarity, this code currently makes the following simplifications relative to the full model described in the paper:
- Currently it only calculates the gap likelihoods for internal chop zones (i.e. not the regions at the ends of the sequence, nor the probability of deleting/inserting tbhe whole entire sequence). This is equivalent to considering only infinitely long sequences.
- Currently it does not use a mixture of geometric distributions for the indel lengths, just a single geometric distribution.

Both these modifications are relatively straightforward to implement, as they only involve changing the calculation of zone exit probabilities and transition rates.
The "work" of the implementation is largely in counting trajectories and calculating their probabilities.

## Installation

The [Boost](https://www.boost.org/) library is required for the command-line tool (though not if you just want to link to the C++ code).

~~~~
make
./trajeclike --gamma .99 --mu .049 -r .543 --maxlen 20 --maxevents 3 --time 1 -v1
~~~~

The values of the mu and r parameters here are (from the MLH 2004 paper) the maximum likelihood fit of a single-component geometric indel distribution to protein alignments in the [HOMSTRAD](https://www.ncbi.nlm.nih.gov/pubmed/9828015) database.


### JavaScript version

~~~~
npm install
~~~~

...then (in node)...

~~~~
var tl = require("./index")
console.log (tl.chopZoneLikelihoods ({ gamma: .99, r: .1, mu: .1 },
	    			     1,  // time
				     { maxLen: 5 })
~~~~

