# trajectory-likelihood

Calculate indel probabilities using the trajectory enumeration method of Miklos, Lunter &amp; Holmes (2004),
[A "Long Indel" model for evolutionary sequence alignment](https://www.ncbi.nlm.nih.gov/pubmed/14694074).

Currently this code only calculates the gap likelihoods for internal chop zones
(i.e. infinitely long sequences).

To find the number of trajectories consistent with a given chop zone, this code uses a finite state machine approach.
This is not optimized for short trajectories, and so is (even) slower than the original method of Miklos _et al_.

## Installation

### C++ version (preferred)

First download and build [Machine Boss](https://github.com/evoldoers/machineboss).
Install the Machine Boss library globally with `make install-lib`.

Then run `make`

And finally, to calculate the chop zone likelihoods

~~~~
./trajeclike --gamma .99 --mu .1 -r .2 --maxlen 6 --maxevents 3 --time 1 -v1
~~~~


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

