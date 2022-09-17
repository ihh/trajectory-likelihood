# trajectory-likelihood

This is a reference implementation of several methods to calculate alignment gap probabilities using several different approximations to
continuous-time Markov chains models of sequence evolution with insertions and deletions of geometrically distributed lengths
(the "general geometric indel" model of [De Maio, 2020](https://pubmed.ncbi.nlm.nih.gov/32653921/);
a special case of the "long indel" model of [Miklós, Lunter &amp; Holmes, 2004](https://www.ncbi.nlm.nih.gov/pubmed/14694074)).

The methods implemented by this code include
- trajectory enumeration, as described by Miklós, Lunter &amp; Holmes in the paper [A "Long Indel" model for evolutionary sequence alignment](https://www.ncbi.nlm.nih.gov/pubmed/14694074) (2004);
- manual moment-matching, as described by De Maio in the paper [The Cumulative Indel Model: Fast and Accurate Statistical Evolutionary Alignment](https://pubmed.ncbi.nlm.nih.gov/32653921/) (2020);
- renormalized moment-matching, as described by Holmes in the paper [A Model of Indel Evolution by Finite-State, Continuous-Time Machines
](https://pubmed.ncbi.nlm.nih.gov/33020189/) (2020).

## Notes on the implementation of trajectory enumeration from MLH 2004

The code includes an implementation of Algorithm 1 from the MLH 2004 paper, which computes the likelihood of a state trajectory in an arbitrary continuous-time discrete-state Markov chain,
with event times integrated out.

For clarity, this code currently makes the following simplifications relative to the full "Long Indel" model described in the paper:
- Currently it only calculates the gap likelihoods for internal chop zones (i.e. not the regions at the ends of the sequence, nor the probability of deleting/inserting the entire sequence). This is equivalent to considering only infinitely long sequences.
- Currently it does not use a mixture of geometric distributions for the indel lengths, just a single geometric distribution.

Generalizing the code to cover these cases would be relatively straightforward, and only involve changing the calculation of zone exit and transition rates.
The "work" of the implementation is largely in counting trajectories and calculating their probabilities.

## Installation

The [Boost](https://www.boost.org/) library is required for the command-line tool (though not if you just want to link to the C++ code).

Type `make` to compile the command-line tool.

## Running

~~~~
./trajeclike --mu .049 --y .543 --time 1 --maxlen 20 --maxevents 3 -v1
~~~~

The values of the mu and y parameters here are (from the MLH 2004 paper) the maximum likelihood fit of a single-component geometric indel distribution to protein alignments in the [HOMSTRAD](https://www.ncbi.nlm.nih.gov/pubmed/9828015) database.

To compare to simulation results, run with the `--simulate` switch
(note that this simulates directly from the underlying continuous-time Markov chain)

~~~~
./trajeclike --mu .049 --y .543 --time 1 --maxlen 20 --simulate -v1
~~~~

For more verbose logging, increase the verbosity (`-v`).
For example, `-v1` will display a legend for the output probability matrix. Higher verbosity levels will display more information, such as progress messages.

For help, use `-h`.

## Analysis

If you turn on verbose logging (`-v1`), the program will report summary statistics including the moments of the gap length distribution.

The node script `relent.js` can be used to compute the relative entropy between gap length distributions computed using different methods.
