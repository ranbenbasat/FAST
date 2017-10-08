# FAST
A FAST algorithm for approximate volume (weighted frequency) estimation over data streams

This is an implementation of the FAST algorithm for efficient estimation of flow volumes (weighted frequency estimation)
the algorithm is presented in the paper Fast Flow Volume Estimation (ACM ICDCN 2018) 
by  Ben Basat, Einziger and Friedman

To run, you will need a trace in which each line is a flow identifier, space, and the size of the packet. For example:
0 167
0 1460
4 40

It is built on top of the open source library released by Cormode and Hadjieleftheriou 
(see http://hadjieleftheriou.com/frequent-items/index.html)

Most parameters are tunable from the command line, e.g.:
frequent-items.exe -phi 0.03125 -M 65528 -reps 1 -gamma 0.0078125 -f myTrace.txt

Here, 0.03125 is the error parameter, 65528 is a bound on the packet size, and 1 is the number of times we repeat the time measurement
gamma is a parameter that affects a space-time tradeoff as discussed in the paper

The main() function is inside a file called hh.cc, and FAST is implemented in FAST.*.

Inside hh.cc you'll need to determine which algorithms you wish to run (we recommend running one at a time when measuring speed).
Currently, it is configured for running FAST alone.

Contact us at: sran@cs.technion.ac.il / gil.einziger@nokia.com

