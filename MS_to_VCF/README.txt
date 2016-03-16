This script takes as input the output of the hudson simulator.
This is generated with a command similar to

ms 4 10 -t 5.0 -s 500

where we make 4 haplotypes that each contain 500 alleles.

Note that it expects at least 4 haplotypes in the file in order to simulate
a family. If more than 4 is provided, it only takes the first 4.
