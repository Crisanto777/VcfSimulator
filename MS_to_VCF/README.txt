This script takes as input the output of the hudson simulator.
This is generated with a command similar to

ms 4 10 -t 5.0 -s 500 > simdata.txt

where we make 4 haplotypes that each contain 500 alleles and 
pipe the output to simdata.txt

Note that it expects at least 4 haplotypes in the file in order to simulate
a family. If more than 4 is provided, it only takes the first 4.


This can be run with the command

python ms_to_vcf.py -i ../data/<filename.txt> -o ../output

mandatory arguments:
	-i --input <(relative) path to .txt containing output of ms simulator>
	-o --output <(relative) path to folder to put output files> 

optional arguments:
	-fathername (String, default=NA0001)
	-mothername (String, default=NA0002)
	-childname (String, default=NA0003)
	-phased (True/False, default=False))
	-blockfile (True/False, default=False)
	-blockdensity (int, default=5)