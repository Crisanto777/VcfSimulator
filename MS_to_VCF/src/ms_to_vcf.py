import argparse
import datetime
import random
from scipy.stats import poisson


def setup_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', '-i', help="The file containing the output of a call to ms that has at least 4 haplotypes")
    parser.add_argument('-output', '-o', help="The name of the folder to place the output.")
    parser.add_argument('-phased', '-ph', default=False, help="Whether you want the vcf to be phased or unphased")
    parser.add_argument('-blockfile', '-b', default=False, help="Whether you want the program to spoof a hapcompass blockfile")
    parser.add_argument('-blockdensity', '-bd', default=5, help="Expected number of blocks per 100 SNPs")
    parser.add_argument('-fathername', '-f', default="NA0001", help="name of father files, without file extension")
    parser.add_argument('-mothername', '-m', default="NA0002", help="name of mother files, without file extension")
    parser.add_argument('-childname', '-c', default="NA0003", help="name of child files, without file extension")
    
    return parser


# Takes in a file that is the output of the ms hudson simulator
def load_haps(filename):
    haps = []
    with open(filename, 'rb') as f:
        # the start of the file contains information unrelated
        # to the haplotypes
        start_processing = False
        for line in f:
            if start_processing:
                haps.append(line.strip())
            else:
                # the line that specifies positions precedes the
                # start of the list of haplotypes
                if len(line) > 0:
                    try:
                        if line.split()[0] == "positions:":
                            start_processing = True
                    except:
                        # empty lines will throw an exception here
                        # we ignore this
                        pass
    return haps

# Takes in a list of sequences where each list is itself a list of individual characters and outputs 
# the entire block
def print_sequence(geno):
    for gen in geno:
        print ''.join(gen).rstrip()
        

def writeVcf(path, genome_name, snp_identity_hash, hapseq1, hapseq2, is_phased):

    # genome_name is something of the form NA0001 etc
    # can set this to false if need to debug
    if True:
        newPath = path + "/" + genome_name + ".vcf"
        print "creating file:", newPath
        with open(newPath, 'w') as f:
            f.write("##fileformat=VCFv4.1\n")
            f.write("##filedate=" 
                    + str(datetime
                        .date
                        .today())
                    .replace("-", "") + "\n")
            f.write("##reference=simulated data\n")
            f.write("##phasing=phased\n")
            f.write("##FILTER=<ID=q10,Description=\"Quality below 10\">\n")
            f.write("##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n")
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            f.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
            f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
            f.write("##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n")
            f.write("#CHROM POS\tID\tREF\tALT\tQUAL FILTER INFO\tFORMAT\t" + genome_name + "\n")

            # header is written. Now we generate the data


            for i in xrange(len(hapseq1)):

                chrom_num = snp_identity_hash[i]["chrom_num"]
                snp_pos = snp_identity_hash[i]["snp_pos"]
                snp_id = snp_identity_hash[i]["snp_id"]
                snp_var_0 = snp_identity_hash[i]["variant"][0]
                snp_var_1 = snp_identity_hash[i]["variant"][1]
                hap1 = hapseq1[i]
                hap2 = hapseq2[i]
                vcf_row = generateVcfRow(chrom_num, snp_pos, snp_id, snp_var_0, snp_var_1, hap1, hap2, is_phased)
                f.write(vcf_row)


def generateVcfRow(chrom_num, snp_pos, input_id, ref, alt, hap1, hap2, is_phased):
    # for more info look at http://www.1000genomes.org/wiki/Analysis/vcf4.0

    chrom_num = str(chrom_num)
    snp_pos = str(snp_pos)
    my_id = str(input_id)
    ref = str(ref)
    qual = "99"
    my_filter = "PASS"
    info = "."
    my_format = "GT:GQ:DP:HQ"

    if is_phased:
        alleles = str(hap1) + "|" + str(hap2) 
    else:
        alleles = str(hap1) + "/" + str(hap2)

    formatted_info = str(alleles + ":"
                      + "48" + ":"
                      + "1" + ":" 
                      + "51,51")

    row_data = [chrom_num, snp_pos, my_id, ref, alt, qual, my_filter, info, my_format, formatted_info]
    return "\t".join(row_data) + "\n"


def generate_snp_identity_hash_entry(chrom_num, rand_num):
    snp_info = {}
    snp_info["chrom_num"] = str(chrom_num)
    snp_info["snp_pos"] = rand_num
    snp_info["snp_id"] = "rs" + str(rand_num * 7)
    snp_info["variant"] = simulate_variants()
    
    return snp_info
    
def simulate_variants():
    option_list = ["A", "C", "T", "G"]

    del option_list[random.randint(0, 3)]
    del option_list[random.randint(0, 2)]

    return option_list

def write_blockfile(path, file_name, hap1, hap2, snp_identity_hash, block_density):
    #choose a random number of blocks, distributed poisson, dependent on the length of the haplotype
    #uses block_density (expected number of blocks per 100 SNPs)
    #this number is fudged - i am unsure about the actual mean or distribution
    num_blocks = poisson.rvs((float(len(hap1))/100)*block_density)
    
    #for each block,
    #    choose a position in the genome, distributed at random
    positions = []
    for i in xrange(int(num_blocks)):
        positions.append(random.randint(0, len(hap1)-5))
        
    #sort the positions, since the output file has them in increasing order
    position_set = set(positions)
    positions = list(position_set)
    #eliminate duplicate positions
    positions = sorted(positions)
    
    #for each position:
    #    choose a random number of snps,
    #       distributed poisson(1) + 2 (since min block length should be 2) 
    #    this number is fudged - i am unsure about the actual mean or distribution
    
    newPath = path + "/" + file_name + ".txt"
    print "creating file:", newPath
    with open(newPath, 'w') as f:
        for position in positions:
            block_len = poisson.rvs(1) + 2

            first_block_line = []
            first_block_line.append("BLOCK")
            first_block_line.append(str(snp_identity_hash[position]['snp_pos']))
            #can't phase off the end of the chromosome
            if ((position + block_len) >= len(hap1)):
                block_len = len(hap1)-position-2
            first_block_line.append(str(snp_identity_hash[position+block_len-1]['snp_pos']))
            first_block_line.append("*") #start position in hapcompass out vcf
            first_block_line.append("*") #end position in hapcompass out vcf
            first_block_line.append("*") #unclear 
            first_block_line.append(str(snp_identity_hash[position]['chrom_num']))

            f.write("\t".join(first_block_line) + "\n")

            for block_snp in xrange(block_len):
                snp_num = position+block_snp
                line_output = []
                line_output.append("VAR_POS_"+str(snp_identity_hash[snp_num]["snp_pos"]))
                line_output.append(str(snp_identity_hash[snp_num]["snp_pos"]))
                line_output.append("*")
                line_output.append(str(hap1[snp_num]))
                line_output.append(str(hap2[snp_num]))
                f.write('\t'.join(line_output) + "\n")
            f.write('\n')

def simulate_family(args, haps, num_de_novo_mutations=0):
    father_file_name = args.fathername
    mother_file_name = args.mothername
    child_file_name = args.childname
    
    output = args.output
    is_phased = (args.phased == "True")
    make_blockfiles = (args.blockfile == "True")

    # generate and keep track of snp positions, id's, etc for consistency across all vcf's
    #    this is the snp_identity_hash
    rng = random.randint
    curr_snp_num = rng(0, 1000)
    snp_identity_hash = {}

    chromosome_num = 1
    for i in xrange(len(haps[0])):

        # if you want to switch chromosomes, you can handle it here
        # ie if i > thresh chrom = chrom + 1
        snp_identity_hash[i] = generate_snp_identity_hash_entry(chromosome_num, curr_snp_num)
        curr_snp_num = curr_snp_num + rng(0, 1000)
        
        
    # this takes the first four haplotypes and assumes that 1 and 2 are father, and 3 and 4 are mother

    # first, arbitrarily assume the first two are the father
    writeVcf(output, father_file_name, snp_identity_hash, haps[0], haps[1], is_phased)

    # then, assume the next two are the other
    # generate the mother using the same id's generated by the father
    writeVcf(output, mother_file_name, snp_identity_hash, haps[2], haps[3], is_phased)

    # it then simulates a child by choosing one haplotype from the father and one from the mother
    c1_index = rng(0, 1) 
    c2_index = rng(2, 3) 
    childhaps = [list(haps[c1_index]), list(haps[c2_index])]

    # record events in a file so that we can verify the accuracy of any inference methods
    child_changes = ("father index: " + str(c1_index) + "\n"
                     + "mother index: " + str(c2_index) + "\n")

    with open((output + "/child_info.txt"), 'w') as f:
    
        # This code adds biological events to the child
        # simulate de novo mutations
        de_novo_muts = set()
        for i in xrange(num_de_novo_mutations):
            hap = rng(0, 1)
            pos = rng(0, len(childhaps[0])-1)
            while (hap, pos) in de_novo_muts:
                hap = rng(0, 1)
                pos = rng(0, len(childhaps[0]))
            if childhaps[hap][pos] == "0":
                childhaps[hap][pos] = "1"
            else:
                childhaps[hap][pos] = "0" 
            f.write("de novo mutation at\t" + str(hap) + "," + 
                    " POS = " + str(snp_identity_hash[pos]["snp_pos"]) + 
                    " ID = " + str(snp_identity_hash[pos]["snp_id"]) + "\n");
            de_novo_muts.add((hap, pos))

        # simulate gene conversion
        # simulate crossover event

        f.write(child_changes)
        writeVcf(output, child_file_name, snp_identity_hash, childhaps[0], childhaps[1], is_phased)
        
        if make_blockfiles:
            """
            for father, mother, and child
                make block file
            """
            write_blockfile(output, father_file_name, haps[0], haps[1], snp_identity_hash, float(args.blockdensity))
            write_blockfile(output, mother_file_name, haps[2], haps[3], snp_identity_hash, float(args.blockdensity))
            write_blockfile(output, child_file_name, childhaps[0], childhaps[1], snp_identity_hash, float(args.blockdensity))
    
if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args()
    haps = load_haps(args.input)
    num_de_novo_mutations = 2
    simulate_family(args, haps, num_de_novo_mutations)    
