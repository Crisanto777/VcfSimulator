import argparse
import datetime
import random
import derive_genotypes


def setup_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', '-i', help="The file containing the output of a call to ms that has at least 4 haplotypes")
    parser.add_argument('-output', '-o', help="The name of the folder to place the output.")
    
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
        

def writeVcf(path, genome_name, snp_identity_hash, hapseq1, hapseq2):

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
                vcf_row = generateVcfRow(chrom_num, snp_pos, snp_id, snp_var_0, snp_var_1, hap1, hap2)
                f.write(vcf_row)


def generateVcfRow(chrom_num, snp_pos, id, ref, alt, hap1, hap2):
    # for more info look at http://www.1000genomes.org/wiki/Analysis/vcf4.0

    chrom_num = str(chrom_num)
    snp_pos = str(snp_pos)
    my_id = str(id)
    ref = str(ref)
    qual = "99"
    my_filter = "PASS"
    info = "."
    my_format = "GT:GQ:DP:HQ"
    alleles = str(hap1) + "|" + str(hap2) 
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



def simulate_family(output, haps, num_de_novo_mutations=0):
        
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
    writeVcf(output, "NA0001", snp_identity_hash, haps[0], haps[1])

    # then, assume the next two are the other
    # generate the mother using the same id's generated by the father
    writeVcf(output, "NA0002", snp_identity_hash, haps[2], haps[3])

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
            pos = rng(0, len(childhaps[0]))
            while (hap, pos) in de_novo_muts:
                hap = rng(0, 1)
                pos = rng(0, len(childhaps[0]))
            if childhaps[hap][i] == "0":
                childhaps[hap][i] = "1"
            else:
                childhaps[hap][i] = "0" 
            f.write("de novo mutation at\t" + str(hap) + "," + str(pos) + "\n")
            de_novo_muts.add((hap, pos))

        # simulate gene conversion
        # simulate crossover event

        f.write(child_changes)
        writeVcf(output, "NA0003", snp_identity_hash, childhaps[0], childhaps[1])

if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args()
    haps = load_haps(args.input)
    num_de_novo_mutations = 2
    simulate_family(args.output, haps, num_de_novo_mutations)    
