import argparse
import datetime
import random

# Takes in a file that is the output of the ms hudson simulator
def load_haps(filename):
    haps = []
    with open(filename, 'rb') as f:
        for line in f:
            if line[0] == ">":
                pass
            else:
                haps.append(line)
    return haps

# takes in a list of haplotypes and merges the first two into the first genotype etc
def derive_genotypes(haps):
    if len(haps) % 2 != 0:
        raise Exception("Number of haplotypes not divisible by 2!")
    else:
        genotypes = []
        for i in xrange(len(haps) / 2):
            hap1 = haps[i - 1]
            hap2 = haps[i]
            if len(hap1) != len(hap2):
                raise Exception("Length of haplotypes not equal!")
            genotype = []
            for j in xrange(len(hap1)):
                if hap1[j] == hap2[j]:
                    genotype.append(hap1[j])
                else:
                    genotype.append("2")
            genotypes.append(genotype)
    return genotypes

# Takes in a list of genotypes where each list is itself a list of individual characters and outputs 
# the entire block
def print_sequence(geno):
    for gen in geno:
        print ''.join(gen).rstrip()
        
def generateVcfRow(name, hap1, hap2):
    pass

def writeVcf(path):

    rng = random.randint
    print rng(0,1000)

    if True:
        with open(path, 'w') as f:
            f.write("##fileformat=VCFv4.0\n")
            f.write("##filedate=" 
                    + str(datetime
                        .date
                        .today())
                    .replace("-", ""))
            f.write("##reference=simulated data")
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', '-i', help="The file containing the output of a call to ms")
    parser.add_argument('-output', '-o', help="The name of the file to output. Please include the extension.")
    args = parser.parse_args()
    haps = load_haps(args.input)
    writeVcf(args.output)
    # print_sequence(haps)
    # geno = derive_genotypes(haps)
    # print_sequence(geno)
