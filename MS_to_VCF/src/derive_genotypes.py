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

