import pandas as pd
import sys
import gzip

_, vcf_file_name, chrom, *args = sys.argv

first_non_par_pos, last_non_par_pos = None, None
if args:
    assert len(args) == 2, args
    first_non_par_pos, last_non_par_pos = args
    first_non_par_pos = int(first_non_par_pos)
    last_non_par_pos = int(last_non_par_pos)

males = []
females = []

het_counts = dict()



if chrom == 'chrX':
    with gzip.open(vcf_file_name, 'rt') as f:

        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *individuals = line.split()
                    break
        for line in f:

            if not line.startswith('#'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()

                for i, call in enumerate(calls):   
                    vars = call.split('|')            
                    het_counts[individuals[i]] += vars[0] != vars[-1] # also false if only one var

        l = []
        for i in range(len(individuals)):
            l.append(het_counts[individuals[i]])
        import numpy as np
        a = np.array(l)
        cutoff = np.mean(a[a!=0]) / 100
        with open('sexes.txt', 'w') as sex_file:                        
            for i, sample in enumerate(individuals):

                if het_counts[sample] > cutoff:
                    females.append(sample)
                    print(sample, 'F', file=sex_file)
                else:
                    males.append(individuals[i])
                    print(individuals[i], 'M', file=sex_file)

        assert len(individuals) == len(males) + len(females)

        males = set(males)
        females = set(females)

    with gzip.open(vcf_file_name, 'rt') as f:
        
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *individuals = line.split()
                    ids = []
                    for sample in individuals:
                        if chrom == 'chrX' and sample in males:
                            ids.append(sample)
                        else:
                            ids.append(sample + "_1")
                            ids.append(sample + "_2")

                    start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
                    start_line.extend(ids)
                    print('\t'.join(start_line))
                else:
                    print(line, end='')
            else:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()

                haploid_calls = []
                for i, call in enumerate(calls):
                    for x in call.split('|'):
                        haploid_calls.append(x)
                        if CHROM == 'chrX' and individuals[i] in males:
                            break

                assert len(haploid_calls) == len(ids), (len(haploid_calls), len(ids), len(males), len(females))

                start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
                start_line.extend(haploid_calls)
                print('\t'.join(start_line))

else:
    with gzip.open(vcf_file_name, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *individuals = line.split()
                    print('\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT] + [f'{sample}_{i}' for sample in individuals for i in [1, 2]]))
                else:
                    print(line, end='')
            else:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
                print('\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT] + [x for call in calls for x in call.split('|')]))
