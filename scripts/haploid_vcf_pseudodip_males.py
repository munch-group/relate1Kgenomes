import pandas as pd
import sys
import gzip

_, vcf_file_name, first_non_par_pos, last_non_par_pos = sys.argv

first_non_par_pos = int(first_non_par_pos)
last_non_par_pos = int(last_non_par_pos)

males = []
females = []

het_counts = dict()

print('sexing samples', file=sys.stderr)
with gzip.open(vcf_file_name, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
        else:
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            if any('|' not in call for call in calls):
                first_non_par_pos = int(POS)
                print('first non-par pos', first_non_par_pos, file=sys.stderr)
                with open('sexes.txt', 'w') as sex_file:
                    for i, call in enumerate(calls):
                        if '|' in call:
                            females.append(all_samples[i])
                            print(all_samples[i], 'F', file=sex_file)
                        else:
                            males.append(all_samples[i])
                            print(all_samples[i], 'M', file=sex_file)
                assert len(all_samples) == len(males) + len(females)
                sexing_done = True
                break

            if first_non_par_pos <= int(POS) <= last_non_par_pos:
                for i, call in enumerate(calls):               
                    het_counts[all_samples[i]] += call[0] != call[2]

if not sexing_done:
    l = []
    for i in range(len(all_samples)):
        l.append(het_counts[all_samples[i]])
    import numpy as np
    a = np.array(l)
    cutoff = np.mean(a[a!=0]) / 100
    with open('sexes.txt', 'w') as sex_file:                        
        for i, sample in enumerate(all_samples):

            if het_counts[sample] > cutoff:
                females.append(sample)
                print(sample, 'F', file=sex_file)
            else:
                males.append(all_samples[i])
                print(all_samples[i], 'M', file=sex_file)

assert len(all_samples) == len(males) + len(females)


print('writing haploids', file=sys.stderr)
with gzip.open(vcf_file_name, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
                ids = []
                for sample in all_samples:
                    if CHROM == '#chrX' and sample in males:
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
                    if CHROM == '#chrX' and all_samples[i] in males:
                        break

            start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
            start_line.extend(haploid_calls)
            print('\t'.join(start_line))
