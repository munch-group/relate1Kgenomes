# Title

*In the link below, replace analysis-template with your repository name and then delete this text*

Description...

![example event parameter](https://github.com/munch-group/relate1Kgenomes/actions/workflows/quarto-publish.yml/badge.svg?event=push)

for CHROM in steps/relate/* ; do for POP in $CHROM/* ; do chrom=$(basename $CHROM) ; pop=$(basename $POP) ; pop=$(basename $POP) && cut -f '1,34,35' -d ' ' $CHROM/$pop/haplotypes_demog_sele.sele | awk "{print \"$pop $chrom \" \$0}" | tail -n +2 | awk '$4 < -6 || $5 < -6' | tr ' ' ',' ; done ; done > results/all_snps_below_p_6.csv

grep chrX results/all_snps_below_p_6.csv > results/all_chrX_snps_below_p_6.csv 
