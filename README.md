# Title

*In the link below, replace analysis-template with your repository name and then delete this text*

Description...

![example event parameter](https://github.com/munch-group/relate1Kgenomes/actions/workflows/quarto-publish.yml/badge.svg?event=push)

echo chrom,pos,p_half_freq,p_two_alleles > ../../results/relate_snps_p_vals.csv ; for CHROM in chr* ; do for POP in $CHROM/* ; do pop=$(basename $POP) ; pop=$(basename $POP) && cut -f '1,34,35' -d ' ' $CHROM/$pop/haplotypes_demog_sele.sele | awk "{print \"$CHROM \" \"$pop \" \$0}" | tail -n +2 | awk '$4 < -2 || $5 < -2' | tr ' ' ',' ; done ; done >> ../../results/relate_snps_p_vals.csv


