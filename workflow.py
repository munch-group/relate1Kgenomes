

# %% [markdown]
# ---
# title: GWF workflow
# execute:
#   eval: false
# ---

# %%

import os
os.environ['NUMEXPR_MAX_THREADS'] = '16'

from gwf import Workflow
import re
from collections import defaultdict
from pathlib import Path
import pandas as pd
from gwf import Workflow, AnonymousTarget
import os, re
from collections import defaultdict
from os.path import abspath, dirname, basename, relpath
from os.path import join as joinpath

def stripsuf(p, strip=-1):
    return p.rsplit('.', strip)[0]

def collect(target, label):
    """
    function to collect fields from a list of dictionaries
    """
    res = defaultdict(list)    
    for out in target.outputs:
        for key, val in out.items():
            if key == label:
                res[key].append(val)
    return dict(res)


def modpath(p, parent=None, base=None, suffix=None):
    """
    function that modifies file path
    """
    
    os.path.basename(p).split('.', 1)

    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = joinpath(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


def combine(*args, only=None):
    """
    function to combine 2 target outputs as an input
    """
    assert all(len(args[0]) == len(args[i]) for i in range(len(args)))
    combined = []
    for j in range(len(args[0])):
        output_group = {}
        for i in range(len(args)):
            if only:
                output_group.update({k: v for k, v in args[j].items() if k in only})
            else:
                output_group.update(args[i][j])
        combined.append(output_group)
    return combined

def download_data(config, chrom):
    """
    download vcf files for a given chromosome
    """
    inputs = []
    outputs = {f'sample_vcf': f"steps/data/{basename(config['sample_vcf'][chrom])}",
               f'sample_vcf_index': f"steps/data/{basename(config['sample_vcf_index'][chrom])}",
               f'mask': f"steps/data/{basename(config['mask'][chrom]).replace('.gz', '')}", # remove gz suffix because it is unpacked
#               f'ancestral_fa': f'steps/data/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom[3:]}.fa'

    }

    options = {'memory': '8g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p steps/data
    wget --directory-prefix steps/data {config['sample_vcf'][chrom]}
    wget --directory-prefix steps/data {config['sample_vcf_index'][chrom]}
    wget --directory-prefix steps/data {config['mask'][chrom]}
    gzip -d steps/data/{basename(config['mask'][chrom])}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def download_meta(config, chromosomes):
    inputs = []
    outputs = {'1000G_2504_seq_index': f"steps/data/{basename(config['1000G_2504_seq_index'])}",
               '1000G_698_seq_index': f"steps/data/{basename(config['1000G_698_seq_index'])}",
    }
    for chrom in chromosomes:
        outputs[f'ancestral_fa_{chrom}'] = f'steps/data/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom[3:]}.fa'

    options = {'memory': '8g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p steps/data
    wget --directory-prefix steps/data {config['ancestral_fa']}
    wget --directory-prefix steps/data {config['1000G_2504_seq_index']}
    wget --directory-prefix steps/data {config['1000G_698_seq_index']}
    cd steps/data/
    tar xfvz {basename(config['ancestral_fa'])}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def download_data(config, chromosomes):

#     inputs = []
#     outputs = {f'sample_vcf': f"steps/data/1000G_2504_seq_chr{chrom}.vcf.gz",
#                f'sample_vcf_index': f"steps/data/1000G_2504_seq_chr{chrom}.vcf.gz.tbi",
#                f'mask': f"steps/data/{basename(config['mask'][chrom]).replace('.gz', '')}", # remove gz suffix because it is unpacked
#     }
#     for chrom in chromosomes:
#         outputs['ancestral_vcf'] = f"steps/data/{basename(config['ancestral_vcf'])}"
#         outputs['1000G_2504_seq_index'] = f"steps/data/{basename(config['1000G_2504_seq_index'])}"
#         outputs['1000G_698_seq_index'] = f"steps/data/{basename(config['1000G_698_seq_index'])}"
#         outputs['ancestral_fa'] = f'steps/data/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom[3:]}.fa'

#     options = {'memory': '8g', 'walltime': '10:00:00'}
#     spec = f'''
#     mkdir -p steps/data
#     '''
#     for chrom in chromosomes:
#         spec += f'''  
#         wget --directory-prefix steps/data {config['sample_vcf'][chrom]}
#         wget --directory-prefix steps/data {config['sample_vcf_index'][chrom]}
#         wget --directory-prefix steps/data {config['mask'][chrom]}
#         '''
#     spec += f'''
#     wget --directory-prefix steps/data {config['ancestral_vcf']}
#     wget --directory-prefix steps/data {config['1000G_2504_seq_index']}
#     wget --directory-prefix steps/data {config['1000G_698_seq_index']}
#     cd steps/data/
#     tar xfvz {basename(config['ancestral_vcf'])}
#     '''
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def download_data(config, chrom):
#     inputs = []
#     outputs = {
#                'sample_vcf': f"steps/data/{basename(config['sample_vcf'][chrom])}",
#                'sample_vcf_index': f"steps/data/{basename(config['sample_vcf_index'][chrom])}",
#                '1000G_2504_seq_index': f"steps/data/{basename(config['1000G_2504_seq_index'][chrom])}",
#                '1000G_698_seq_index': f"steps/data/{basename(config['1000G_698_seq_index'][chrom])}",
#                'mask': f"steps/data/{basename(config['mask']).replace('.gz', '')}", # remove gz suffix because it is unpacked
#                'ancestral_fa': 'steps/data/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_X.fa',
#     }
#     options = {'memory': '8g', 'walltime': '10:00:00'}
#     spec = f'''
#     mkdir -p steps/data
#     wget --directory-prefix steps/data {config['sample_vcf'][chrom]}
#     wget --directory-prefix steps/data {config['sample_vcf_index'][chrom]}
#     wget --directory-prefix steps/data {config['1000G_2504_seq_index'][chrom]}
#     wget --directory-prefix steps/data {config['1000G_698_seq_index'][chrom]}

#     wget --directory-prefix steps/data {config['mask'][chrom]}
#     cd steps/data/
#     gzip -d {basename(config['mask'][chrom])}
#     tar xfvz {basename(config['ancestral_vcf'])}
#     '''
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map, chrom):
    """
    map of recombination rate across the X chromosome made by DECODE genetics
    """
    inputs = [decode_hg38_sexavg_per_gen]
    outputs = [genetic_map]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    grep = rf'{chrom}\s'
    spec = f'''
    mkdir -p {dirname(genetic_map)}
    cat {decode_hg38_sexavg_per_gen} | tail -n +2 | grep "{grep}" | cut -f 2,4,5 | (echo pos COMBINED_rate Genetic_Map ; cat - ; ) > {genetic_map}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def female_haploid(chrX_filtered_eagle2_phased, chrom):
    """
    turn diploid females (XX) into two individual haplotypes (haploid individuals) like males
    """
    phased_haplotypes = f'steps/relate/{chrom}/haplotypes.vcf.gz'
    inputs = [chrX_filtered_eagle2_phased]
    outputs = {'haplotypes': phased_haplotypes}
    options = {'memory': '10g', 'walltime': '01:20:00'}
    spec = f'''
    python scripts/haploid_vcf.py {chrX_filtered_eagle2_phased} | gzip > {phased_haplotypes}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def haplotype_id(phased_haplotypes, chrom):
    """
    construct files with haplotype IDs
    """
    phased_haplotypes_id = f'steps/relate/{chrom}/haplotypes_ids_{chrom}.txt'
    inputs = [phased_haplotypes]
    outputs = {'ids': phased_haplotypes_id}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    # conda install -c bioconda bcftools
    # conda install openssl   ## to install libcrypto.so.1.0.0 library
    bcftools query -l {phased_haplotypes} > {phased_haplotypes_id}
    sleep 5
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def all_pop_labels(phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index, chrom):
    """
    construct pops labels mapping each haplotype to a pop
    (group haplotypes according to the pop to which the individuals carrying those haplotypes belong)
    """
    phased_haplotypes_poplabels = f'steps/relate/{chrom}/poplabels.txt'
    inputs = [phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index]
    outputs = {'poplabels': phased_haplotypes_poplabels}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    python scripts/make_poplabels.py {phased_haplotypes_id} {high_coverage_seq_index} {related_high_coverage_seq_index} > {phased_haplotypes_poplabels} 
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def convert_vcf(phased_haplotypes, phased_haplotypes_poplabels, chrom):
    """
    Define the function to convert VCF to haps/sample format
    """
    phased_haplotypes_haps = modpath(phased_haplotypes, suffix=f'_{chrom}.haps')
    phased_haplotypes_sample = modpath(phased_haplotypes, suffix=f'_{chrom}.sample')
    inputs = [phased_haplotypes_poplabels, phased_haplotypes]
    outputs = {'haps': phased_haplotypes_haps+'.gz', 'sample': phased_haplotypes_sample+'.gz'}
    options = {'memory': '10g', 'walltime': '01:00:00'}
    spec = f'''
    {config['relate_dist_dir']}/bin/RelateFileFormats --mode ConvertFromVcf --haps {phased_haplotypes_haps} --sample {phased_haplotypes_sample} -i {phased_haplotypes.replace('.vcf.gz', '')} --poplabels {phased_haplotypes_poplabels}
    sleep 20
    gzip --force {phased_haplotypes_haps}
    gzip --force {phased_haplotypes_sample}
    '''
    # Returning outputs as well
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


## start with specific pop ##

def exclude_related(path, pop, chrom):
    """
    exclude related individuals to avoid biases arising from shared genetic material
    """
    output_dir = f'steps/relate/{chrom}/{pop}/excluded'
    output_path = modpath(path, parent=output_dir, suffix=f'_related.txt')
    inputs = {'path' : path}
    outputs = {'path' : output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v '#' {path} | cut -f 10 > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def ids_other_pop(path, pop, chrom):
    """
    find IDs of haplotypes from all other pops so we can exclude them
    """
    output_dir = f'steps/relate/{chrom}/{pop}/excluded'
    output_path = modpath(path, parent=output_dir, suffix='_non_pop.txt')
    inputs = {'path' : path}
    outputs = {'path' : output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v {pop} {path} | cut -f 1 -d ' ' > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def combine_files(path, pop=None, related=None, chrom=None):
    """
    combine excluded files: both related and non pop individuals
    """
    # output_dir = modpath(output_path, base='', suffix='')
    output_dir = f'steps/relate/{chrom}/{pop}/combined'
    output_path = modpath(path, parent=output_dir, base='', suffix=f'excluded_combined_{chrom}.txt')
    inputs = {'path': path, 'related': related}
    outputs = {'path': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    cat {path} {related} | sort | uniq > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def excluded_list(path, pop=None, haplotype_id=None, chrom=None):
    """
    construct a list of excluded individuals
    """
    # output_dir = modpath(output_path, base='', suffix='')
    output_dir = f'steps/relate/{chrom}/{pop}/combined'
    output_path = modpath(path, parent=output_dir, base='', suffix=f'excluded_list.txt')
    inputs = {'path': path, 'haplotype_id': haplotype_id}
    outputs = {'exclude_list': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -f {path} {haplotype_id} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def pop_labels(exclude_list, pop=None, poplabels=None, chrom=None):
    """
    construct a list of only individuals from the pop of interest
    """
    output_dir = f'steps/relate/{chrom}/{pop}/included'
    output_path = joinpath(output_dir, f'included_pop_labels.txt')
    inputs = {'exclude_list': exclude_list, 'poplabels': poplabels}
    outputs = {'poplabels': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v -f {exclude_list} {poplabels} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def prepare_files(exclude_list, pop=None, haps=None, sample=None, ancestor=None, mask=None, poplabels=None, chrom=None):
    """
    prepare input files for RELATE
    """
    output_dir = f'steps/relate/{chrom}/{pop}'
    inputs = {'haps': haps, 
              'sample': sample, 
              'ancestor': ancestor, 
              'mask':mask, 
              'poplabels':poplabels, 
              'exclude_list':exclude_list}
    output_path = joinpath(output_dir, 'haplotypes')
    outputs = {'haps': f'{output_path}.haps.gz', 
               'sample': f'{output_path}.sample.gz', 
               'dist': f'{output_path}.dist.gz', 
               'poplabels': f'{output_path}.poplabels',
               'annot': f'{output_path}.annot'} 
    options = {'memory': '20g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    {config['relate_dist_dir']}/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps {haps} --sample {sample} --ancestor {ancestor} --mask {mask} --remove_ids {exclude_list} --poplabels {poplabels} -o {output_path}
    sleep 20
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# compute sfs to make sure singletons are not missing (sanity check)
# zcat 1000g_LWK_phased_haplotypes.haps.gz | cut -d ' ' -f 4- | tr -d -c '1\n' | awk '{ print length; }' | sort -n | uniq -c

def relate(genetic_map, pop=None, sample=None, haps=None, annot=None, dist=None, chrom=None):
    """
    run the inference of tree sequences using RELATE
    """
    file_name_output = abspath(f'steps/relate/{chrom}/{pop}/haplotypes')
    inputs = {'sample': sample, 
              'haps': haps, 
              'annot': annot, 
              'dist': dist,
              'map': genetic_map,
              }
    outputs = {'anc': file_name_output + '.anc', 
               'mut': file_name_output + '.mut'
               }
    options = {'memory': '24g', 'walltime': '10:00:00'}
    # program creates a temporary folder for temporary files and if it already exists relate won't run
    spec= f'''
    mkdir -p {dirname(file_name_output)}
    cd {dirname(file_name_output)}
    rm -rf {basename(file_name_output)}
    {config['relate_dist_dir']}/bin/Relate --mode All -m 1.25e-8 -N 20000 \
        --sample {relpath(sample, dirname(file_name_output))} \
        --haps {relpath(haps, dirname(file_name_output))} \
        --map {relpath(genetic_map, dirname(file_name_output))} \
        --annot {relpath(annot, dirname(file_name_output))} \
        --dist {relpath(dist, dirname(file_name_output))} \
        --memory 20 -o {basename(file_name_output)}
    sleep 90
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def estimate_pop_sizes(anc, mut=None, poplabels=None):
    """
    estimate historical pop size trajectory from initially inferred tree sequences
    setting --threshold 0. This is so that the branch lengths in all trees are updated for the estimated pop size history. 
    inputs: inferred .anc/.mut files and a .poplabels file
    outputs: two versions of coalescence rates/pop sizes are outputted
    .coal --> contains coalescence rates and cross-coalescence rates, treating all samples as one pop
    *.pairwise.coal/.bin --> coalescence rate file and corresponding binary file containing coalescence rates between pairs of samples
    """
    input_base_path = stripsuf(anc)
    output_base_path = joinpath(dirname(anc), basename(input_base_path) + '_demog')
    inputs = {'anc': anc, 
              'mut': mut, 
              'poplabels_size': poplabels}
    outputs = {'anc': output_base_path + '.anc.gz', 
               'mut': output_base_path + '.mut.gz', 
               'coal': output_base_path + '.coal', 
               'pairwise.coal': output_base_path + '.pairwise.coal', 
               'pairwise.bin': output_base_path + '.pairwise.bin'}
    options = {'memory': '8g', 'walltime': '08:00:00'}
    # program creates a temporary folder for temporary files and if it already exists relate won't run
    # remove *.gz becuase relate will not overwrite existing gz files
    spec = f'''
    mkdir -p {dirname(output_base_path)}
    cd {dirname(output_base_path)}
    rm -rf {basename(output_base_path)}
    rm -f {basename(output_base_path)}.anc.gz
    rm -f {basename(output_base_path)}.mut.gz
    {config['relate_dist_dir']}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i {basename(input_base_path)} --poplabels {relpath(poplabels, dirname(output_base_path))} -o {basename(output_base_path)} --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
    sleep 20
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def detect_selection(anc, pop=None, mut=None, poplabels=None):
    """
    detect selection using RELATEs builtin statistic
    .freq --> Records the frequency of the derived allele at generations genN .. gen1
    .lin --> Records the number of lineages in the tree at generations genN .. gen1 as well as the number of lineages when the mutation had frequency 2
    .sele --> Records the log10 p-value for selection evidence at generations genN .. gen1 as well as the log10 p-value when the
    mutation had frequency 2. Log10 p-value is set to 1 if mutation had frequency <= 1 at a generation. 
    """
    # output_dir = f'results/relate/{pop}'
    # file_name_input = joinpath(dirname(anc), 'haplotypes_demog')
    # file_name_output = f'steps/relate/{pop}/haplotypes'

    input_base_path = stripsuf(anc)
    output_base_path = joinpath(dirname(anc), basename(input_base_path) + '_sele')
    inputs = {'anc': anc, 
              'mut': mut, 
              'poplabels': poplabels}
    outputs = {'freq': output_base_path + '.freq', 
               'lin': output_base_path + '.lin', 
               'sele': output_base_path + '.sele'}
    options = {'memory': '20g', 'walltime': '10:00:00'}
    # program creates a temporary folder for temporary files and if it already exists relate won't run
    # remove *.gz becuase relate will not overwrite existing gz files    
    spec = f'''
    mkdir -p {dirname(output_base_path)}
    cd {dirname(output_base_path)}
    rm -rf {basename(output_base_path)}
    {config['relate_dist_dir']}/scripts/DetectSelection/DetectSelection.sh -i {basename(input_base_path)} -m 1.25e-8 --poplabels {relpath(poplabels, dirname(output_base_path))} -o {basename(output_base_path)}
    sleep 80
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def tree_seq(anc, pop=None, mut=None):
    """
    convert to tree sequence file format (tskit)
    this function converts anc/mut files inferred by Relate into the tree sequence file format used by tskit. In the current
    implementation, each tree is stored with new nodes in the tree sequence file format, leading to no compression. In addition,
    information about how long branches persist, and how many mutations map to a branch are lost by this conversion.
    """

    input_base_path = stripsuf(anc)
    output_base_path = joinpath(dirname(anc), basename(input_base_path) + '_trees')
    inputs = {'anc': anc, 'mut': mut}
    outputs = {'trees': output_base_path + '.trees'}
    options = {'memory': '8g', 'walltime': '04:00:00'}
    spec = f'''
    mkdir -p {dirname(output_base_path)}
    cd {dirname(output_base_path)}
    {config['relate_dist_dir']}/bin/RelateFileFormats --mode ConvertToTreeSequence -i {basename(input_base_path)} -o {basename(output_base_path)}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def workflow(working_dir=os.getcwd(), defaults={}, config=None):

    gwf = Workflow(working_dir=working_dir, defaults=defaults)

    # collect targets for use as submodule
    tgt = defaultdict(list)

    # chromosomes = ['chr21', 'chr22', 'chrX'] 
    chromosomes = [f'chr{x}' for x in range(1, 23)] + ['chrX'] 

    tgt[f'download_meta'] = gwf.target_from_template(
        f'download_meta',
        download_meta(config, chromosomes)
    )
    
    for chrom in chromosomes:

        tgt[f'download_data_{chrom}'] = gwf.target_from_template(
            f'download_data_{chrom}',
            download_data(config, chrom)
        )
        
        # tgt['download'] = gwf.target_from_template(f'download',
        #     download_data(config))

        tgt[f'maps_{chrom}'] = gwf.target_from_template(
            f'maps_{chrom}',
            decode_genetic_maps(
                config['decode_hg38_sexavg_per_gen'], 
                f'steps/relate/{chrom}/genetic_map.tsv', chrom
                )
            )

        tgt[f'haploids_{chrom}'] = gwf.target_from_template(
            f'haploids_{chrom}',
            female_haploid(
                tgt[f'download_data_{chrom}'].outputs[f'sample_vcf'], chrom
                )
            )

        tgt[f'ids_{chrom}'] = gwf.target_from_template(
            f'ids_{chrom}', 
            haplotype_id(
                tgt[f'haploids_{chrom}'].outputs[f'haplotypes'], chrom
                )
            )
        tgt[f'all_poplab_{chrom}'] = gwf.target_from_template(
            f'all_poplab_{chrom}',
            all_pop_labels(
                tgt[f'ids_{chrom}'].outputs['ids'], 
                tgt[f'download_meta'].outputs['1000G_2504_seq_index'], 
                tgt[f'download_meta'].outputs['1000G_698_seq_index'],
                chrom
                )
            )
        tgt[f'conv_vcf_{chrom}'] = gwf.target_from_template(
            f'conv_vcf_{chrom}', 
            convert_vcf(
                tgt[f'haploids_{chrom}'].outputs['haplotypes'], 
                tgt[f'all_poplab_{chrom}'].outputs['poplabels'],
                chrom
                )
            )
        """
        African Ancestry in SW U                [ASW]	 62
        African Caribbean in Barbado            [ACB]	120
        Bengali in Banglades                    [BEB]	144
        British From England and Scotlan        [GBR]	100
        Chinese Dai in Xishuangbanna, China     [CDX]	102
        Colombian in Medellín, Colombia         [CLM]	136
        Esan in Nigeria                         [ESN]	173
        Finnish in Finla                        [FIN]	103
        Gambian in Western Division – Mandin    [GWD]	179
        Gujarati Indians in Houston, Texas, USA [GIH]	109
        Han Chinese in Beijing, Chin            [CHB]	120
        Han Chinese Sout                        [CHS]	163
        Iberian Populations in Spain            [IBS]	157
        Indian Telugu in the U.K                [ITU]	118
        Japanese in Tokyo, Japan                [JPT]	120
        Kinh in Ho Chi Minh City, Vietna        [KHV]	124
        Luhya in Webuye, Ken                    [LWK]	120
        Mende in Sierra Leon                    [MSL]	128
        Mexican Ancestry in Los Angeles CA U    [MXL]	 71
        Peruvian in Lima Per                    [PEL]	122
        Puerto Rican in Puerto Rico             [PUR]	139
        Punjabi in Lahore, Pakistan             [PJL]	158
        Sri Lankan Tamil in the                 [STU]	128
        Toscani in Itali                        [TSI]	114
        Yoruba in Ibadan, Nigeri                [YRI]   120
        """

        populations = ['YRI']
        # populations = ['ASW', 'ACB', 'BEB', 'GBR', 'CDX', 'CLM', 'ESN', 'FIN', 'GWD', 
        #             'GIH', 'CHB', 'CHS', 'IBS', 'ITU', 'JPT', 'KHV', 'LWK', 'MSL', 
        #             'MXL', 'PEL', 'PUR', 'PJL', 'STU', 'TSI', 'YRI']

        for pop in populations:
            
            # exlcude related
            related_target = gwf.map(exclude_related,
                                    [(tgt['download_meta'].outputs['1000G_698_seq_index'], pop, chrom)],
                                    name=f"excl_rel_{chrom}_{pop}")
            tgt[f'exclude_related_{chrom}_{pop}'] = related_target
            related = related_target.outputs[0]  # list

            # get ids for other pops
            input_other_pop = [(tgt[f'all_poplab_{chrom}'].outputs['poplabels'], pop, chrom)]
            tgt[f'ids_other_{chrom}_{pop}'] = gwf.map(ids_other_pop, 
                                                    input_other_pop, 
                                                    name=f"ids_other_{chrom}_{pop}")

            # combine related and other pops
            tgt[f'comb_files_{chrom}_{pop}'] = gwf.map(combine_files, 
                                                    tgt[f'ids_other_{chrom}_{pop}'].outputs,
                                                    extra = {'pop': pop, 
                                                             'chrom': chrom,
                                                            'related':related,
                                                            'chrom': chrom
                                                            }, 
                                                    name=f"combine_files_{chrom}_{pop}")

            # list of excluded
            tgt[f'excl_{chrom}_{pop}'] = gwf.map(excluded_list, 
                                                tgt[f'comb_files_{chrom}_{pop}'].outputs, 
                                                extra = {'pop': pop, 
                                                         'chrom': chrom,
                                                        'haplotype_id':tgt[f'ids_{chrom}'].outputs['ids'],
                                                        'chrom': chrom
                                                        }, 
                                                name=f"excl_{chrom}_{pop}")

            # list of included
            tgt[f"poplabels_{chrom}_{pop}"] = gwf.map(pop_labels, 
                                        tgt[f'excl_{chrom}_{pop}'].outputs, 
                                        extra = {'pop': pop, 
                                                 'chrom': chrom,
                                                'poplabels': tgt[f'all_poplab_{chrom}'].outputs['poplabels'],
                                                'chrom': chrom
                                                }, 
                                        name=f"poplabels_{chrom}_{pop}")

            # prepare input for relate
            tgt[f"prepare_{chrom}_{pop}"] = gwf.map(prepare_files, 
                                                tgt[f'excl_{chrom}_{pop}'].outputs, 
                                                extra = {'pop': pop,   
                                                         'chrom': chrom,                                        
                                                            'haps': tgt[f'conv_vcf_{chrom}'].outputs['haps'],
                                                            'sample': tgt[f'conv_vcf_{chrom}'].outputs['sample'],
                                                            'ancestor': tgt[f'download_meta'].outputs[f'ancestral_fa_{chrom}'], 
                                                            'mask':tgt[f'download_data_{chrom}'].outputs['mask'],
                                                            'poplabels': f'steps/relate/{chrom}/poplabels.txt' #######
                                                            },
                                                name=f"prepare_{chrom}_{pop}")

            # run relate
            tgt[f"relate_{chrom}_{pop}"] = gwf.map(relate, 
                                                tgt[f'maps_{chrom}'].outputs, 
                                                extra = {'pop': pop, 
                                                         'chrom': chrom,
                                                        'haps': tgt[f"prepare_{chrom}_{pop}"].outputs[0]['haps'],
                                                        'sample': tgt[f"prepare_{chrom}_{pop}"].outputs[0]['sample'], 
                                                        'annot': tgt[f"prepare_{chrom}_{pop}"].outputs[0]['annot'], 
                                                        'dist': tgt[f"prepare_{chrom}_{pop}"].outputs[0]['dist'],
                                                        }, 
                                                name=f"relate_{chrom}_{pop}")

            
            # estimate pop sizes
            tgt[f"demog_{chrom}_{pop}"] = gwf.map(estimate_pop_sizes, 
                                        # collect(tgt[f"relate_{pop}"], 'anc'),
                                                [tgt[f"relate_{chrom}_{pop}"].outputs[0]['anc']],                                    
                                                extra = {#'pop': pop,  
                                                        'mut': tgt[f"relate_{chrom}_{pop}"].outputs[0]['mut'],
                                                        'poplabels': tgt[f"prepare_{chrom}_{pop}"].outputs[0]['poplabels'],
                                                        }, 
                                                name=f"demog_{chrom}_{pop}")

            # detect selection
            tgt[f"sel_{chrom}_{pop}"] = gwf.map(detect_selection, 
                                            [tgt[f"demog_{chrom}_{pop}"].outputs[0]['anc']],                                           
                                            #    [f'steps/relate/{pop}/haplotypes_demog.anc'], 
                                            extra = {'pop': pop, 
                                                        'mut':  tgt[f"demog_{chrom}_{pop}"].outputs[0]['mut'], # f'steps/relate/{pop}/haplotypes_demog.mut', 
                                                        # 'poplabels': tgt[f"demog_{pop}"].outputs[0]['poplabels'], # f'steps/relate/{pop}/haplotypes_demog.mut'
                                                        'poplabels': tgt[f"prepare_{chrom}_{pop}"].outputs[0]['poplabels'], # f'steps/relate/{pop}/haplotypes_demog.mut'
                                                        }, 
                                            name=f"sel_{chrom}_{pop}")


            # convert to tree sequence file format (tskit),
            tgt[f"trees_{chrom}_{pop}"] = gwf.map(tree_seq, 
                                                #  [f'steps/relate/{pop}/haplotypes.anc'], 
                                                #  tgt[f"relate_{pop}"].outputs[0]['anc'],
                                                [tgt[f"demog_{chrom}_{pop}"].outputs[0]['anc']],
                                                extra = {'pop': pop, 
                                                        'mut': tgt[f"demog_{chrom}_{pop}"].outputs[0]['mut'], #f'steps/relate/{pop}/haplotypes.mut'
                                                        }, 
                                                name=f"trees_{chrom}_{pop}")

    return gwf, tgt



####################################################################
# Use code like this to run this as standalone workflow: 
####################################################################

import yaml
with open('config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

# # reorganize by chrom
# chromosomes = [f'chr{x}' for x in range(1, 23)] + ['chrX'] 
# _config = defaultdict(dict)
# for chrom in chromosomes:
#     for key1, value1 in config_yaml.items():
#         if type(value1) is dict:
#             if chrom in value1:
#                 _config[chrom][key1]  = value1[chrom]
#         else:
#             _config[chrom][key1] = value1
# config = _config # gloabl variable for use in targets
# # is redifined in the workflow function

gwf, targets  = workflow(working_dir=os.getcwd(), 
                        defaults={'account': 'xy-drive'},
                        config=config)

####################################################################
# Use code like this to run this as a submodule workflow: 
####################################################################

# data_dir = '/home/kmt/Primategenomes/data/final_tables'
# state_posterior_files = sorted(Path(data_dir).glob('**/*.HDF'))
# ils = importlib.import_module('primate-ils.workflow')
# gwf, codeml_targets  = ils.workflow(working_dir=working_dir, 
#                                     defaults={'account': 'xy-drive'},
#                                     input_files=state_posterior_files)
# globals()['ils'] = gwf
