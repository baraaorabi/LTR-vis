configfile: 'config.yaml'
import itertools

def get_abs_path(path):
    import os
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")

outpath = get_abs_path(config['outpath'])

download_d   = '{}/downloads'.format(outpath)

species = list(set(itertools.chain.from_iterable([config['reference'][x].keys() for x in config['reference']])))

rule all:
    input:
        expand('{}/{{species}}.{{ref_type}}'.format(download_d), ref_type=config['reference'], species=species),
        expand('{}/{{sample}}.fastq'.format(download_d), sample=config['samples']),
        expand('{}/{{sample}}.{{species}}.{{acids}}.paf'.format(download_d), sample=config['samples'], species=species, acids=['cdna']),

rule download_ref:
    output:
        '{}/{{species}}.{{ref_type}}'.format(download_d)
    wildcard_constraints:
        species='|'.join(species),
        ref_type='|'.join(config['reference']),
    params:
        url=lambda wildcards : config['reference'][wildcards.ref_type][wildcards.species]
    conda:
        'conda.env'
    shell:
        'wget {params.url} -O {output}.gz;'
        '  zcat {output}.gz > {output};'
        '  rm {output}.gz'

rule download_sample:
    output:
        fastq='{}/{{sample}}.fastq'.format(download_d)
    wildcard_constraints:
        sample='|'.join(config['samples']),
    params:
        url=lambda wildcards : config['samples'][wildcards.sample]
    conda:
        'conda.env'
    shell:
        'wget {params.url} -O {output.fastq}'

rule map_reads:
    input:
        reads  = '{}/{{sample}}.fastq'.format(download_d),
        target = '{}/{{species}}.{{acids}}.fasta'.format(download_d),
    output:
        paf = '{}/{{sample}}.{{species}}.{{acids}}.paf'.format(download_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        species='|'.join(species),
        acids='cdna|dna',
    params:
        mapping_settings = lambda wildcards: config['mapping_settings'][wildcards.acids]
    conda:
        'conda.env'
    threads:
        32
    shell:
        'minimap2 -x {params.mapping_settings} -t {threads} {input.target} {input.reads} > {output.paf} '

rule gene_reads:
    input:
        script = config['gene_reads'],
        reads  = '{}/{{sample}}.fastq'.format(download_d),
        gtf    = '{}/{{species}}.gtf'.format(download_d),
        genome = '{}/{{species}}.dna.fasta'.format(download_d),
        paf    = '{}/{{sample}}.{{species}}.cdna.paf'.format(download_d),
    output:
        gene = '{}/{{species}}.{{gene}}.fasta'.format(download_d),
        reads  = '{}/{{sample}}.{{species}}.{{gene}}.fastq'.format(download_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        gene='|'.join(config['genes']),
        species='|'.join(species),
    conda:
        'conda.env'
    shell:
        ''
