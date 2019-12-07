configfile: 'config.yaml'
import itertools

def get_abs_path(path):
    import os
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")

outpath = get_abs_path(config['outpath'])

download_d = '{}/downloads'.format(outpath)
working_d  = '{}/{{species}}-{{gene}}'.format(outpath)

species = list(set(itertools.chain.from_iterable([config['reference'][x].keys() for x in config['reference']])))

rule all:
    input:
        expand('{}/{{species}}.{{ref_type}}'.format(download_d), ref_type=config['reference'], species=species),
        expand('{}/{{sample}}.fastq'.format(download_d), sample=config['samples']),
        expand('{}/{{sample}}.{{species}}.{{acids}}.paf'.format(download_d), sample=config['samples'], species=species, acids=['cdna']),
        expand('{}/gene.fastq'.format(working_d), gene=config['genes'], species=species),
        expand('{}/gene.on-motifs.tsv'.format(working_d), gene=config['genes'], species=species),
        expand('{}/{{sample}}.reads.fastq'.format(working_d), gene=config['genes'], sample=config['samples'], species=species),
        expand('{}/{{sample}}.reads.aln.tsv'.format(working_d), gene=config['genes'], sample=config['samples'], species=species),
        expand('{}/{{sample}}.reads-to-gene.paf'.format(working_d), gene=config['genes'], sample=config['samples'], species=species),
        expand('{}/{{sample}}.reads-to-gene.json'.format(working_d), gene=config['genes'], sample=config['samples'], species=species),
        expand('{}/{{sample}}.reads-to-gene.html'.format(working_d), gene=config['genes'], sample=config['samples'], species=species),
        expand('{}/{{sample}}.reads.targets.done'.format(working_d), gene=config['genes'], sample=config['samples'], species=species),

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
        '  chmod -w {output};'
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
        'wget {params.url} -O {output.fastq}; chmod -w {output.fastq}'

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
        'minimap2 {params.mapping_settings} -t {threads} {input.target} {input.reads} > {output.paf} '

rule gene:
    input:
        script = config['exec']['gene'],
        gtf    = '{}/{{species}}.gtf'.format(download_d),
        genome = '{}/{{species}}.dna.fasta'.format(download_d),
    output:
        gene  = '{}/gene.fastq'.format(working_d),
    wildcard_constraints:
        gene='|'.join(config['genes']),
        species='|'.join(species),
    conda:
        'conda.env'
    shell:
        '{input.script} -g {wildcards.gene} -ar {input.gtf} -gr {input.genome} -go {output.gene}'

rule gene_reads:
    input:
        script = config['exec']['gene_reads'],
        reads  = '{}/{{sample}}.fastq'.format(download_d),
        gtf    = '{}/{{species}}.gtf'.format(download_d),
        paf    = '{}/{{sample}}.{{species}}.cdna.paf'.format(download_d),
    output:
        reads = '{}/{{sample}}.reads.fastq'.format(working_d),
        alns  = '{}/{{sample}}.reads.aln.tsv'.format(working_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        gene='|'.join(config['genes']),
        species='|'.join(species),
    conda:
        'conda.env'
    shell:
        '{input.script} -g {wildcards.gene}  -ar {input.gtf} '
        '  -p {input.paf} -r {input.reads}'
        '  -ro {output.reads} -ao {output.alns} '

rule gene_targets:
    input:
        reads = '{}/{{sample}}.reads.fastq'.format(working_d),
    output:
        targets_done ='{}/{{sample}}.reads.targets.done'.format(working_d),
    params:
        prefix ='{}/{{sample}}.read'.format(working_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        gene='|'.join(config['genes']),
        species='|'.join(species),
    run:
        outpath = '{}.{{}}.fasta'.format(params.prefix)
        rids = list()
        print(outpath)
        for idx,line in enumerate(open(input.reads)):
            line = line.rstrip()
            if idx%4 == 0:
                rid = line[1:]
            elif idx%4 == 1:
                seq = line
            elif idx%4 == 2:
                rids.append(rid)
                outfile = open(outpath.format(rid),'w+')
                print('>{}'.format(rid), file=outfile)
                print('{}'.format(seq), file=outfile)
                outfile.close()
        open(output.targets_done, 'w+')


rule gene_read_mapping:
    input:
        reads = '{}/{{sample}}.reads.fastq'.format(working_d),
        gene  = '{}/gene.fastq'.format(working_d),
    output:
        paf = '{}/{{sample}}.reads-to-gene.paf'.format(working_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        gene='|'.join(config['genes']),
        species='|'.join(species),
    params:
        mapping_settings = lambda wildcards: config['mapping_settings']['read-to-gene']
    conda:
        'conda.env'
    shell:
        'minimap2 {params.mapping_settings} -t {threads}  {input.gene} {input.reads} > {output.paf}'

rule target_motifs:
    input:
        script = config['exec']['motifs'],
        gene  = '{}/gene.fastq'.format(working_d),
        motifs  = config['motifs'],
    output:
        tsv = '{}/gene.on-motifs.tsv'.format(working_d),
    wildcard_constraints:
        gene='|'.join(config['genes']),
        species='|'.join(species),
    conda:
        'conda.env'
    shell:
        '{input.script} -m {input.motifs} -g {input.gene} -o {output.tsv}'

rule paf_to_json:
    input:
        script  = config['exec']['paf_to_json'],
        paf     = '{}/{{sample}}.reads-to-gene.paf'.format(working_d),
        target  = '{}/gene.fastq'.format(working_d),
        target_motifs = '{}/gene.on-motifs.tsv'.format(working_d),
        queries = '{}/{{sample}}.reads.fastq'.format(working_d),
        motifs = config['motifs'],
    output:
        json = '{}/{{sample}}.reads-to-gene.json'.format(working_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        gene='|'.join(config['genes']),
        species='|'.join(species),
    conda:
        'conda.env'
    shell:
        '{input.script} -t {input.target} -tm {input.target_motifs} -q {input.queries} -p {input.paf} -m {input.motifs} -jo {output.json}'

rule vis:
    input:
        html = 'vis/genevis.html',
        json = '{}/{{sample}}.reads-to-gene.json'.format(working_d),
    output:
        html = '{}/{{sample}}.reads-to-gene.html'.format(working_d),
    shell:
        'cat '
        '   <(head -n -4 {input.html}) '
        '   {input.json}'
        '   <(tail -n 3 {input.html}) '
        '       > {output.html}'
