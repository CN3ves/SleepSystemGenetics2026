'''
This Snakefile contains rules to align processed reads to reference genome:

bowtie_index: Downloads the appropriate annotation
align_atac: Aligned reads to genome
'''
localrules: bowtie_index

rule bowtie_index:
    '''
    Downloads bowtie2 genome index for alignment
    '''

    output:
        idx  = 'results/2-BXD_align/index/mm10.1.bt2',
    log:
        'logs/2-BXD_align/index.log'
    benchmark:
        'benchmarks/2-BXD_align/index.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/2-BXD_align/index',
        genome = config['atac_genome']

    shell:
        '''
        mkdir -p {params.dir}

        echo "Downloading index for {params.genome}" > {log}
        wget https://genome-idx.s3.amazonaws.com/bt/{params.genome}.zip -P {params.dir}  &>> {log}

        echo "Unpacking genome index" >> {log}
        unzip {params.dir}/{params.genome}.zip -d {params.dir}
        rm {params.dir}/{params.genome}.zip

        echo "Logs saved in <{log}>" >> {log}
        '''

rule align_atac:
    '''
    bowtie2: Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences 
    Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/

    Performs: 
    - Read aligment to reference genome
    '''

    input:
        reads_full = ancient(rules.subsample.output.full),
        reads_sub1 = ancient(rules.subsample.output.sub1),
        reads_sub2 = ancient(rules.subsample.output.sub2),
        idx =  ancient(rules.bowtie_index.output.idx)
    output:
        full = protected('results/2-BXD_align/aligment/{sample}_full.sam'),
        sub1 = protected('results/2-BXD_align/aligment/{sample}_sub1.sam'),
        sub2 = protected('results/2-BXD_align/aligment/{sample}_sub2.sam')
    log:
        'logs/2-BXD_align/align_{sample}.log'
    benchmark:
        'benchmarks/2-BXD_align/align_{sample}.txt'
    resources:
        mem_mb = 5000,
        time = '01:30:00'
    threads: 20
    params:
        dir='results/2-BXD_align/aligment',
        index=rules.bowtie_index.params.dir,
        genome = rules.bowtie_index.params.genome
    shell:
        '''
        module load bowtie2/2.5.4

        mkdir -p {params.dir}

        echo "Aligning full sample <{input.reads_full}> to <{params.genome}> genome annotation" > {log}
        bowtie2 --time --phred33 --local --very-sensitive-local --threads 20 --seed 87 -x {params.index}/{params.genome} -U {input.reads_full} -S {output.full} 2>> {log}

        echo "Aligning subsample <{input.reads_sub1}> to <{params.genome}> genome annotation" >> {log}
        bowtie2 --time --phred33 --local --very-sensitive-local --threads 20 --seed 87 -x {params.index}/{params.genome} -U {input.reads_sub1} -S {output.sub1} 2>> {log}

        echo "Aligning subsample <{input.reads_sub2}> to <{params.genome}> genome annotation" >> {log}
        bowtie2 --time --phred33 --local --very-sensitive-local --threads 20 --seed 87 -x {params.index}/{params.genome} -U {input.reads_sub2} -S {output.sub2} 2>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

