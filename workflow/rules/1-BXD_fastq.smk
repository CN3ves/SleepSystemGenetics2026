'''
This Snakefile contains rules to QC and process the raw reads Fastq files:

readsQC: Perform QC and process raw reads
readsQC_summary: Aggregates reads QC results
subsample (optional): subsamples samples to have the same starting number of reads
'''

localrules: readsQC_summary

rule readsQC:
    '''
    Fastp: A tool designed to provide fast all-in-one preprocessing for FastQ files.
    Reference: https://doi.org/10.1093/bioinformatics/bty560; https://github.com/opengene/fastp
    
    Performs: 
    -Read quality
    -Read Trimming (adaptors, quality)
    -Read deduplication
    -Post-trimming quality checks
    '''

    input:
        reads = ancient('rawdata/atac/fastq/{sample}.fastq.gz'),
    output:
        trimmed  = 'results/1-BXD_fastq/QC/fastq/{sample}.fq.gz',
        html = 'results/1-BXD_fastq/QC/html/{sample}.html',
        json = 'results/1-BXD_fastq/QC/json/{sample}.json'
    log:
        'logs/1-BXD_fastq/QC_{sample}.log'
    benchmark:
        'benchmarks/1-BXD_fastq/QC_{sample}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    params:
        dir='results/1-BXD_fastq/QC'
    threads: 16
    shell:
        '''
        module load fastp/1.0.1
        mkdir -p {params.dir}

        echo "Running Fastp on <{input.reads}>" > {log}
        
        fastp -w 16 -D -i {input.reads} -o {output.trimmed} -h {output.html} -j {output.json} &>> {log}
        
        echo "Processed fastq files saved in <{output.trimmed}>" >> {log}
        echo "Reports saved in <{output.html}> and <{output.json}>" >> {log}
        echo "Logs saved in <{log}>" >> {log}
        '''

rule readsQC_summary:
    '''
    R scripts used to agreggate and summarise the QC results from readsQC
    '''

    input:
        ancient(get_QCjson)
    output:
        plot = 'results/1-BXD_fastq/summaries/QC_plot_density.png',
    log:
        'logs/1-BXD_fastq/QCsummary.log'
    benchmark:
        'benchmarks/1-BXD_fastq/QCsummary.txt'
    resources:
        mem_mb = 1000,
        time = '00:05:00'
    params:
        dir = 'results/1-BXD_fastq/summaries',
        source=rules.readsQC.params.dir,
        metadata = config['atac_metadata']
    threads: 1
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}

        echo "Summarising QC results" >> {log}
        Rscript workflow/scripts/2.1-QC_reads_table.R -m {params.metadata} -d {params.source}/json/ -o {params.dir}

        echo "Plotting QC summaries" >> {log}
        Rscript workflow/scripts/2.2-QC_reads_plots.R -m {params.metadata} -f {params.dir}/QC_fastq.csv -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule subsample:
    '''
    Seqtk: Toolkit for processing sequences in FASTA/Q formats
    Reference: https://github.com/lh3/seqtk

    Performs: 
    - Fastq subsampling
    '''

    input:
        reads = ancient(rules.readsQC.output.trimmed),
        qc_check = ancient(rules.readsQC_summary.output.plot)
    output:
        full  = protected('results/1-BXD_fastq/subsampling/{sample}_full.fq.gz'),
        sub1 = protected('results/1-BXD_fastq/subsampling/{sample}_sub1.fq.gz'),
        sub2 = protected('results/1-BXD_fastq/subsampling/{sample}_sub2.fq.gz')
    log:
        'logs/1-BXD_fastq/sampling_{sample}.log'
    benchmark:
        'benchmarks/1-BXD_fastq/sampling_{sample}.txt'
    resources:
        mem_mb = 15000,
        time = '00:30:00'
    threads: 1
    params:
        dir = 'results/1-BXD_fastq/subsampling/',
        #cat logs/1-BXD_QC/*| grep -i 'reads passed filter'| sed 's/.*: //g' | sort -un |head 
        #shows min reads 23429582 
        sub1_reads= 20000000,# 85.36% of minimum
        sub2_reads= 15000000 # 64.02% of minimum

    shell:
        '''
        module load seqtk/1.4
        mkdir -p {params.dir}

        echo "Copy full fastq file <{input.reads}> to {output.full}" >> {log}
        cp {input.reads} {output.full}

        echo "Subseting <{input.reads}> to <{output.sub1}" >> {log}
        seqtk sample -s10 {input.reads} {params.sub1_reads} | gzip > {output.sub1} 2>> {log}
 
        echo "Subseting <{input.reads}> to <{output.sub2}" >> {log}
        seqtk sample -s100 {input.reads} {params.sub2_reads} | gzip > {output.sub2} 2>> {log}    

        echo "Logs saved in <{log}>" >> {log}
        '''
