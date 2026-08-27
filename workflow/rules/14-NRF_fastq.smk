'''
This Snakefile contains rules to QC and process the raw reads Fastq files:

readsQC_nrf: Perform QC and process raw reads
readsQC_summary_nrf: Aggregates reads QC results
'''

localrules: readsQC_summary_nrf

rule readsQC_nrf:
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
        reads = ancient('rawdata/nrf1/fastq/{sample}.fastq.gz'),
    output:
        trimmed  = 'results/14-NRF_fastq/QC/fastq/{sample}.fq.gz',
        html = 'results/14-NRF_fastq/QC/html/{sample}.html',
        json = 'results/14-NRF_fastq/QC/json/{sample}.json'
    log:
        'logs/14-NRF_fastq/QC_{sample}.log'
    benchmark:
        'benchmarks/14-NRF_fastq/QC_{sample}.txt'
    resources:
        mem_mb = 5000,
        time = '00:5:00'
    params:
        dir='results/14-NRF_fastq/QC'
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

rule readsQC_summary_nrf:
    '''
    R scripts used to agreggate and summarise the QC results from readsQC
    '''

    input:
        ancient(get_QCjson_nrf)
    output:
        plot = 'results/14-NRF_fastq/summaries/QC_plot_density.png',
    log:
        'logs/14-NRF_fastq/QCsummary.log'
    benchmark:
        'benchmarks/14-NRF_fastq/QCsummary.txt'
    resources:
        mem_mb = 1000,
        time = '00:05:00'
    params:
        dir = 'results/14-NRF_fastq/summaries',
        source=rules.readsQC_nrf.params.dir,
        metadata = config['nrf_metadata']
    threads: 1
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}

        echo "Summarising QC results" >> {log}
        Rscript workflow/scripts/14.1-QC_reads_table.R -m {params.metadata} -d {params.source}/json/ -o {params.dir}

        echo "Plotting QC summaries" >> {log}
        Rscript workflow/scripts/14.2-QC_reads_plots.R -m {params.metadata} -f {params.dir}/QC_fastq.csv -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

