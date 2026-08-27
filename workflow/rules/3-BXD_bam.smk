 '''
This Snakefile contains rules to process and QC the alignment files:

sort_aligment: Sort and index BAM files
bamQC: Extract QC information from the BAM files
alignQC_summary: Aggregates aligment QC results
make_wig (optional): produces coverage files for browser visualisation
'''

localrules: alignQC_summary

rule sort_aligment:
    '''
    samtools: SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format
    Reference: https://www.htslib.org/

    Performs: 
    - Convert .sam to binary .bam files
    - Sort files by name and by coordinate
    - Index .bam files
    '''

    input:
        reads_full = ancient(rules.align_atac.output.full),
        reads_sub1 = ancient(rules.align_atac.output.sub1),
        reads_sub2 = ancient(rules.align_atac.output.sub2),
    output:
        full_co = protected('results/3-BXD_bam/sort/coord/{sample}_full.co.bam'),
        sub1_co = protected('results/3-BXD_bam/sort/coord/{sample}_sub1.co.bam'),
        sub2_co = protected('results/3-BXD_bam/sort/coord/{sample}_sub2.co.bam'),
        full_nm = protected('results/3-BXD_bam/sort/name/{sample}_full.nm.bam'),
        sub1_nm = protected('results/3-BXD_bam/sort/name/{sample}_sub1.nm.bam'),
        sub2_nm = protected('results/3-BXD_bam/sort/name/{sample}_sub2.nm.bam')
    log:
        'logs/3-BXD_bam/sort_{sample}.log'
    benchmark:
        'benchmarks/3-BXD_bam/sort_{sample}.txt'
    resources:
        mem_mb = 20000,
        time = '02:00:00'
    threads: 20
    params:
        dir='results/3-BXD_bam/sort'
    shell:
        '''
        module load samtools/1.21

        mkdir -p {params.dir}/name {params.dir}/coord

        echo "Sorting sample <{input.reads_full}> by name" > {log}
        samtools view -S -b {input.reads_full} | samtools sort -n -o {params.dir}/name/{wildcards.sample}_full.nm.bam -@ 20 &>> {log}
        echo "Sorting sample <{input.reads_full}> by coordinate" >> {log}
        samtools sort -o {params.dir}/coord/{wildcards.sample}_full.co.bam -@ 20 {params.dir}/name/{wildcards.sample}_full.nm.bam &>> {log}
        echo "Indexing sample <{input.reads_full}> " >> {log}
        sleep 60
        samtools index -@ 20 {params.dir}/coord/{wildcards.sample}_full.co.bam 2>> {log}
        

        echo "Sorting sample <{input.reads_sub1}> by name" >> {log}
        samtools view -S -b {input.reads_sub1} | samtools sort -n -o {params.dir}/name/{wildcards.sample}_sub1.nm.bam -@ 20 &>> {log}
        echo "Sorting sample <{input.reads_sub1}> by coordinate" >> {log}
        samtools sort -o {params.dir}/coord/{wildcards.sample}_sub1.co.bam -@ 20 {params.dir}/name/{wildcards.sample}_sub1.nm.bam &>> {log}
        echo "Indexing sample <{input.reads_sub1}> " >> {log}
        sleep 60
        samtools index -@ 20 {params.dir}/coord/{wildcards.sample}_sub1.co.bam 2>> {log}
        
        echo "Sorting sample <{input.reads_sub2}> by name" >> {log}
        samtools view -S -b {input.reads_sub2} | samtools sort -n -o {params.dir}/name/{wildcards.sample}_sub2.nm.bam -@ 20 &>> {log}
        echo "Sorting sample <{input.reads_sub2}> by coordinate" >> {log}
        samtools sort -o {params.dir}/coord/{wildcards.sample}_sub2.co.bam -@ 20 {params.dir}/name/{wildcards.sample}_sub2.nm.bam &>> {log}
        echo "Indexing sample <{input.reads_sub2}> " >> {log}
        sleep 60
        samtools index -@ 20 {params.dir}/coord/{wildcards.sample}_sub2.co.bam 2>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

rule make_wig:
    '''
    deeptools: A suite of python tools particularly developed for the efficient analysis of high-throughput sequencing data
    Reference: https://deeptools.readthedocs.io/en/develop/index.html

    Performs: 
    - Convertes BAM to coverage bigwig files
    '''
    input:
        bam = ancient(rules.sort_aligment.output.full_co)
    output:
        bw = 'results/3-BXD_bam/wigs/{sample}.bw',
    log:
        'logs/3-BXD_bam/wigs_{sample}.log'
    benchmark:
        'benchmarks/3-BXD_bam/wigs_{sample}.txt'
    resources:
        mem_mb = 1000,
        time = '00:10:00'
    params:
        dir='results/3-BXD_bam/wigs',
        size=2652783500 #https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    threads: 10
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}
        
        echo "Estimating coverage for <{input.bam}>" >> {log}
        bamCoverage --bam {input.bam} -o {output.bw} --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize {params.size} --numberOfProcessors 10 --ignoreForNormalization chrX &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule bamQC:
    '''
    ATACseqQC: R scripts used to agreegate and summarise the QC results from align_atac
    Reference:  https://doi.org/10.1186/s12864-018-4559-3; https://jianhong.github.io/ATACseqQC/articles/ATACseqQC.html

    Performs: 
    - Calls sRsamtools to read BAM file
    - Performs recomended QC from package
    '''

    input:
        bam = ancient(rules.sort_aligment.output.full_co)
    output:
        plot = protected('results/3-BXD_bam/bamQC/{sample}_full.RData'),
    log:
        'logs/3-BXD_bam/bamQC_{sample}.log'
    benchmark:
        'benchmarks/3-BXD_bam/bamQC_{sample}.txt'
    resources:
        mem_mb = 50000,
        time = '02:00:00'
    threads: 1
    params:
        dir='results/3-BXD_bam/bamQC',
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}

        echo "Summarising QC results" >> {log}
        Rscript workflow/scripts/3.1-QC_align_bamQC.R -b {input.bam} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule alignQC_summary:
    '''
    R scripts used to agreegate and summarise the QC results from align_atac and bamQC

    Performs: 
    - Agreegates QC results
    - Plots summary QC
    '''

    input:
        ancient(get_QCbams)
    output:
        plot = protected('results/3-BXD_bam/alignQC/QCbam_align.csv'),
    log:
        'logs/3-BXD_bam/alignQC.log'
    benchmark:
        'benchmarks/3-BXD_bam/alignQC.txt'
    resources:
        mem_mb = 1000,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/3-BXD_bam/alignQC',
        alignQC='logs/2-BXD_align/',
        bamQC=rules.bamQC.params.dir,
        metadata = config['atac_metadata']
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}

        echo "Summarising alignment results" >> {log}
        Rscript workflow/scripts/3.2-QC_align_QC.R -d {params.alignQC} -m {params.metadata} -o {params.dir}

        echo "Ploting bam summaries" >> {log}
        Rscript workflow/scripts/3.3-QC_align_plots.R -d {params.bamQC} -m {params.metadata} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''
