 '''
This Snakefile contains rules to process and QC the alignment files:

sort_aligment_nrf: Sort and index BAM files
bamQC_nrf: Extract QC information from the BAM files
alignQC_summary_nrf: Aggregates aligment QC results
eda_nrf: Normalization of read counts
qualimap: QC dignostics per sample
qualimap_aggregate: Aggregated QC diagnostics
multi_QC2: Re-run multiQC for summary including RNA-seq
'''

localrules: alignQC_summary_nrf, qualimap, qualimap_aggregate, multi_QC2

rule sort_aligment_nrf:
    '''
    samtools: SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format
    Reference: https://www.htslib.org/

    Performs: 
    - Convert .sam to binary .bam files
    - Sort files by name and by coordinate
    - Index .bam files
    '''

    input:
        sam = ancient(rules.align_rna.output.full)
    output:
        bam = protected('results/16-NRF_bam/sort/{sample}.co.bam')
    log:
        'logs/16-NRF_bam/sort_{sample}.log'
    benchmark:
        'benchmarks/16-NRF_bam/sort_{sample}.txt'
    resources:
        mem_mb = 20000,
        time = '00:10:00'
    threads: 20
    params:
        dir='results/16-NRF_bam/sort'
    shell:
        '''
        module load samtools/1.21

        mkdir -p {params.dir}

        echo "Sorting sample <{input.sam}> by coordinate" > {log}
        samtools sort -@ 20 -o {params.dir}/{wildcards.sample}.co.bam {input.sam} &>> {log}
        
        echo "Indexing sample <{input.sam}> " >> {log}
        samtools index -@ 20 {params.dir}/{wildcards.sample}.co.bam 2>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

rule bamQC_nrf:
    '''
    ATACseqQC: R scripts used to agreegate and summarise the QC results from align_atac
    Reference:  https://doi.org/10.1186/s12864-018-4559-3; https://jianhong.github.io/ATACseqQC/articles/ATACseqQC.html

    Performs: 
    - Calls sRsamtools to read BAM file
    - Performs recomended QC from package
    '''

    input:
        bam = ancient(rules.sort_aligment_nrf.output.bam)
    output:
        plot = protected('results/16-NRF_bam/bamQC/{sample}.RData'),
    log:
        'logs/16-NRF_bam/bamQC_{sample}.log'
    benchmark:
        'benchmarks/16-NRF_bam/bamQC_{sample}.txt'
    resources:
        mem_mb = 50000,
        time = '02:00:00'
    threads: 1
    params:
        dir='results/16-NRF_bam/bamQC',
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}

        echo "Summarising QC results" >> {log}
        Rscript workflow/scripts/16.1-QC_align_bamQC.R -b {input.bam} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule alignQC_summary_nrf:
    '''
    R scripts used to agreegate and summarise the QC results from align_atac and bamQC

    Performs: 
    - Agreegates QC results
    - Plots summary QC
    '''

    input:
        ancient(get_QCbams_nrf)
    output:
        plot = protected('results/16-NRF_bam/alignQC/QCbam_align.csv'),
    log:
        'logs/16-NRF_bam/alignQC.log'
    benchmark:
        'benchmarks/16-NRF_bam/alignQC.txt'
    resources:
        mem_mb = 1000,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/16-NRF_bam/alignQC',
        alignQC='results/15-NRF_align/aligment/',
        bamQC=rules.bamQC_nrf.params.dir,
        metadata = config['nrf_metadata']
    shell:
        '''
        module load  r-light/4.5.2
        mkdir -p {params.dir}

        echo "Summarising alignment results" > {log}
        Rscript workflow/scripts/16.2-QC_align_QC.R -d {params.alignQC} -m {params.metadata} -o {params.dir}

        echo "Ploting bam summaries" >> {log}
        Rscript workflow/scripts/16.3-QC_align_plots.R -d {params.bamQC} -m {params.metadata} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule eda_nrf:
    '''
    EDAseq: Numerical and graphical summaries of RNA-Seq read data.
    Reference: https://bioconductor.org/packages/release/bioc/vignettes/EDASeq/inst/doc/EDASeq.html

    Performs:
    - withinLaneNormalization and betweenLaneNormalization quantile normalization
    '''
    input:
        fastqs = ancient(lambda wildcards: [f'results/14-NRF_fastq/QC/fastq/{sample}.fq.gz' for sample in config['nrf_samples']]),
        bams = ancient(lambda wildcards: [f'results/16-NRF_bam/sort/{sample}.co.bam' for sample in config['nrf_samples']]),
        counts= ancient(lambda wildcards: [f'results/15-NRF_align/aligment/{sample}_ReadsPerGene.out.tab' for sample in config['nrf_samples']]),
    output:
        counts=protected('results/16-NRF_bam/EDA/EDA_normalised.RData'),
    log:
        'logs/16-NRF_bam/eda.log'
    benchmark:
        'benchmarks/16-NRF_bam/eda.txt'
    resources:
        mem_mb = 30000,
        time = '03:15:00'
    threads: 1
    params:
        dir='results/16-NRF_bam/EDA',
        metadata = config['nrf_metadata']
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Running QC diagnostics and normalizing samples" > {log}
        Rscript workflow/scripts/16.4_QC_eda.R -f '{input.fastqs}' -b '{input.bams}' -c '{input.counts}' -m {params.metadata} -o {params.dir} 2>> {log}
 
        echo "Logs saved in <{log}>" >> {log}
        '''

rule qualimap:
    '''
    Qualimap: Quality control of alignment sequencing data
    Reference: http://qualimap.conesalab.org/doc_html/index.html

    Performs:
    - fast analysis across the reference genome of mapping coverage and nucleotide distribution;
    - easy-to-interpret summary of the main properties of the alignment data;
    - analysis of the reads mapped inside/outside of the regions defined in an annotation reference;
    - computation and analysis of read counts obtained from intersting of read alignments with genomic features;
    - analysis of the adequacy of the sequencing depth in RNA-seq experiments;
    - support for multi-sample comparison for alignment data and counts data;
    - clustering of epigenomic profiles.
    '''
    input:
        bam=ancient(rules.sort_aligment_nrf.output.bam),
        reference=ancient(rules.star_index.output.gtf)
    output:
        report='results/16-NRF_bam/qualimap/{sample}/rnaseq_qc_results.txt',
    log:
        'logs/16-NRF_bam/qualimap_{sample}.log'
    benchmark:
        'benchmarks/16-NRF_bam/qualimap_{sample}.txt'
    resources:
        mem_mb = 10000,
        time = '00:45:00'
    threads: 10
    params:
        dir='results/16-NRF_bam/qualimap',
    shell:
        '''
        module load qualimap/2.3

        mkdir -p {params.dir}

        echo "Qualimap bamqc for <{input.bam}>" > {log}
        qualimap bamqc -bam {input.bam} -c -gff {input.reference} -outdir {params.dir}/{wildcards.sample} -gd MOUSE -nt 10 -os &>> {log}

        echo "Qualimap rnaseq for <{input.bam}>" >> {log}
        qualimap rnaseq -bam {input.bam} -gtf {input.reference} -outdir {params.dir}/{wildcards.sample} --java-mem-size=20G &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule qualimap_aggregate:
    '''
    Continuation of rule qualimap for all samples instead of per sample
    '''
    input:
        qualimap= ancient(lambda wildcards: [f'results/16-NRF_bam/qualimap/{sample}/rnaseq_qc_results.txt' for sample in config['nrf_samples']]),
        bams = ancient(lambda wildcards: [f'results/16-NRF_bam/sort/{sample}.co.bam' for sample in config['nrf_samples']]),
        counts= ancient(lambda wildcards: [f'results/15-NRF_align/aligment/{sample}_ReadsPerGene.out.tab' for sample in config['nrf_samples']]),
    output:
        report='results/16-NRF_bam/qualimap/bamqc_multi/multisampleBamQcReport.html',
    log:
        'logs/16-NRF_bam/qualimap.log'
    benchmark:
        'benchmarks/16-NRF_bam/qualimap.txt'
    resources:
        mem_mb = 20000,
        time = '01:00:00'
    threads: 10
    params:
        dir='results/16-NRF_bam/qualimap',
        metadata = config['nrf_metadata']
    shell:
        '''
        module load qualimap/2.3 r-light/4.5.2

        mkdir -p {params.dir}

        echo "Exonic reads: $(cat {params.dir}/*/rnaseq_qc_results.txt | grep exonic | cut -d '(' -f 2 | cut -d ')' -f1 | sort -u | head -n3 )" >> {log}
        Rscript workflow/scripts/16.5_qualimap.R -b '{input.bams}' -c '{input.counts}' -m {params.metadata} -o {params.dir} 

        export R_LIBS='/scratch/csousane/atac/Rlibs/'
        echo "Qualimap counts per group" >> {log}

        for geno in CT FV FT
        do
            qualimap counts -c -d {params.dir}/"qualimap_"$geno"samples.txt"  -k 5 -outdir {params.dir}/$geno -s MOUSE &>> {log}
        done

        echo "Qualimap counts" >> {log}
        qualimap counts -c -d {params.dir}/"qualimap_samples.txt"  -k 5 -outdir {params.dir}/all -s MOUSE &>> {log}

        echo "Qualimap multi-bamqc" >> {log}
        qualimap multi-bamqc -c  -d {params.dir}/qualimap_multi.txt -outdir {params.dir}/bamqc_multi &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule multi_QC2:
    '''
    MultiQC: Tool to automatically summarize analysis results for multiple tools and samples in a single report
    Reference: https://doi.org/10.1093/bioinformatics/btw354

    Performs:
    - Check log files and identifies tool used
    - Collect QC results from identified steps
    - Produces a single report summarising all QC
    '''
    input:
        check=rules.alignQC_summary_nrf.output.plot,
        qualimaps=rules.qualimap_aggregate.output.report
    output:
        report='results/16-NRF_bam/MultiQC/multiqc_report.html',
    log:
        'logs/16-NRF_bam/MultiQC.log'
    benchmark:
        'benchmarks/16-NRF_bam/MultiQC.txt'
    resources:
        mem_mb = 2000,
        time = '01:00:00'
    threads: 1
    params:
        dir='results/16-NRF_bam/MultiQC',
    shell:
        '''
        mkdir -p {params.dir}
   
        echo "Running MultiQC on all logs" > {log}
        multiqc --interactive ./ --ignore Rlibs --ignore libs --ignore workflow --ignore rawdata &>> {log}

        echo "Moving results to <{params.dir}>" >> {log}
        mv -f multiqc* {params.dir}/
        echo "Logs saved in <{log}>" >> {log}
        '''

