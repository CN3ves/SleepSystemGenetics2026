'''
This Snakefile contains rules to normalised and filter counts:

atac_counts: Builds count R object for downstream analyses using ATAC-seq count files
rna_counts: Builds count R object using RNA-seq count matrix to match ATAC-seq metadata
filter_counts: Filter counts by expression
eda: Exploratory data analyses
'''

localrules: atac_counts, rna_counts, filter_counts, eda

rule atac_counts:
    '''
    edgeR: Bioconductor package for differential expression analysis of digital gene expression data
    Reference: 10.1093/bioinformatics/btp616

    Performs:
    - Collects all ATAC-seq counts results into an R object
    - Adds relevant metadata to the count R object
    '''
    input:
        counts=ancient(get_counts),
        meta= ancient(rules.check_variants.output.meta),
        check=ancient(rules.multi_QC.output.report)
    output:
        object=protected('results/7-BXD_normalization/count_objects/atac_counts.RData'),
    log:
        'logs/7-BXD_normalization/counts_atac.log'
    benchmark:
        'benchmarks/7-BXD_normalization/counts_atac.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/7-BXD_normalization/count_objects'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Loading ATAC-seq read counts into R object" > {log}
        Rscript workflow/scripts/7.0_Normalization_atac.R -c '{input.counts}' -m {input.meta} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule rna_counts:
    '''
    Same as atac_counts, but using the RNA-seq count matrix from Diessler, 2018
    Reference: 10.1371/journal.pbio.2005750 
    '''
    input:
        counts=ancient('rawdata/rna_bxd/RNAseq.txt'),
        meta= ancient('rawdata/rna_bxd/RNAseq_metadata.txt'),
        check=rules.multi_QC.output.report
    output:
        object=protected('results/7-BXD_normalization/count_objects/rna_counts.RData'),
    log:
        'logs/7-BXD_normalization/counts_rna.log'
    benchmark:
        'benchmarks/7-BXD_normalization/counts_rna.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/7-BXD_normalization/count_objects'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Loading RNA-seq count matrix into R object" > {log}
        Rscript workflow/scripts/7.1_Normalization_rna.R -c {input.counts} -m {input.meta} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule filter_counts:
    '''
    edgeR: Bioconductor package for differential expression analysis of digital gene expression data
    Reference: 10.12688/f1000research.8987.2

    Performs:
    - Filters features by expression levels base on sample group (BXD line and treatment)
    - Default values: At least 10 count-per-million in 70% of the samples (2 samples) and more than 15 count-per-million cross all samples
    '''
    input:
        counts=ancient(lambda wildcards: f"results/7-BXD_normalization/count_objects/{wildcards.exp}_counts.RData")
    output:
        object=protected('results/7-BXD_normalization/filter/{exp}_filtered_counts.RData'),
    log:
        'logs/7-BXD_normalization/filter_{exp}.log'
    benchmark:
        'benchmarks/7-BXD_normalization/filter_{exp}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/7-BXD_normalization/filter'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Loading RNA-seq count matrix into R object" > {log}
        Rscript workflow/scripts/7.2_Normalization_filter.R -c {input.counts} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule eda:
    '''
    Exploratory analyses of the datasets

    Performs:
    - TMM normalization
    - Dimension reduction analyses
    - Sample clustering
    - Correlation analyses for confounding factor 
    '''
    input:
        counts=ancient(rules.filter_counts.output.object),
        meta=ancient(rules.check_variants.output.meta)
    output:
        counts=protected('results/7-BXD_normalization/EDA/{exp}_filtered_normalised_counts.RData'),
        plot=protected('results/7-BXD_normalization/EDA/{exp}_QC_covars.png')
    log:
        'logs/7-BXD_normalization/eda_{exp}.log'
    benchmark:
        'benchmarks/7-BXD_normalization/eda_{exp}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/7-BXD_normalization/EDA',
        annotation="Strain,Treatment"
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Running multidimentional scaling on <{wildcards.exp}>" > {log}
        Rscript workflow/scripts/7.4_Normalization_MDS.R -c {input.counts} -o {params.dir} 2>> {log}

        echo "Running clustering heatmaps on <{wildcards.exp}>" >> {log}
        Rscript workflow/scripts/7.5_Normalization_heatmap.R -c {output.counts} -m {input.meta} -a {params.annotation} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''
