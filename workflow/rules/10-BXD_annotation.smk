'''
This Snakefile contains rules to annotate regions:

merge_regions: Aggregates contiguous regions into continuous blocks
plot_regions: Get QC statistics for regions
annot2gene: Annotate regions to nearest transcript
get_encode: Download ENCODE cCRE data
annot2encode: Annotate regions to cCRE
plot_volcanos: Plot statistical results using annotation
'''

localrules: merge_regions, plot_regions, annot2gene, get_encode, annot2encode, plot_volcanos

rule merge_regions:
    '''
    csaw: R package to process aggregated genomic windows 
    Reference: https://doi.org/10.1038/nmeth.3252

    Performs:
    - Aggregates small significant regions into continuous blocks
    - Processes summary statistics
    '''
    input:
        dar = ancient('results/8-BXD_differential/DA/atac_peaks_treament.csv'),
        qtl = ancient(rules.qtl_aggregate.output.atac),
        int = ancient('results/9-BXD_qtl/ttest/qtl_ttests_atac.csv'), 
        features = ancient(rules.make_features.output.gtf),
        counts = ancient(lambda wildcards: f"{rules.differential_analysis.params.dir}/atac_counts_disp.RData"),
        metadata =ancient(rules.check_variants.output.meta),
    output:
        dar = protected('results/10-BXD_annotation/reduced/ranges_dar.RData'),
        qtl = protected('results/10-BXD_annotation/reduced/ranges_qtl.RData'),
        int = protected('results/10-BXD_annotation/reduced/ranges_int.RData')
    log:
        'logs/10-BXD_annotation/reduce.log'
    benchmark:
        'benchmarks/10-BXD_annotation/reduce.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 20
    params:
        dir='results/10-BXD_annotation/reduced'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Merge significant regions for annotation" > {log}
        Rscript workflow/scripts/10.1-Annotation_reduce.R -c '{input.counts}' -m {input.metadata} -p {input.features} -d {input.dar} -q {input.qtl} -i {input.int} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule plot_regions:
    '''
    Get and plot some descriptive statistics of the significant regions
    '''
    input:
        dar = ancient(rules.merge_regions.output.dar),
        qtl = ancient(rules.merge_regions.output.qtl),
        int = ancient(rules.merge_regions.output.int), 
        peaks =  ancient(rules.make_features.output.gtf),
        tab = ancient(rules.qtl_aggregate.output.atac),
        tab_fc = ancient(rules.qtl_aggregate.output.atacFC),
    output:
        plot = protected('results/10-BXD_annotation/plots/covplot_DAR.png')
    log:
        'logs/10-BXD_annotation/plots.log'
    benchmark:
        'benchmarks/10-BXD_annotation/plots.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/10-BXD_annotation/plots'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Getting statistics for regions" > {log}
        Rscript workflow/scripts/10.2-Annotation_plot.R -p {input.peaks} -d {input.dar} -q {input.qtl} -i {input.int} -b {input.tab} -f {input.tab_fc} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule annot2gene:
    '''
    ChIPseeker: R package for peak annotation, comparison and visualization
    Reference: https://doi.org/10.1093/bioinformatics/btv145

    Performs:
    - Peak annotation using reference gene models
    '''
    input:
        dar = ancient(rules.merge_regions.output.dar),
        qtl = ancient(rules.merge_regions.output.qtl),
        int = ancient(rules.merge_regions.output.int), 
        peaks =  ancient(rules.make_features.output.gtf),
        check =ancient(rules.plot_regions.output.plot)
    output:
        annotation = protected('results/10-BXD_annotation/gene/ranges_annotated.RData')
    log:
        'logs/10-BXD_annotation/gene.log'
    benchmark:
        'benchmarks/10-BXD_annotation/gene.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/10-BXD_annotation/gene'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Annotating regions to mm10 gene models" > {log}
        Rscript workflow/scripts/10.3-Annotation_2gene.R -p {input.peaks} -d {input.dar} -q {input.qtl} -i {input.int} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule get_encode:
    '''
    Downloads ENCODE candidate cis-regulatory regions
    '''
    output:
        cre  = 'results/10-BXD_annotation/cre/mm10-cCREs.bed',
    log:
        'logs/10-BXD_annotation/cre.log'
    benchmark:
        'benchmarks/10-BXD_annotation/cre.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/10-BXD_annotation/cre'
    shell:
        '''
        mkdir -p {params.dir}

        echo "Downloading ENCODE cCRE" >> {log}
        wget https://downloads.wenglab.org/V3/mm10-cCREs.bed -P {params.dir}/  &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''
rule annot2encode:
    '''
    SCREEN: Search Candidate cis-Regulatory Elements by ENCODE
    Reference: https://doi.org/10.1038/s41586-020-2493-4

    Performs:
    - Peak annotation using ENCODE cCRE data
    '''
    input:
        encode = ancient(rules.get_encode.output.cre),
        annotation = ancient(rules.annot2gene.output.annotation)
    output:
        annotation = protected('results/10-BXD_annotation/encode/cCRE_annotated.RData')
    log:
        'logs/10-BXD_annotation/encode.log'
    benchmark:
        'benchmarks/10-BXD_annotation/encode.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/10-BXD_annotation/encode'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Annotating regions to mm10 ENCODE cCRE" > {log}
        Rscript workflow/scripts/10.4-Annotation_2cre.R -a {input.annotation} -e {input.encode} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule plot_volcanos:
    '''
    Get annotations for volcano plot statistical summaty
    '''
    input:
        deg = ancient('results/8-BXD_differential/DA/rna_stats_treament.RData'),
        qtl = ancient('results/9-BXD_qtl/tables/QTL_table_filter_rna.csv'),
        annotation = ancient(rules.annot2encode.output.annotation),
    output:
        plot = protected('results/10-BXD_annotation/plots/RNA_volcano_genes.png')
    log:
        'logs/10-BXD_annotation/volcano.log'
    benchmark:
        'benchmarks/10-BXD_annotation/volcano.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/10-BXD_annotation/plots'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Plotting statistical results" > {log}
        Rscript workflow/scripts/10.5-Annotation_volcanos.R -a {input.annotation} -r {input.deg} -q {input.qtl} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''
