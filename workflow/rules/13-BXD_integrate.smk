'''
This Snakefile contains rules to correlate RNA-seq and ATAC-seq results

correlation_annotation: Calculate Pearson correlation betweeen ATAC and RNA based on annotated genes
correlation_qtl: Calculate Pearson correlation  betweeen ATAC and RNA based on QTL
correlation_plots: Plot correlation results
granie: Run Gene Regulatory Network Inference including Enhancers
'''

localrules: correlation_qtl, correlation_plots

rule correlation_annotation:
    '''
    Calculate Pearson correlation betweeen chromatin accessibility and the transcripts of the annotated genes 
    '''
    input:
        rna=ancient(lambda wildcards: f"{rules.differential_analysis.params.dir}/rna_counts_disp.RData"),
        rna_diff=ancient(lambda wildcards: f"{rules.differential_analysis.params.dir}/rna_peaks_treament.csv"),
        atac=ancient(lambda wildcards: f"{rules.differential_analysis.params.dir}/atac_counts_disp.RData"),
        atac_diff=ancient(lambda wildcards: f"{rules.differential_analysis.params.dir}/atac_peaks_treament.csv"),
        regions=rules.make_features.output.gtf
    output:
        cors=protected("results/13-BXD_integrate/cors/gene_region_cors_annotation.csv"),
    log:
        'logs/13-BXD_integrate/cors.log'
    benchmark:
        'benchmarks/13-BXD_integrate/cors.txt'
    resources:
        mem_mb = 50000,
        time = '0:20:00'
    threads: 40
    params:
        dir='results/13-BXD_integrate/cors'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Calculating RNA/ATAC correlations based on annotation" > {log}
        Rscript workflow/scripts/13.1-Integrate_cors.R -a {input.atac} -c {input.atac_diff} -r {input.rna} -t {input.rna_diff} -p {input.regions} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule correlation_qtl:
    '''
    Calculate Pearson correlation betweeen chromatin accessibility and transcriptssharing interaction QTL peaks 
    '''
    input:
        rna=ancient("results/8-BXD_differential/DA/rna_counts_disp.RData"),
        rna_qtl=ancient('results/9-BXD_qtl/ttest/qtl_ttests_rna.csv'),
        atac=ancient("results/8-BXD_differential/DA/atac_counts_disp.RData"),
        atac_qtl=ancient('results/9-BXD_qtl/ttest/qtl_ttests_atac.csv'),
        footprint_qtl=ancient(rules.footprint_qtl.output.qtls),
        regions=rules.make_features.output.gtf
    output:
        cors=protected("results/13-BXD_integrate/qtl/gene_region_cors_qtl.csv")
    log:
        'logs/13-BXD_integrate/qtl.log'
    benchmark:
        'benchmarks/13-BXD_integrate/qtl.txt'
    resources:
        mem_mb = 10000,
        time = '0:10:00'
    threads: 1
    params:
        dir='results/13-BXD_integrate/qtl'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Calculating RNA/ATAC correlations based on annotation" > {log}
        Rscript workflow/scripts/13.2-Integrate_qtl.R -a {input.atac} -c {input.atac_qtl} -r {input.rna} -t {input.rna_qtl} -p {input.regions} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule correlation_plots:
    '''
    Plot overlaps for ATAC and RNA-seq QTL and differential results.
    '''
    input:
        rna_diff=ancient("results/8-BXD_differential/DA/rna_peaks_treament.csv"),
        rna_qtl=ancient(rules.qtl_aggregate.output.rna),
        rnaFC_qtl=ancient(rules.qtl_aggregate.output.rnaFC),
        atac_annot=ancient('results/10-BXD_annotation/encode/cCRE_annotated.RData'),
        atac_qtl=ancient(rules.qtl_aggregate.output.atac),
        atacFC_qtl=ancient(rules.qtl_aggregate.output.atacFC),
        sleep_qtl=ancient(rules.qtl_aggregate.output.sleep),
        annot=ancient(rules.correlation_annotation.output.cors),
        qtl=ancient(rules.correlation_qtl.output.cors),
    output:
        stats=protected('results/13-BXD_integrate/plots/cors_int.png'),
    log:
        'logs/13-BXD_integrate/plots.log'
    benchmark:
        'benchmarks/13-BXD_integrate/plots.txt'
    resources:
        mem_mb = 10000,
        time = '0:10:00'
    threads: 1
    params:
        dir='results/13-BXD_integrate/plots'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Calculating RNA/ATAC correlations based on annotation" > {log}
        Rscript workflow/scripts/13.3-Integrate_plot.R -r {input.rna_diff} -q {input.rna_qtl} -i {input.rnaFC_qtl}  -y {input.atac_annot} -a {input.atac_qtl} -c {input.atacFC_qtl} -s {input.sleep_qtl}  -z {input.annot} -x {input.qtl} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''
 
rule granie:
    '''
    GRaNIE: Gene Regulatory Network Inference including Enhancers
    Reference: https://grp-zaugg.embl-community.io/GRaNIE/; https://doi.org/10.15252/msb.202311627

    Performs: 
    - Inference of the Gene Regulatory Network based on ATAC-seq and RNA-seq raw data
    '''
    input:
        atac=ancient(rules.atac_counts.output.object),
        rna=ancient(rules.rna_counts.output.object),
        regions=ancient(rules.make_features.output.gtf)
    output:
        stats=protected('results/13-BXD_integrate/granie/graph.RData'),
    log:
        'logs/13-BXD_integrate/granie.log'
    benchmark:
        'benchmarks/13-BXD_integrate/granie.txt'
    resources:
        mem_mb = 100000,
        time = '4:00:00'
    threads: 20
    params:
        dir='results/13-BXD_integrate/granie'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Calculating RNA/ATAC correlations based on annotation" > {log}
        Rscript workflow/scripts/13.4-Integrate_granie.R -r {input.rna} -a {input.atac} -p {input.regions} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''
 

 
