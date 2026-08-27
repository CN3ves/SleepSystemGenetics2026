'''
This Snakefile contains rules to perform statistical differential analyses:

differential_analysis: Estimate differential accessibility (ATAC) or expression (RNA)
differential_plots: Plot statistical results from differential analyses
'''

localrules: differential_analysis, differential_plots

rule differential_analysis:
    '''
    edgeR: Bioconductor package for differential expression analysis of digital gene expression data
    Reference: 10.1093/bioinformatics/btp616

    Performs:
    - Differential analyses for effect of SD using BXD strain as covariate (~0 + lines + treatment) using C57xCTRL as baseline.
    - Fits quasi-likelihood (QL) generalization of the negative binomial (NB) model
    - Approximate conditional likelihood is used to estimate feature dispertion using Cox-Reid profile-adjusted likelihood (CR)
    '''
    input:
        counts=ancient(rules.eda.output.counts),
        frip=ancient("results/6-BXD_features/frip/frip_mito.csv")
    output:
        stats=protected('results/8-BXD_differential/DA/{exp}_stats_treament.RData'),
    log:
        'logs/8-BXD_differential/DA_{exp}.log'
    benchmark:
        'benchmarks/8-BXD_differential/DA_{exp}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/8-BXD_differential/DA'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Differential analysis on <{input.counts}>" > {log}
        Rscript workflow/scripts/8.1_Differential_analysis.R -c '{input.counts}' -o {params.dir} 2>> {log}

        echo "Differential analysis on <{input.counts}> excluding low FRiP samples" > {log}
        Rscript workflow/scripts/8.1_Differential_analysis.R -c '{input.counts}' -f {input.frip} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule differential_plots:
    '''
    Plots statistical results from differential analyses 
    '''
    input:
        stats=ancient(rules.differential_analysis.output.stats),
        counts=ancient(rules.eda.output.counts)
    output:
        plot=protected('results/8-BXD_differential/plots/{exp}_coeff_heatFC.png'),
    log:
        'logs/8-BXD_differential/plots_{exp}.log'
    benchmark:
        'benchmarks/8-BXD_differential/plots_{exp}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/8-BXD_differential/plots'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Ploting results from <{input.stats}>" > {log}
        Rscript workflow/scripts/8.2_Differential_plots.R -s '{input.stats}' -c {input.counts} -o {params.dir} 2>> {log}


        frip=$(echo {input.stats} | sed 's/stats/frip_stats/g')
        echo "Ploting results from <$frip>" > {log}
        Rscript workflow/scripts/8.2_Differential_plots.R -s $frip -c {input.counts} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

