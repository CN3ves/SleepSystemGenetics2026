'''
This Snakefile contains rules to produce the figures for the manuscript
'''
localrules: Figure1, FigureS1, Figure2, FigureS2, Figure3, FigureS3, FigureS4

rule Figure1:
    '''
    Figure 1: Acute sleep deprivation leads to genome wide changes in chromatin accessibility
    '''
    input:
        S1=ancient(rules.tableS1.output.cre),
        S2=ancient(rules.tableS2.output.table),
        S3=ancient(rules.tableS3.output.enrich)
    output:
        fig1a=protected('manuscript/figures/data/1a.csv'),
        fig1b=protected('manuscript/figures/data/1b.csv'),
        fig1c=protected('manuscript/figures/data/1c.csv'),
    log:
        'logs/20-Figures/F1.log'
    benchmark:
        'benchmarks/20-Figures/F1.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        data='results/19-Tables/data/',
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure 1a" > {log}
        Rscript workflow/scripts/20.Figure_1a.R  -a {input.S2} -o {params.dir}

        echo "Make Figure 1b" >> {log}
        Rscript workflow/scripts/20.Figure_1b.R  -a {input.S2} -b {input.S1} -o {params.dir}

        echo "Make Figure 1c" >> {log}
        Rscript workflow/scripts/20.Figure_1c.R  -a {input.S3} -o {params.dir}


        echo "Logs saved in <{log}>" >> {log}
        '''

rule FigureS1:
    '''
    Figure S1
    '''
    input:
        S1=ancient(rules.tableS1.output.cre),
        S2=ancient(rules.tableS2.output.table),
        S3=ancient(rules.tableS3.output.table),
        f=ancient(rules.Figure1.output.fig1a)
    output:
        figS1a=protected('manuscript/figures/data/S1a.csv'),
        figS1b=protected('manuscript/figures/data/S1b.csv'),
        figS1c=protected('manuscript/figures/data/S1c.csv'),
        figS1d=protected('manuscript/figures/data/S1d.csv'),
        figS1e=protected('manuscript/figures/data/S1e.csv')
    log:
        'logs/20-Figures/FS1.log'
    benchmark:
        'benchmarks/20-Figures/FS1.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        vars='results/4-BXD_genotype',
        meta=config['atac_metadata']
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure S1a" > {log}
        Rscript workflow/scripts/20.Figure_1Sa.R  -a {params.vars}/variants/ -b {params.vars}/genotypes/ -c {params.meta} -o {params.dir}

        echo "Make Figure S1b" >> {log}
        Rscript workflow/scripts/20.Figure_1Sb.R  -a {input.S2} -o {params.dir}

        echo "Make Figure S1c" >> {log}
        Rscript workflow/scripts/20.Figure_1Sc.R  -a {input.S2} -o {params.dir}
    
        echo "Make Figure S1d" >> {log}
        Rscript workflow/scripts/20.Figure_1Sd.R  -a {input.S2} -b {input.S3} -o {params.dir}

        echo "Make Figure S1e" >> {log}
        Rscript workflow/scripts/20.Figure_1Se.R -a {input.S1} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule Figure2:
    '''
    Figure 2: Transcriptional and signaling changes upstream of the chromatin-accessibility changes.
    '''
    input:
        S4=ancient(rules.tableS4.output.table),
        S5=ancient(rules.tableS5.output.table)
    output:
        fig2a=protected('manuscript/figures/data/2a.csv'),
        fig2b=protected('manuscript/figures/data/2b.csv'),
        fig2c=protected('manuscript/figures/data/2c.csv'),
    log:
        'logs/20-Figures/F2.log'
    benchmark:
        'benchmarks/20-Figures/F2.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        footprints='results/12-BXD_footprints'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure 2a" > {log}
        Rscript workflow/scripts/20.Figure_2a.R -a {params.footprints}/plots/footprint_analysis.csv -b {params.footprints}/heatmap/differential_statistics.txt -o {params.dir}

        echo "Make Figure 2b" >> {log}
        Rscript workflow/scripts/20.Figure_2b.R -a {params.footprints}/plots/footprint_analysis.csv -o {params.dir}

        echo "Make Figure 2c" >> {log}
        Rscript workflow/scripts/20.Figure_2c.R -a {params.footprints}/plots/footprint_analysis.csv -b {input.S5} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule FigureS2:
    '''
    Figure S2
    '''
    input:
        S2=ancient(rules.tableS2.output.table),
        f=ancient(rules.Figure2.output.fig2a)
    output:
        figS2=protected('manuscript/figures/data/S2.csv'),
    log:
        'logs/20-Figures/FS2.log'
    benchmark:
        'benchmarks/20-Figures/FS2.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        counts='results/8-BXD_differential/DA/'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure S2" > {log}
        Rscript workflow/scripts/20.Figure_2S.R -a {params.counts}/atac_counts_disp.RData -b {params.counts}/rna_counts_disp.RData -c {input.S2} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule Figure3:
    '''
    Figure 3: Genetic influence on chromatin accessibility and transcriptome.
    '''
    input:
        S1=ancient(rules.tableS1.output.cre),
        S2=ancient(rules.tableS2.output.table),
        S6=ancient(rules.tableS6.output.table)
    output:
        fig3a=protected('manuscript/figures/data/3a.csv')
    log:
        'logs/20-Figures/F3.log'
    benchmark:
        'benchmarks/20-Figures/F3.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        qtls='results/9-BXD_qtl'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure 3a" > {log}
        Rscript workflow/scripts/20.Figure_3a.R -a {input.S6} -b {input.S1} -o {params.dir}

        echo "Make Figure 3b" >> {log}
        Rscript workflow/scripts/20.Figure_3b.R -a {params.qtls}/data/ -b {params.qtls}/QTL -c {input.S2} -d {input.S6} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule FigureS3:
    '''
    Figure S3
    '''
    input:
        S2=ancient(rules.tableS2.output.table),
        f=ancient(rules.Figure3.output.fig3a)
    output:
        figS3a=protected('manuscript/figures/data/S3a.csv'),
        figS3b=protected('manuscript/figures/data/S3b.csv'),
        figS3c=protected('manuscript/figures/data/S3c.csv'),
        figS3d=protected('manuscript/figures/data/S3d.csv'),
        figS3e=protected('manuscript/figures/data/S3e.csv'),
        figS3f=protected('manuscript/figures/data/S3f.csv')
    log:
        'logs/20-Figures/FS3.log'
    benchmark:
        'benchmarks/20-Figures/FS3.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        norm='results/7-BXD_normalization/EDA/'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure S3a" > {log}
        Rscript workflow/scripts/20.Figure_3Sa.R -a {params.norm}/atac_filtered_normalised_counts.RData -o {params.dir}

        echo "Make Figure S3b" >> {log}
        Rscript workflow/scripts/20.Figure_3Sb.R -a {params.norm}/atac_filtered_normalised_counts.RData -o {params.dir}

        echo "Make Figure S3c" >> {log}
        Rscript workflow/scripts/20.Figure_3Sc.R -a {params.norm}/rna_filtered_normalised_counts.RData -o {params.dir}

        echo "Make Figure S3d" >> {log}
        Rscript workflow/scripts/20.Figure_3Sd.R -a {params.norm}/rna_filtered_normalised_counts.RData -o {params.dir}

        echo "Make Figure S3e" >> {log}
        Rscript workflow/scripts/20.Figure_3Se.R -a {params.norm}/atac_filtered_normalised_counts.RData -b {input.S2} -o {params.dir}

        echo "Make Figure S3f" >> {log}
        Rscript workflow/scripts/20.Figure_3Sf.R -a {params.norm}/rna_filtered_normalised_counts.RData -b {input.S2} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule FigureS4:
    '''
    Figure S4
    '''
    input:
        S6=ancient(rules.tableS6.output.table),
        f=ancient(rules.Figure3.output.fig3a)
    output:
        figS3a=protected('manuscript/figures/data/S4a.csv'),
        figS3b=protected('manuscript/figures/data/S4b.csv')
    log:
        'logs/20-Figures/FS4.log'
    benchmark:
        'benchmarks/20-Figures/FS4.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/figures',
        qtls='results/9-BXD_qtl',
        counts='results/8-BXD_differential/DA/'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Figure S4a" > {log}
        Rscript workflow/scripts/20.Figure_4Sa.R -a {input.S6} -o {params.dir}

        echo "Make Figure S4b" >> {log}
        Rscript workflow/scripts/20.Figure_4Sb.R -a {params.qtls}/data/ -b {params.counts}/atac_counts_disp.RData -c {params.counts}/rna_counts_disp.RData  -d {input.S6} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''
