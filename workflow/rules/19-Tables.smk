'''
This Snakefile contains rules to produce the tables for the manuscript
'''
localrules: tableS1, tableS2, tableS3, tableS4, tableS5, tableS6, tableS7, tableS8, tableS9

rule tableS1:
    '''
    BXD ATAC-seq QC summary
    '''
    input:
        check_point=ancient(rules.BXD_checkpoint.output.check),
    output:
        table=protected('manuscript/tables/TableS1-ATAC_QC_summary.xlsx'),
        cre=protected('manuscript/tables/data/cCRE.RData'),
    log:
        'logs/19-Tables/S1.log'
    benchmark:
        'benchmarks/19-Tables/S1.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S1" > {log}
        Rscript workflow/scripts/19.1-TableS1.R  \
            -a results/1-BXD_fastq/summaries/QC_fastq.csv \
            -b results/3-BXD_bam/alignQC/QC_align.csv \
            -c results/3-BXD_bam/alignQC/QCbam_align.csv \
            -d logs/5-BXD_peaks/ \
            -e results/6-BXD_features//QC/QC_peaks.csv \
            -f results/6-BXD_features/list//genomic_features.bed \
            -g results/6-BXD_features/list//genomic_features.gtf \
            -i results/10-BXD_annotation/cre/mm10-cCREs.bed \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS2:
    '''
    BXD Differential analyses and correlations
    '''
    input:
        S1=ancient(rules.tableS1.output.table)
    output:
        table=protected('manuscript/tables/TableS2-Differential_Correlation.xlsx'),
    log:
        'logs/19-Tables/S2.log'
    benchmark:
        'benchmarks/19-Tables/S2.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S2" > {log}
        Rscript workflow/scripts/19.2-TableS2.R  \
            -a results/8-BXD_differential/DA/atac_stats_treament.RData \
            -b results/8-BXD_differential/DA/rna_stats_treament.RData \
            -c results/13-BXD_integrate/cors/gene_region_cors_annotation.csv \
            -d {input.S1} \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS3:
    '''
    BXD differential GO and KEGG enrichment analyses
    '''
    input:
        S2=ancient(rules.tableS2.output.table)
    output:
        table=protected('manuscript/tables/TableS3-Enrichment_BXD.xlsx'),
        enrich=protected('manuscript/tables/data/enrichment.RData'),
    log:
        'logs/19-Tables/S3.log'
    benchmark:
        'benchmarks/19-Tables/S3.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S3" > {log}
        Rscript workflow/scripts/19.3-TableS3.R  \
            -a {input.S2} \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS4:
    '''
    BXD footprints analyses
    '''
    input:
        check_point=ancient(rules.BXD_checkpoint.output.check)
    output:
        table=protected('manuscript/tables/TableS4-Footprints.xlsx'),
    log:
        'logs/19-Tables/S4.log'
    benchmark:
        'benchmarks/19-Tables/S4.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S4" > {log}
        Rscript workflow/scripts/19.4-TableS4.R  \
            -a results/12-BXD_footprints/motifs/ \
            -b results/12-BXD_footprints/plots/footprint_analysis.csv \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS5:
    '''
    KEA3: Improved kinase enrichment analysis via data integration 
    reference: 10.1093/nar/gkab359 
    '''
    input:
        check_point=ancient(rules.BXD_checkpoint.output.check)
    output:
        table=protected('manuscript/tables/TableS5-KEA_enrichment.xlsx'),
    log:
        'logs/19-Tables/S5.log'
    benchmark:
        'benchmarks/19-Tables/S5.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S5" > {log}
        Rscript workflow/scripts/19.5-TableS5.R  \
            -a results/12-BXD_footprints/KEA/ \
            -o {params.dir} 

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS6:
    '''
    BXD QTL results
    '''
    input:
        S2=ancient(rules.tableS2.output.table)
    output:
        table=protected('manuscript/tables/TableS6-QTL.xlsx'),
    log:
        'logs/19-Tables/S6.log'
    benchmark:
        'benchmarks/19-Tables/S6.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S6" > {log}
        Rscript workflow/scripts/19.6-TableS6.R  \
            -a results/9-BXD_qtl/tables/ \
            -b results/9-BXD_qtl/ttest/ \
            -c rawdata/sleep_bxd \
            -d {input.S2} \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS7:
    '''
    Sleep Deprivation response network
    '''
    input:
        S2=ancient(rules.tableS2.output.table),
        S4=ancient(rules.tableS4.output.table),
        S5=ancient(rules.tableS5.output.table),
        S6=ancient(rules.tableS6.output.table)
    output:
        table=protected('manuscript/tables/TableS7-GRN.xlsx'),
    log:
        'logs/19-Tables/S7.log'
    benchmark:
        'benchmarks/19-Tables/S7.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S7" > {log}
        Rscript workflow/scripts/19.7-TableS7.R  \
            -a results/13-BXD_integrate/granie/graph.RData \
            -b results/8-BXD_differential/DA/atac_counts_disp.RData \
            -c results/8-BXD_differential/DA/rna_counts_disp.RData \
            -d results/12-BXD_footprints/motifs/ \
            -e {input.S2} \
            -f {input.S4} \
            -g {input.S5} \
            -i {input.S6} \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS8:
    '''
    Nrf1 RNA-seq QC summary
    '''
    input:
        check_point=ancient(rules.NRF_checkpoint.output.check),
    output:
        table=protected('manuscript/tables/TableS8-RNAQC_summary.xlsx'),
    log:
        'logs/19-Tables/S8.log'
    benchmark:
        'benchmarks/19-Tables/S8.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S8" > {log}
        Rscript workflow/scripts/19.8-TableS8.R  \
            -a rawdata/nrf1/nrf_metadata.csv \
            -b results/14-NRF_fastq/summaries/QC_fastq.csv  \
            -c results/16-NRF_bam/alignQC/QC_align.csv \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule tableS9:
    '''
    Nrf1 Differential analysis
    '''
    input:
        check_point=ancient(rules.NRF_checkpoint.output.check),
    output:
        table=protected('manuscript/tables/TableS9-Differential_Expression.xlsx'),
    log:
        'logs/19-Tables/S9.log'
    benchmark:
        'benchmarks/19-Tables/S9.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='manuscript/tables'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}/data

        echo "Make Table S9" > {log}
        Rscript workflow/scripts/19.9-TableS9.R  \
            -a results/17-NRF_differential/DA/nrf_stats.RData \
            -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''
