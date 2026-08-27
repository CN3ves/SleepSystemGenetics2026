'''
This Snakefile contains rules to perform statistical differential analyses:

differential_analysis_nrf: Estimate differential expression 
enrichment_rna_nrf: Enrichment for transcripts
'''

localrules: differential_analysis, enrichment_rna_nrf

rule differential_analysis_nrf:
    '''
    edgeR: Bioconductor package for differential expression analysis of digital gene expression data
    Reference: 10.1093/bioinformatics/btp616

    Performs:
    - Differential analyses for effect of SD using BXD strain as covariate (~0 + lines + treatment) using C57xCTRL as baseline.
    - Fits quasi-likelihood (QL) generalization of the negative binomial (NB) model
    - Approximate conditional likelihood is used to estimate feature dispertion using Cox-Reid profile-adjusted likelihood (CR)
    '''
    input:
        counts=ancient(rules.eda_nrf.output.counts),
        check=ancient(rules.multi_QC2.output.report)
    output:
        stats=protected('results/17-NRF_differential/DA/nrf_stats.RData'),
    log:
        'logs/17-NRF_differential/DA.log'
    benchmark:
        'benchmarks/17-NRF_differential/DA.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/17-NRF_differential/DA'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Differential analysis on <{input.counts}>" > {log}
        Rscript workflow/scripts/17.1-Differential_analysis.R -c '{input.counts}' -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule enrichment_rna_nrf:
    '''
    clusterProfiler: R package for Biomedical Knowledge Mining
    Reference: https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

    Performs:
    - Enrichment analysis using the Gene Ontology knowledge base
    - Enrichment analysis using the KEGG knowledge base
    - Enrichment analysis using the Wikipaths knowledge base
    - Enrichment analysis using the Reactome knowledge base
    '''
    input:
        deg = ancient(rules.differential_analysis_nrf.output.stats)
    output:
        enrichment = protected('results/17-NRF_differential/enrichment/Reactome_GSEA_CT_enrichment.csv')
    log:
        'logs/17-NRF_differential/enrichment.log'
    benchmark:
        'benchmarks/17-NRF_differential/enrichment.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/17-NRF_differential/enrichment'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Running GO enrichment for RNA transcripts" > {log}
        Rscript workflow/scripts/17.2-Enrichment.R -d {input.deg} -t GO -o {params.dir} 2>> {log}

        echo "Running KEGG enrichment for RNA transcripts" >> {log}
        Rscript workflow/scripts/17.2-Enrichment.R -d {input.deg} -t KEGG -o {params.dir} 2>> {log}

        echo "Running Wikipaths enrichment for RNA transcripts" >> {log}
        Rscript workflow/scripts/17.2-Enrichment.R -d {input.deg} -t Wikipaths -o {params.dir} 2>> {log}

        echo "Running Reactome enrichment for RNA transcripts" >> {log}
        Rscript workflow/scripts/17.2-Enrichment.R -d {input.deg} -t Reactome -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''