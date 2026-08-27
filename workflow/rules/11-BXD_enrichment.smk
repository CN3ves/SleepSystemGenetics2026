'''
This Snakefile contains rules to perform enrichment analises:

enrichment_atac: Enrichment for ATAC regions
enrichment_rna: Enrichment for transcripts

'''

localrules: enrichment_atac, enrichment_rna

rule enrichment_atac:
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
        annotation = ancient(rules.annot2encode.output.annotation),
        features =  ancient(rules.make_features.output.gtf),
        stats =ancient('results/8-BXD_differential/DA/atac_stats_treament.RData'),
        check=ancient(rules.plot_volcanos.output.plot)
    output:
        annotation = protected('results/11-BXD_enrichment/atac/QTLxSD_ORA_GO_genes_all_enrichment.csv')
    log:
        'logs/11-BXD_enrichment/atac.log'
    benchmark:
        'benchmarks/11-BXD_enrichment/atac.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/11-BXD_enrichment/atac'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Running GO enrichment for ATAC regions" > {log}
        Rscript workflow/scripts/11.1-Enrichment_atac.R -a {input.annotation} -t GO -s {input.stats} -p {input.features} -o {params.dir} 2>> {log}

        echo "Running KEGG enrichment for ATAC regions" >> {log}
        Rscript workflow/scripts/11.1-Enrichment_atac.R -a {input.annotation} -t KEGG -s {input.stats} -p {input.features} -o {params.dir} 2>> {log}

        echo "Running Wikipaths enrichment for ATAC regions" >> {log}
        Rscript workflow/scripts/11.1-Enrichment_atac.R -a {input.annotation} -t Wikipaths -s {input.stats} -p {input.features} -o {params.dir} 2>> {log}

        echo "Running Reactome enrichment for ATAC regions" >> {log}
        Rscript workflow/scripts/11.1-Enrichment_atac.R -a {input.annotation} -t Reactome -s {input.stats} -p {input.features} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule enrichment_rna:
    '''
    Same as enrichment_atac, but using the RNA-seq results
    '''
    input:
        deg = ancient('results/8-BXD_differential/DA/rna_stats_treament.RData'),
        qtl = ancient(rules.qtl_aggregate.output.rna),
        int =  ancient(rules.qtl_aggregate.output.rnaFC),
        check=ancient(rules.plot_volcanos.output.plot)
    output:
        annotation = protected('results/11-BXD_enrichment/rna/QTLxSD_ORA_GO_genes_all_enrichment.csv')
    log:
        'logs/11-BXD_enrichment/rna.log'
    benchmark:
        'benchmarks/11-BXD_enrichment/rna.txt'
    resources:
        mem_mb = 500,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/11-BXD_enrichment/rna'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Running GO enrichment for RNA transcripts" > {log}
        Rscript workflow/scripts/11.2-Enrichment_rna.R -d {input.deg} -t GO -q {input.qtl} -i {input.int} -o {params.dir} 2>> {log}

        echo "Running KEGG enrichment for RNA transcripts" >> {log}
        Rscript workflow/scripts/11.2-Enrichment_rna.R -d {input.deg} -t KEGG -q {input.qtl} -i {input.int} -o {params.dir} 2>> {log}

        echo "Running Wikipaths enrichment for RNA transcripts" >> {log}
        Rscript workflow/scripts/11.2-Enrichment_rna.R -d {input.deg} -t Wikipaths -q {input.qtl} -i {input.int} -o {params.dir} 2>> {log}

        echo "Running Reactome enrichment for RNA transcripts" >> {log}
        Rscript workflow/scripts/11.2-Enrichment_rna.R -d {input.deg} -t Reactome -q {input.qtl} -i {input.int} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''