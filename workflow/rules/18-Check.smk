'''
This Snakefile contains rules to check and enforce all missing rules to run

BXD_checkpoint: Gather all BXD results for figures
NRF_checkpoint: Gather all Nrf results for figures
'''
localrules: BXD_checkpoint, NRF_checkpoint

rule BXD_checkpoint:
    '''
    Enforces all "loose-ends" (rules not pulled downstream by another rule) of BXD integrative analyses to run:
    - Differential analyses plots
    - Enrichment analyses plots
    - Footprint analyses
    - Correlation analyses
    - Regulatory network inference 

    Then lists as outputs all the files required to produce the manuscript figures
    '''
    input:
        check_cors=ancient('results/13-BXD_integrate/plots/cors_int.png'),
        check_grn=ancient('results/13-BXD_integrate/granie/graph.RData'),
        check_diff=ancient(lambda wildcards: [f"results/8-BXD_differential/plots/{exp}_coeff_heatFC.png" for exp in ['rna','atac']]),
        check_enrich=ancient(lambda wildcards: [f'results/11-BXD_enrichment/{exp}/QTLxSD_ORA_GO_genes_all_enrichment.csv' for exp in ['rna','atac']])
    output:
        check=protected('results/18-Checkpoint/bxd.check'),
    log:
        'logs/18-Checkpoint/bxd.log'
    benchmark:
        'benchmarks/18-Checkpoint/bxd.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='results/18-Checkpoint/'
    shell:
        '''
        mkdir -p {params.dir}

        echo "List of all outputs for BXD analyses:" > {log}
        find results/ -type f | grep BXD | sed 's/results\\///g' | sort -n >> {log}
        touch {params.dir}/bxd.check
        '''

rule NRF_checkpoint:
    '''
    Enforces all "loose-ends" (rules not pulled downstream by another rule) of NRF integrative analyses to run:
    - Enrichment analyses plots

    Then lists as outputs all the files required to produce the manuscript figures
    '''
    input:
        check_enrich=ancient('results/17-NRF_differential/enrichment/Reactome_GSEA_CT_enrichment.csv'),
    output:
        check=protected('results/18-Checkpoint/nrf.check'),
    log:
        'logs/18-Checkpoint/bxd.log'
    benchmark:
        'benchmarks/18-Checkpoint/bxd.txt'
    resources:
        mem_mb = 1000,
        time = '0:10:00'
    threads: 1
    params:
        dir='results/18-Checkpoint/'
    shell:
        '''
        mkdir -p {params.dir}

        echo "List of all outputs for NRF analyses:" > {log}
        find results/ -type f | grep NRF | sed 's/results\\///g' | sort -n >> {log}
        touch {params.dir}/nrf.check
        '''

       
       
