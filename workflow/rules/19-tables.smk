'''
This Snakefile contains rules to produce the tables for the manuscript
'''
localrules: tableS1

rule tableS1:
    '''
    
    '''
    input:
        check_cors=ancient('results/13-BXD_integrate/plots/cors_int.png'),
        check_grn=ancient('results/13-BXD_integrate/granie/graph.RData'),
        check_diff=ancient(lambda wildcards: [f"results/8-BXD_differential/plots/{exp}_coeff_heatFC.png" for exp in ['rna','atac']]),
        check_enrich=ancient(lambda wildcards: [f'results/11-BXD_enrichment/{exp}/QTLxSD_ORA_GO_genes_all_enrichment.csv' for exp in ['rna','atac']]),
        check_footprints=ancient(lambda wildcards: [f'results/12-BXD_footprints/plots/agreggated_diff{sub}.svg' for sub in ['','_sub1','_sub2']])
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


       
       
