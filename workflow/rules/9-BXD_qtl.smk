'''
This Snakefile contains rules to perform QTL analyses:

get_info: Download BXD information fro QTL analysos
sleep_aspects: Prepares sleep phenotype data for use
prepare_data: Prepares all required files for QTL analysis
qtl_model: Run QTL analysis
qtl_perms: Runs permitation analyses
qtl_peaks: Calls significant LOD peaks
qtl_aggregate: Merges the results split due to numebr of permutations into a individual tables
qtl_ttest: Rund post-hoc tests to identify baseline and SD variant-dependent response 
'''

localrules: get_info, sleep_aspects, prepare_data, qtl_aggregate

rule get_info:
    '''
    Download [BXD information](https://github.com/rqtl/qtl2data) 
    previously processed from [GeneNetworks](https://github.com/kbroman/Teaching_CTC2019):
    '''
    output:
        check  = 'results/9-BXD_qtl/info/bxd.json',
    log:
        'logs/9-BXD_qtl/info.log'
    benchmark:
        'benchmarks/9-BXD_qtl/info.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/9-BXD_qtl/info'
    shell:
        '''
        mkdir -p {params.dir}

        echo "Downloading BXD line information" >> {log}
        wget https://raw.githubusercontent.com/rqtl/qtl2data/main/BXD/bxd.zip -P {params.dir}/
        unzip {params.dir}/bxd.zip -d {params.dir}/
        rm {params.dir}/bxd.zip 

        echo "Logs saved in <{log}>" >> {log}
        '''

rule sleep_aspects:
    '''
    Converts sleep phenotyping from Diessler, 2018 to edgeR object.
    Reference: 10.1371/journal.pbio.2005750 
    '''
    input:
        counts='rawdata/sleep_bxd/phenotypes.txt',
        info=rules.get_info.output.check
    output:
        pheno=protected('results/9-BXD_qtl/aspects/sleep_counts_disp.RData'),
        rna=protected('results/9-BXD_qtl/aspects/rna_counts_disp.RData'),
        atac=protected('results/9-BXD_qtl/aspects/atac_counts_disp.RData'),
    log:
        'logs/9-BXD_qtl/aspects.log'
    benchmark:
        'benchmarks/9-BXD_qtl/aspects.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/9-BXD_qtl/aspects',
        data=rules.differential_analysis.params.dir
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Loading sleep data <{input.counts}> into R object" > {log}
        Rscript workflow/scripts/9.1_QTL_sleep.R -c '{input.counts}'  -o {params.dir} 

        echo "Copying other results to <{params.dir}>" >> {log}
        echo "cp {params.data}/*disp* {params.dir}" >> {log}
        cp {params.data}/*disp* {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule prepare_data:
    '''
    Prepares the required files for QTL analysis and defines the configuration json file 
    '''
    input:
        counts=ancient(lambda wildcards: f"results/9-BXD_qtl/aspects/{wildcards.exp}_counts_disp.RData")
    output:
        data=protected('results/9-BXD_qtl/data/{exp}_geno.csv'),
    log:
        'logs/9-BXD_qtl/data_{exp}.log'
    benchmark:
        'benchmarks/9-BXD_qtl/data_{exp}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/9-BXD_qtl/data',
        source=rules.get_info.params.dir
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Prepare <{input.counts}> for QTL analysis" > {log}
        Rscript workflow/scripts/9.2_QTL_prepare.R -c '{input.counts}' -d {params.source} -o {params.dir} 2>> {log}

        cp {params.source}/*map* {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule qtl_model:
    '''
    r/qtl2: R package for mapping quantitative trait loci (QTL) in experimental crosses
    Reference: https://kbroman.org/qtl2/; https://doi.org/10.1534/genetics.118.301595

    Performs:
    - Main model with treatment as covariate (scan1)
    '''
    input:
        check=ancient('results/9-BXD_qtl/data/{exp}_geno.csv')
    output:
        data=protected('results/9-BXD_qtl/QTL/{exp, atac|rna|sleep}{chr, [chr]*[0-9XY]*}{fc, [FC]*}_QTL_model.RData'),
    log:
        'logs/9-BXD_qtl/QTL_{exp}_{chr}_{fc}_model.log'
    benchmark:
        'benchmarks/9-BXD_qtl/QTL_{exp}_{chr}_{fc}_model.txt'
    resources:
        mem_mb = 20000,
        time = '00:40:00'
    threads: 40
    params:
        dir='results/9-BXD_qtl/QTL',
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Prepare for QTL analysis" > {log}
        json="results/9-BXD_qtl/data/{wildcards.exp}{wildcards.chr}{wildcards.fc}.json"

        echo "Running model on <$json>" >> {log}
        Rscript workflow/scripts/9.3_QTL_run.R -j $json -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule qtl_perms:
    '''
    r/qtl2: R package for mapping quantitative trait loci (QTL) in experimental crosses
    Reference: https://kbroman.org/qtl2/; https://doi.org/10.1534/genetics.118.301595

    Performs:
    - Permutation analysis to generate null distribution (scan1perm)
    '''
    input:
        qtl=ancient('results/9-BXD_qtl/QTL/{exp}{chr}{fc}_QTL_model.RData')
    output:
        data=protected('results/9-BXD_qtl/QTL/{exp, atac|rna|sleep}{chr, [chr]*[0-9XY]*}{fc, [FC]*}_QTL_perm{i}.RData'),
    log:
        'logs/9-BXD_qtl/QTL_{exp}_{chr}_{fc}_perm_{i}.log'
    benchmark:
        'benchmarks/9-BXD_qtl/QTL_{exp}_{chr}_{fc}_perm_{i}.txt'
    resources:
        mem_mb = 170000,
        time = '6:15:00'
    threads: 10
    params:
        dir='results/9-BXD_qtl/QTL',
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "QTL permutations number {wildcards.i}" > {log}
        echo "<{input.qtl}>" >> {log}
        Rscript workflow/scripts/9.4_QTL_perms.R -q {input.qtl} -i {wildcards.i} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule qtl_peaks:
    '''
    r/qtl2: R package for mapping quantitative trait loci (QTL) in experimental crosses
    Reference: https://kbroman.org/qtl2/; https://doi.org/10.1534/genetics.118.301595

    Performs:
    - Significance calling of LOD peaks using different statistical cut-offs (find_peaks)
    '''
    input:
        model=ancient(rules.qtl_model.output.data),
        perms=ancient(get_permutations)
    output:
        qtl=protected('results/9-BXD_qtl/sigs/{exp, atac|rna|sleep}{chr, [chr]*[0-9XY]*}{fc, [FC]*}_QTL_sigs.csv'),
    log:
        'logs/9-BXD_qtl/sigs_{exp}_{chr}_{fc}.log'
    benchmark:
        'benchmarks/9-BXD_qtl/sigs_{exp}_{chr}_{fc}.txt'
    resources:
        mem_mb = 40000,
        time = '0:40:00'
    threads: 1
    params:
        dir='results/9-BXD_qtl/sigs',
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Searching significant QTL peaks for <{input.model}>" > {log}
        Rscript workflow/scripts/9.5_QTL_sig.R -m {input.model} -p '{input.perms}' -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule qtl_aggregate:
    '''
    Agregates all separate results into individual tables and filters results avois genetic imbalanced results
    '''
    input:
        qtls=ancient(get_qtls),
    output:
        atac=protected('results/9-BXD_qtl/tables/QTL_table_filter_atac.csv'),
        atacFC=protected('results/9-BXD_qtl/tables/QTL_table_filter_atacFC.csv'),
        rna=protected('results/9-BXD_qtl/tables/QTL_table_filter_rna.csv'),
        rnaFC=protected('results/9-BXD_qtl/tables/QTL_table_filter_rnaFC.csv'),
        sleep=protected('results/9-BXD_qtl/tables/QTL_table_filter_sleep.csv'),
    log:
        'logs/9-BXD_qtl/aggregate.log'
    benchmark:
        'benchmarks/9-BXD_qtl/aggregate.txt'
    resources:
        mem_mb = 500,
        time = '0:10:00'
    threads: 1
    params:
        dir='results/9-BXD_qtl/tables',
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Aggregating QTL results into a single table" > {log}
        echo "<{input.qtls}>" >> {log}
        Rscript workflow/scripts/9.6_QTL_aggregate.R -t '{input.qtls}' -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule qtl_ttest:
    '''
    Performs pot-hoc t-test to identify different trends between variants in basline and response to SD
    '''
    input:
        qtls=ancient(lambda wildcards: f"results/9-BXD_qtl/tables/QTL_table_filter_{wildcards.exp}FC.csv"),
    output:
        classification=protected('results/9-BXD_qtl/ttest/qtl_ttests_{exp}.csv'),
    log:
        'logs/9-BXD_qtl/ttest_{exp}FC.log'
    benchmark:
        'benchmarks/9-BXD_qtl/ttest_{exp}.txt'
    resources:
        mem_mb = 100000,
        time = '8:00:00'
    threads: 20
    params:
        dir='results/9-BXD_qtl/ttest',
        models=rules.qtl_model.params.dir
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Classifying interaction QTL results from <{input.qtls}>" > {log}
        Rscript workflow/scripts/9.7_QTL_classify.R -q {params.models} -t {input.qtls} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''