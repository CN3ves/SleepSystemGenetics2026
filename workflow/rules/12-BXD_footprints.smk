'''
This Snakefile contains rules to perform footprint analyses

merge_bams: Merges replicate reads into one high depth sample
scan_footprints: Scan reads for footprings
scan_motifs: Scan footprints for binding motifs
differential_activity: Differential footprint analyses for each line (HINT-ATAC produces p-values)
differential_heatmap: Simultaneous differential footprint analyses for all lines (HINT-ATAC normalizes equaly)
footprint_plots: Plot the results from the differential analyses
footprint_model: Run QTL model on footprints
footprint_perms: Runs permitation analyses
footprint_qtl: Run all QTL analyses for footprints
'''
localrules: scan_motifs, differential_heatmap, footprint_plots, footprint_qtl

rule merge_bams:
    '''
    Merge files from all BXD line and SD treatment groups into a single, high seqeuncing deapth sample
    (including downsampled files)
    '''
    input:
        bams = ancient(get_bams)
    output:
        full = protected('results/12-BXD_footprints/bams/{group}.bam'),
        sub1 = protected('results/12-BXD_footprints/bams/{group}_sub1.bam'),
        sub2 = protected('results/12-BXD_footprints/bams/{group}_sub2.bam')
    log:
        'logs/12-BXD_footprints/bams_{group}.log'
    benchmark:
        'benchmarks/12-BXD_footprints/bams_{group}.txt'
    resources:
        mem_mb = 10000,
        time = '1:00:00'
    threads: 10
    params:
        dir='results/12-BXD_footprints/bams'
    shell:
        '''
        module load  samtools/1.21
        mkdir -p {params.dir}

        echo "Merge bam files {input.bams}" > {log}
        samtools merge {params.dir}/{wildcards.group}.bam {input.bams} &>> {log}
        
        echo "Indexing file" >> {log}
        samtools index -@ 10 {params.dir}/{wildcards.group}.bam  &>> {log}

        echo "Checking downsampled files" >> {log}
        sub1=$(echo {input.bams} | sed 's/full/sub1/g' | tr ' ' '\n' | shuf --random-source=<(yes 1) -n 3 | tr '\n' ' ')
        sub2=$(echo {input.bams} | sed 's/full/sub2/g' | tr ' ' '\n' | shuf --random-source=<(yes 2) -n 3 | tr '\n' ' ')

        echo "Merge bam files $sub1" >> {log}
        samtools merge {params.dir}/{wildcards.group}_sub1.bam $sub1 &>> {log}
        samtools index -@ 10 {params.dir}/{wildcards.group}_sub1.bam &>> {log}

        echo "Merge bam files $sub2" >> {log}
        samtools merge {params.dir}/{wildcards.group}_sub2.bam $sub2 &>> {log}
        samtools index -@ 10 {params.dir}/{wildcards.group}_sub2.bam &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule scan_footprints:
    '''
    HINT-ATAC: A framework that uses open chromatin data to identify the active transcription factor binding sites.
    Reference: https://doi.org/10.1186/s13059-019-1642-2; https://reg-gen.readthedocs.io/en/latest/hint/tutorial-single-cell.html

    Performs:
    - Scans BAM files to identify footprings
    - Produces wig files of the identified motifs
    '''
    input:
        full = ancient(rules.merge_bams.output.full),
        sub1 = ancient(rules.merge_bams.output.sub1),
        sub2 = ancient(rules.merge_bams.output.sub2),
        bed = ancient(lambda wildcards: f'{rules.merge_regions.output.dar}'.replace(".RData",".bed"))
    output:
        bed = protected('results/12-BXD_footprints/footprints/{group}.bed'),
        sub1 = protected('results/12-BXD_footprints/footprints/{group}_sub1.bed'),
        sub2 = protected('results/12-BXD_footprints/footprints/{group}_sub2.bed')
    log:
        'logs/12-BXD_footprints/footprints_{group}.log'
    benchmark:
        'benchmarks/12-BXD_footprints/footprints_{group}.txt'
    resources:
        mem_mb = 1000,
        time = '0:40:00'
    conda:
        config['hint_env']
    threads: 20
    params:
        dir='results/12-BXD_footprints/footprints'
    shell:
        '''
        mkdir -p {params.dir}

        echo "Calling footprint for {input.full}" > {log}
        rgt-hint footprinting --atac-seq  --organism=mm10 --output-location={params.dir} --output-prefix={wildcards.group} {input.full} {input.bed} &>> {log}

        echo "Make bigwigs" >> {log}
        rgt-hint tracks --bc  --organism=mm10 --output-location={params.dir} --output-prefix={wildcards.group} {input.full} {input.bed} &>> {log}

        echo "Calling footprint for {input.sub1}" >> {log}
        rgt-hint footprinting --atac-seq  --organism=mm10 --output-location={params.dir} --output-prefix={wildcards.group}_sub1 {input.sub1} {input.bed} &>> {log}

        echo "Calling footprint for {input.sub2}" >> {log}
        rgt-hint footprinting --atac-seq  --organism=mm10 --output-location={params.dir} --output-prefix={wildcards.group}_sub2 {input.sub2} {input.bed} &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule scan_motifs:
    '''
    HINT-ATAC: A framework that uses open chromatin data to identify the active transcription factor binding sites.
    Reference: https://doi.org/10.1186/s13059-019-1642-2; https://reg-gen.readthedocs.io/en/latest/hint/tutorial-single-cell.html

    Performs:
    - Find motif predicted binding sites for motifs using the default JASPAR database
    '''
    input:
        beds = ancient(get_footprints)
    output:
        bed = protected('results/12-BXD_footprints/motifs/DBA_CTRL_mpbs.bed'),
    log:
        'logs/12-BXD_footprints/motifs.log'
    benchmark:
        'benchmarks/12-BXD_footprints/motifs.txt'
    resources:
        mem_mb = 10000,
        time = '0:10:00'
    conda:
        config['hint_env']
    threads: 1
    params:
        dir='results/12-BXD_footprints/motifs'
    shell:
        '''
        mkdir -p {params.dir}
        
        echo "Scan motifs" > {log}
        rgt-motifanalysis matching --organism=mm10 --input-files {input.beds} --output-location {params.dir} &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule differential_activity:
    '''
    HINT-ATAC: A framework that uses open chromatin data to identify the active transcription factor binding sites.
    Reference: https://doi.org/10.1186/s13059-019-1642-2; https://reg-gen.readthedocs.io/en/latest/hint/tutorial-dendritic-cell.html

    Performs:
    - Estimate differential activity between treated and untreads samples for each BXD line
    '''
    input:
        beds = ancient(lambda wildcards: [f'results/12-BXD_footprints/motifs/{wildcards.group}_{sd}{wildcards.sub}_mpbs.bed' for sd in ['CTRL','SD']]),
        bams = ancient(lambda wildcards: [f'results/12-BXD_footprints/bams/{wildcards.group}_{sd}{wildcards.sub}.bam' for sd in ['CTRL','SD']])
    output:
        diff = protected('results/12-BXD_footprints/diff/{group, [a-zA-Z0-9]*}{sub, (|_sub1|_sub2)}/differential_statistics.txt'),
    log:
        'logs/12-BXD_footprints/diff_{group, [a-zA-Z0-9]*}{sub, (_sub1|_sub2)*}.log'
    benchmark:
        'benchmarks/12-BXD_footprints/diff_{group, [a-zA-Z0-9]*}{sub, (_sub1|_sub2)*}.txt'
    resources:
        mem_mb = 15000,
        time = '0:25:00'
    conda:
        config['hint_env']
    threads: 20
    params:
        dir='results/12-BXD_footprints/diff'
    shell:
        '''
        mkdir -p {params.dir}
        
        echo "Estimate differential activity for {wildcards.group} {wildcards.sub}" > {log}
        cond=$(echo {input.beds} | sort | sed 's/results\\/12-BXD_footprints\\/motifs\\///g' | sed 's/_mpbs.bed//g' | sed 's/_sub[12]//g' | sed  's/.*_\\(.*\\) .*_\\(.*\\)/\\1,\\2/')
        mpbs=$(echo {input.beds} | sort | tr ' ' ',' | sed 's/,$//g')
        bams=$(echo {input.bams} | sort | tr ' ' ',' | sed 's/,$//g')

        echo "BED files: $mpbs" >> {log}
        echo "BAM files: $bams" >> {log}
        echo "groups: $cond" >> {log}

        rgt-hint differential --organism=mm10 --bc --nc 20 --mpbs-files=$mpbs --reads-files=$bams --conditions=$cond --output-location={params.dir}/{wildcards.group}{wildcards.sub} &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule differential_heatmap:
    '''
    Same as differential_activity, but considering all samples together.
    Note: [total time: 85h 30m 26s]
    '''
    input:
        beds = ancient(all_bed),
        bams = ancient(all_bam)
    output:
        diff = protected('results/12-BXD_footprints/heatmap/differential_statistics.txt'),
    log:
        'logs/12-BXD_footprints/diff_all.log'
    benchmark:
        'benchmarks/12-BXD_footprints/diff_all.txt'
    resources:
        mem_mb = 50000,
        time = '90:00:00'
    conda:
        config['hint_env']
    threads: 20
    params:
        dir='results/12-BXD_footprints/heatmap'
    shell:
        '''
        mkdir -p {params.dir}
        
        echo "Estimate differential activity for all samples" > {log}
        cond=$(echo {input.beds} | sort | sed 's/_mpbs.bed//g' | sed 's/results\\/12-BXD_footprints\\/motifs\\///g' |  tr ' ' ',' | sed 's/,$//g')
        mpbs=$(echo {input.beds} | sort | tr ' ' ',' | sed 's/,$//g')
        bams=$(echo {input.bams} | sort | tr ' ' ',' | sed 's/,$//g')

        echo "BED files: $mpbs" >> {log}
        echo "BAM files: $bams" >> {log}
        echo "groups: $cond" >> {log}

        rgt-hint differential --organism=mm10 --bc --nc 20 --mpbs-files=$mpbs --reads-files=$bams --conditions=$cond --output-location={params.dir}/ --no-lineplots &>> {log} 
        echo "Logs saved in <{log}>" >> {log}
        '''

rule footprint_plots:
    '''
    Plot the results from the differential footprint analyses

    Performs:
    - Fisher's combined probability test to aggregation statistics for all BXD lines 
    - Plot heatmap with all estimated activity changes 
    '''
    input:
        diffs=diffprints,
        heatmap=ancient(rules.differential_heatmap.output.diff)
    output:
        plot=protected('results/12-BXD_footprints/plots/agreggated_diff{sub, (_sub1|_sub2)*}.svg'),
    log:
        'logs/12-BXD_footprints/plots{sub, (_sub1|_sub2)*}.log'
    benchmark:
        'benchmarks/12-BXD_footprints/plots{sub, (_sub1|_sub2)*}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/12-BXD_footprints/plots'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Ploting results for footprint analyses on <{wildcards.sub}>" > {log}
        Rscript workflow/scripts/12.1-Footprint_plot.R -f '{input.diffs}' -a {input.heatmap} -t '{wildcards.sub}' -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule footprint_model:
    '''
    QTL for footprints
    '''
    input:
        check=ancient(lambda wildcards: [f'results/12-BXD_footprints/plots/agreggated_diff{sub}.svg'for sub in ['','_sub1','_sub2']]),   
    output:
        qtl=protected('results/12-BXD_footprints/model/footprints_QTL_model.RData'),
        int=protected('results/12-BXD_footprints/model/footprintsFC_QTL_model.RData'),
    log:
        'logs/12-BXD_footprints/model.log'
    benchmark:
        'benchmarks/12-BXD_footprints/model.txt'
    resources:
        mem_mb = 20000,
        time = '02:00:00'
    threads: 40
    params:
        dir='results/12-BXD_footprints/model',
        source=rules.get_info.params.dir,
        counts='results/12-BXD_footprints/plots/footprints_counts_disp.RData'
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Prepare <{params.counts}> for QTL analysis" > {log}
        Rscript workflow/scripts/9.2_QTL_prepare.R -c '{params.counts}' -d {params.source} -o {params.dir} 2>> {log}
        cp {params.source}/*map* {params.dir}

        echo "Running model on footprints.json" >> {log}
        Rscript workflow/scripts/9.3_QTL_run.R -j {params.dir}/footprints.json -o {params.dir} 2>> {log}
        echo "Running model on footprintsFC.json" >> {log}
        Rscript workflow/scripts/9.3_QTL_run.R -j {params.dir}/footprintsFC.json -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule footprint_perms:
    '''
    Permutation for footprints
    '''
    input:
        qtl=ancient(rules.footprint_model.output.qtl),
        int=ancient(rules.footprint_model.output.int),
    output:
        qtl=protected('results/12-BXD_footprints/perms/footprints_QTL_perm{i}.RData'),
        inf=protected('results/12-BXD_footprints/perms/footprintsFC_QTL_perm{i}.RData'),
    log:
        'logs/12-BXD_footprints/QTL_footprints_perm_{i}.log'
    benchmark:
        'benchmarks/12-BXD_footprints/QTL_footprints_perm_{i}.txt'
    resources:
        mem_mb = 100000,
        time = '4:00:00'
    threads: 10
    params:
        dir='results/12-BXD_footprints/perms',
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "QTL permutations number {wildcards.i}" > {log}
        echo "<{input.qtl}>" >> {log}
        Rscript workflow/scripts/9.4_QTL_perms.R -q {input.qtl} -i {wildcards.i} -o {params.dir} 2>> {log}

        echo "<{input.int}>" >> {log}
        Rscript workflow/scripts/9.4_QTL_perms.R -q {input.int} -i {wildcards.i} -o {params.dir} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule footprint_qtl:
    '''
    Permutation for footprints
    '''
    input:
        qtl=ancient(rules.footprint_model.output.qtl),
        int=ancient(rules.footprint_model.output.int),
        perms_qtl=ancient(lambda wildcards: [f"results/12-BXD_footprints/perms/footprints_QTL_perm{i}.RData" for i in range(1, int(config['qtl_permutation_split'])+1)]),
        perms_int=ancient(lambda wildcards: [f"results/12-BXD_footprints/perms/footprintsFC_QTL_perm{i}.RData" for i in range(1, int(config['qtl_permutation_split'])+1)]),
    output:
        qtls=protected('results/9-BXD_qtl/ttest/qtl_ttests_footprints.csv'),
    log:
        'logs/12-BXD_footprints/qtl.log'
    benchmark:
        'benchmarks/12-BXD_footprints/qtl.txt'
    resources:
        mem_mb = 200000,
        time = '3:00:00'
    threads: 10
    params:
        dir='results/12-BXD_footprints/qtl',
        model=rules.footprint_model.params.dir,
        tab=rules.qtl_aggregate.params.dir,
        inv=rules.qtl_ttest.params.dir
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}

        echo "Searching significant QTL peaks for <{input.qtl}>" > {log}
        echo "Permutation files <{input.perms_qtl}>" >> {log}
        Rscript workflow/scripts/9.5_QTL_sig.R -m {input.qtl} -p '{input.perms_qtl}' -o {params.dir} 2>> {log}

        echo "Searching significant QTL peaks for <{input.int}>" >> {log}
        echo "Permutation file <{input.perms_int}>" >> {log}
        Rscript workflow/scripts/9.5_QTL_sig.R -m {input.int} -p '{input.perms_int}' -o {params.dir} 2>> {log}

        echo "Aggregating QTL results into a single table" >> {log}
        Rscript workflow/scripts/12.2_QTL_aggregate.R -t '{params.dir}/footprints_QTL_sigs.csv {params.dir}/footprintsFC_QTL_sigs.csv' -o {params.dir} 2>> {log}

        echo "Classifying interaction QTL results from <{input.int}>" >> {log}
        Rscript workflow/scripts/9.7_QTL_classify.R -q {params.model} -t results/12-BXD_footprints/qtl/QTL_table_filter_footprintsFC.csv -o {params.dir} 2>> {log}

        echo "Adding QTL information to <{params.tab}>" >> {log}
        cp {params.dir}/QTL_table* {params.tab}/
        echo "Adding interaction information to <{params.inv}>" >> {log}
        cp {params.dir}/qtl_ttest* {params.inv}/

        echo "Logs saved in <{log}>" >> {log}
        '''
