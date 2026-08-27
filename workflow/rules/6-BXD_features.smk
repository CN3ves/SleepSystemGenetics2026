'''
This Snakefile contains rules to unify genomic features and count reads matching to them:

make_features: Collect called peaks and make a unified list of features for analysis
make_mito (optional): Collect called peaks (including mitochondrial DNA) and make a unified list of features for analysis
make_parents (optional): Collect called peaks from sampled parent triplicate lines (mostly to close the pipeline)
make_subsample (optional): Collect called peaks from downsampled BAM files and make a unified list of features for analysis
get_htseq: Install HTSeq python package
count_features: Count reads matching the genomic feature list
frip (optional): Use the standard featurecount to estimat Fraction of Reads in Peaks 
count_QC: Run QC statistics for read count
multiQC: Run multiQC for summary of all analyses results so far
'''

localrules: make_features, make_mito,make_parents, make_subsample, get_htseq, count_QC, multi_QC

rule make_features:
    '''
    bedtools2: The bedtools utilities provide a wide-range of genomics analysis tasks, the most widely-used enable genome arithmetic (set theory on the genome)
    Reference:

    Performs:
    - Collect all narrowpeak files 
    - Merge all identified peaks into a candidate list
    '''
    input:
        ancient(get_peaks)
    output:
        bed=protected('results/6-BXD_features/list/genomic_features.bed'),
        gtf=protected('results/6-BXD_features/list/genomic_features.gtf'),
    log:
        'logs/6-BXD_features/features.log'
    benchmark:
        'benchmarks/6-BXD_features/features.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/6-BXD_features/list'
    shell:
        '''
        module load bedtools2/2.31.1 
        mkdir -p {params.dir}

        echo "Merging <{input}> narrowpeaks and converting to <merged_features.bed>" >> {log}
        zcat {input} | cut -f 1-6 >  {params.dir}/merged_features.bed 2>> {log}
        echo "Sorting features to <merged_features.sorted.bed>" >> {log}
        sort -k1,1 -k2,2n {params.dir}/merged_features.bed > {params.dir}/merged_features.sorted.bed 2>> {log}
        echo "Reducing features to <{output.bed}>" >> {log}
        bedtools merge -i  {params.dir}/merged_features.sorted.bed >  {output.bed} 2>> {log}

       
        python3  workflow/scripts/6.1_Peaks_bed2gft.py {output.bed} {output.gtf} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule make_mito:
    '''
    Same as make_features, but including mitochondrial DNA.
    '''
    input:
        ancient(get_mito)
    output:
        bed=protected('results/6-BXD_features/mito/genomic_features.bed'),
        gtf=protected('results/6-BXD_features/mito/genomic_features.gtf'),
    log:
        'logs/6-BXD_features/mito.log'
    benchmark:
        'benchmarks/6-BXD_features/mito.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/6-BXD_features/mito'
    shell:
        '''
        module load bedtools2/2.31.1 
        mkdir -p {params.dir}

        echo "Merging <{input}> narrowpeaks and converting to <merged_features.bed>" >> {log}
        zcat {input} | cut -f 1-6 >  {params.dir}/merged_features.bed 2>> {log}
        echo "Sorting features to <merged_features.sorted.bed>" >> {log}
        sort -k1,1 -k2,2n {params.dir}/merged_features.bed > {params.dir}/merged_features.sorted.bed 2>> {log}
        echo "Reducing features to <{output.bed}>" >> {log}
        bedtools merge -i  {params.dir}/merged_features.sorted.bed >  {output.bed} 2>> {log}

        echo "Spliting features into smaller (<250bp) fragments in <{output.gtf}> " >> {log}
        python3  workflow/scripts/6.1_Peaks_bed2gft.py {output.bed} {output.gtf} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule make_parents:
    '''
    Same as make_features, but for random parent samples.
    '''
    input:
        ancient(get_parents)
    output:
        bed=protected('results/6-BXD_features/parent/genomic_features{i}.bed'),
        gtf=protected('results/6-BXD_features/parent/genomic_features{i}.gtf'),
    log:
        'logs/6-BXD_features/parent_{i}.log'
    benchmark:
        'benchmarks/6-BXD_features/parent_{i}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/6-BXD_features/parent'
    shell:
        '''
        module load bedtools2/2.31.1 
        mkdir -p {params.dir}

        echo "Merging <{input}> narrowpeaks and converting to <merged_features{wildcards.i}.bed>" >> {log}
        zcat {input} | cut -f 1-6 >  {params.dir}/merged_features{wildcards.i}.bed 2>> {log}
        echo "Sorting features to <merged_features{wildcards.i}.sorted.bed>" >> {log}
        sort -k1,1 -k2,2n {params.dir}/merged_features{wildcards.i}.bed > {params.dir}/merged_features{wildcards.i}.sorted.bed 2>> {log}
        echo "Reducing features to <{output.bed}>" >> {log}
        bedtools merge -i  {params.dir}/merged_features{wildcards.i}.sorted.bed >  {output.bed} 2>> {log}

        echo "Spliting features into smaller (<250bp) fragments in <{output.gtf}>" >> {log}
        python3  workflow/scripts/6.1_Peaks_bed2gft.py {output.bed} {output.gtf} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule make_subsample:
    '''
    Same as make_features, but for subsamples reads.
    '''
    input:
        ancient(get_subsample)
    output:
        bed=protected('results/6-BXD_features/subsample/genomic_features{i}.bed'),
        gtf=protected('results/6-BXD_features/subsample/genomic_features{i}.gtf'),
    log:
        'logs/6-BXD_features/subsample_{i}.log'
    benchmark:
        'benchmarks/6-BXD_features/subsample_{i}.txt'
    resources:
        mem_mb = 5000,
        time = '00:10:00'
    threads: 1
    params:
        dir='results/6-BXD_features/subsample'
    shell:
        '''
        module load bedtools2/2.31.1 
        mkdir -p {params.dir}

        echo "Merging <{input}> narrowpeaks and converting to <merged_features{wildcards.i}.bed>" >> {log}
        zcat {input} | cut -f 1-6 >  {params.dir}/merged_features{wildcards.i}.bed 2>> {log}
        echo "Sorting features to <merged_features{wildcards.i}.sorted.bed>" >> {log}
        sort -k1,1 -k2,2n {params.dir}/merged_features{wildcards.i}.bed > {params.dir}/merged_features{wildcards.i}.sorted.bed 2>> {log}
        echo "Reducing features to <{output.bed}>" >> {log}
        bedtools merge -i  {params.dir}/merged_features{wildcards.i}.sorted.bed > {output.bed} 2>> {log}

        echo "Spliting features into smaller (<250bp) fragments in <{output.gtf}> " >> {log}
        python3  workflow/scripts/6.1_Peaks_bed2gft.py {output.bed} {output.gtf} &>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule get_htseq:
    '''
    Install Htseq python library to count features
    reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287950/
    '''
    output:
        check  = 'results/6-BXD_features/htseq/check',
    log:
        'logs/6-BXD_features/htseq.log'
    benchmark:
        'benchmarks/6-BXD_features/htseq.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/6-BXD_features/htseq'
    shell:
        '''
        mkdir -p {params.dir}

        echo "Installing HTSeq python library" >> {log}
        pip install HTSeq 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        touch {params.dir}/check

        '''

rule count_features:
    '''
    HTSeq: A Python package that provides infrastructure to process data from high-throughput sequencing assays.
    reference: https://htseq.readthedocs.io/en/release_0.11.1/count.html

    Peforms:
    - Counts reads matching features provided
    - Mode nonunique union 
    '''
    input:
        bam=ancient(rules.correct_samples.output.full_co),
        features=ancient(rules.make_features.output.gtf),
        check=ancient(rules.get_htseq.output.check),
    output:
        counts=protected('results/6-BXD_features/counts/{sample}.count'),
    log:
        'logs/6-BXD_features/counts_{sample}.log'
    benchmark:
        'benchmarks/6-BXD_features/counts_{sample}.txt'
    resources:
        mem_mb = 1000,
        time = '00:45:00'
    threads: 1
    params:
        dir='results/6-BXD_features/counts'
    shell:
        '''
        mkdir -p {params.dir}

        echo "Counting featured for <{input.bam}>" >> {log}
        python3 -m HTSeq.scripts.count -r pos -s no -m union -t Peak -i Peak_ID -q --nonunique all -f bam {input.bam} {input.features} > {output.counts} 2>> {log}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule frip:
    '''
    featureCounts: An efficient general purpose program for assigning sequence reads to genomic features
    Reference: https://doi.org/10.1093/bioinformatics/btt656

    Performs:
    - Count features matching genomic features (peaks)
    - Estimates Fraction of Reads in Peaks
    '''
    input:
        gtf=ancient(lambda wildcards: f"results/6-BXD_features/{wildcards.type}/genomic_features{wildcards.i}.gtf")
    output:
        tab="results/6-BXD_features/frip/frip_{type,[A-Za-z]+}{i, [0-9]*}.csv"
    log:
        'logs/6-BXD_features/frip_{type,[A-Za-z]+}{i, [0-9]*}.log'
    benchmark:
        'benchmarks/6-BXD_features/frip_{type,[A-Za-z]+}{i, [0-9]*}.txt'
    resources:
        mem_mb = 15000,
        time = '01:00:00'
    threads: 30
    params:
        dir='results/6-BXD_features/frip',
        bams=rules.correct_samples.params.dir
    shell:
        '''
        mkdir -p {params.dir}
        module load r-light/4.5.2
        
        echo "Estimate frip using <{input.gtf}>" >> {log}
        Rscript workflow/scripts/6.2_Peaks_frips.R  -b {params.bams}/coord/ -f {input.gtf} -o {params.dir} 

        echo "Logs saved in <{log}>" >> {log}
        '''

rule count_QC:
    '''
    R scripts used to agreggate and summarise the QC results from peak calling and counting reads
    '''
    input:
        counts=ancient(get_counts),
        gtf=ancient(rules.make_features.output.gtf),
        frip=ancient(get_frips),
        metadata = ancient(rules.check_variants.output.meta)
    output:
        plot='results/6-BXD_features/QC/QC_peaks_sizes.png',
    log:
        'logs/6-BXD_features/QC.log'
    benchmark:
        'benchmarks/6-BXD_features/QC.txt'
    resources:
        mem_mb = 2000,
        time = '01:00:00'
    threads: 1
    params:
        dir='results/6-BXD_features/QC',
        peaks='logs/5-BXD_peaks/',
    shell:
        '''
        mkdir -p {params.dir}
        module load r-light/4.5.2

        echo "Running QC for count" >> {log}
        Rscript workflow/scripts/6.3_Peaks_QC.R -m {input.metadata} -l {params.peaks} -c '{input.counts}' -f {input.gtf} -r '{input.frip}' -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}
        '''

rule multi_QC:
    '''
    MultiQC: Tool to automatically summarize analysis results for multiple tools and samples in a single report
    Reference: https://doi.org/10.1093/bioinformatics/btw354

    Performs:
    - Check log files and identifies tool used
    - Collect QC results from identified steps
    - Produces a single report summarising all QC
    '''
    input:
        check=rules.count_QC.output.plot
    output:
        report='results/6-BXD_features/MultiQC/multiqc_report.html',
    log:
        'logs/6-BXD_features/MultiQC.log'
    benchmark:
        'benchmarks/6-BXD_features/MultiQC.txt'
    resources:
        mem_mb = 2000,
        time = '01:00:00'
    threads: 1
    params:
        dir='results/6-BXD_features/MultiQC',
    shell:
        '''
        mkdir -p {params.dir}
        echo "Installing MultiQC python library" > {log}
        pip install multiqc 2>> {log}

        echo "Running MultiQC on all logs" > {log}
        multiqc --interactive ./ --ignore Rlibs --ignore libs --ignore workflow --ignore rawdata &>> {log}

        echo "Moving results to <{params.dir}>" >> {log}
        mv -f multiqc* {params.dir}/
        echo "Logs saved in <{log}>" >> {log}
        '''

