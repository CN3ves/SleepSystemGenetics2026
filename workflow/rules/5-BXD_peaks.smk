'''
This Snakefile contains rules to call peaks:

get_genrich: Downloads Genrich and the mm10 blacklisted regions
call_peaks: Calls peaks from replicate samples
call_mitochondria (optional): Used to test the impact of removing mitochondrial DNA on peaks called
call_parents (optional): Used to test the effect of having more parent line replicates
call_subsampled (optional): Used to test the effect of number of reads 
'''

localrules: get_genrich

rule get_genrich:
    '''
    Downloads mm10 blacklist and Genrich software
    reference: https://doi.org/10.1038/s41598-019-45839-z; https://github.com/Boyle-Lab/Blacklist
    '''
    output:
        blacklist  = protected('results/5-BXD_peaks/caller/mm10-blacklist.v2.bed'),
    log:
        'logs/5-BXD_peaks/caller.log'
    benchmark:
        'benchmarks/5-BXD_peaks/caller.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/5-BXD_peaks/caller',
        blacklist = 'mm10-blacklist.v2.bed.gz'

    shell:
        '''
        mkdir -p {params.dir}

        echo "Downloading <{params.blacklist}>" >> {log}
        wget  https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/{params.blacklist} -P {params.dir}  &>> {log}
        echo "Unpacking {params.blacklist}" >> {log}
        gunzip {params.dir}/{params.blacklist} &>> {log}

        echo "Installing <Genrich>" >> {log}
        git clone https://github.com/jsh58/Genrich.git {params.dir}/Genrich &>> {log}
        (cd {params.dir}/Genrich/ && bash ./Makefile)  &>> {log}

        echo "Logs saved in <{log}>" >> {log}

        '''

rule call_peaks:
    '''
    Genrich: Peak-caller for genomic enrichment assays (e.g. ChIP-seq, ATAC-seq).
    Reference: https://github.com/jsh58/Genrich

    Performs: 
    - Peak calling on replicates statistics
    '''
    input:
        bams=ancient(get_groups),
        blacklist=ancient(rules.get_genrich.output.blacklist)
    output:
        genos=protected('results/5-BXD_peaks/peaks/{group}.narrowPeak.gz'),
        pval=protected('results/5-BXD_peaks/peaks/logs/{group}.pval.log.gz'),
        pileup=protected('results/5-BXD_peaks/peaks/pileups/{group}.pileup.gz'),
        bed=protected('results/5-BXD_peaks/peaks/beds/{group}.bed.gz'),
        dups=protected('results/5-BXD_peaks/peaks/logs/{group}.dups.gz')
    log:
        'logs/5-BXD_peaks/peaks_{group}.log'
    benchmark:
        'benchmarks/5-BXD_peaks/peaks_{group}.txt'
    resources:
        mem_mb = 40000,
        time = '01:45:00'
    threads: 1
    params:
        dir='results/5-BXD_peaks/peaks',
        genrich=rules.get_genrich.params.dir
    shell:
        '''
        mkdir -p {params.dir}/pileups {params.dir}/beds {params.dir}/logs

        echo "Calling peaks for {input.bams}" >> {log}
        {params.genrich}/Genrich/Genrich -t '{input.bams}'  -o {output.genos} -f {output.pval} -k {output.pileup} -b  {output.bed} -m 30 -r -R {output.dups} -q 0.05 -y -j -e chrM -E {input.blacklist} -z -v &>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

rule call_mitochondria:
    '''
    Same as call_peaks, but does not exclude mitochondrial DNA (-e chrM).
    '''
    input:
        bams=ancient(get_groups),
        blacklist=ancient(rules.get_genrich.output.blacklist)
    output:
        genos=protected('results/5-BXD_peaks/mito/{group}.narrowPeak.gz'),
        pval=protected('results/5-BXD_peaks/mito/logs/{group}.pval.log.gz'),
        pileup=protected('results/5-BXD_peaks/mito/pileups/{group}.pileup.gz'),
        bed=protected('results/5-BXD_peaks/mito/beds/{group}.bed.gz'),
        dups=protected('results/5-BXD_peaks/mito/logs/{group}.dups.gz')
    log:
        'logs/5-BXD_peaks/mito_{group}.log'
    benchmark:
        'benchmarks/5-BXD_peaks/mito_{group}.txt'
    resources:
        mem_mb = 40000,
        time = '01:45:00'
    threads: 1
    params:
        dir='results/5-BXD_peaks/mito',
        genrich=rules.get_genrich.params.dir
    shell:
        '''
        mkdir -p {params.dir}/pileups {params.dir}/beds {params.dir}/logs

        echo "Calling peaks with mitochondrial DNA for {input.bams}" >> {log}
        {params.genrich}/Genrich/Genrich -t '{input.bams}'  -o {output.genos} -f {output.pval} -k {output.pileup} -b  {output.bed} -m 30 -r -R {output.dups} -q 0.05 -y -j -E {input.blacklist} -z -v &>> {log}
   
        echo "Logs saved in <{log}>" >> {log}
        '''

rule call_parents:
    '''
    Same as call_peaks, but only considers random triplicates of the parent lines (same N as BXD lines).
    NOTE: i defines a new file, thus forcing a new run with a new random subset
    '''
    input:
        bams=ancient(subset_parents),
        blacklist=ancient(rules.get_genrich.output.blacklist)
    output:
        genos=protected('results/5-BXD_peaks/parents/{group}_rand{i}.narrowPeak.gz'),
        pval=protected('results/5-BXD_peaks/parents/logs/{group}_rand{i}.pval.log.gz'),
        pileup=protected('results/5-BXD_peaks/parents/pileups/{group}_rand{i}.pileup.gz'),
        bed=protected('results/5-BXD_peaks/parents/beds/{group}_rand{i}.bed.gz'),
        dups=protected('results/5-BXD_peaks/parents/logs/{group}_rand{i}.dups.gz')
    log:
        'logs/5-BXD_peaks/parents_{group}_{i}.log'
    benchmark:
        'benchmarks/5-BXD_peaks/parents_{group}_{i}.txt'
    resources:
        mem_mb = 40000,
        time = '01:45:00'
    threads: 1
    params:
        dir='results/5-BXD_peaks/parents',
        genrich=rules.get_genrich.params.dir
    shell:
        '''
        mkdir -p {params.dir}/pileups {params.dir}/beds {params.dir}/logs

        echo "Calling peaks for {input.bams}" >> {log}
        {params.genrich}/Genrich/Genrich -t '{input.bams}'  -o {output.genos} -f {output.pval} -k {output.pileup} -b  {output.bed} -m 30 -r -R {output.dups} -q 0.05 -y -j -e chrM -E {input.blacklist} -z -v &>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

rule call_subsampled:
    '''
    Same as call_peaks, but only considers random triplicates of the parent lines (same N as BXD lines).
    NOTE: i defines a new file, thus forcing a new run with a new random subset
    '''
    input:
        bams=ancient(subset_subsample),
        blacklist=ancient(rules.get_genrich.output.blacklist)
    output:
        genos=protected('results/5-BXD_peaks/subset/{group}_sub{i}.narrowPeak.gz'),
        pval=protected('results/5-BXD_peaks/subset/logs/{group}_sub{i}.pval.log.gz'),
        pileup=protected('results/5-BXD_peaks/subset/pileups/{group}_sub{i}.pileup.gz'),
        bed=protected('results/5-BXD_peaks/subset/beds/{group}_sub{i}.bed.gz'),
        dups=protected('results/5-BXD_peaks/subset/logs/{group}_sub{i}.dups.gz')
    log:
        'logs/5-BXD_peaks/subset_{group}_{i}.log'
    benchmark:
        'benchmarks/5-BXD_peaks/subset_{group}_{i}.txt'
    resources:
        mem_mb = 40000,
        time = '01:45:00'
    threads: 1
    params:
        dir='results/5-BXD_peaks/subset',
        genrich=rules.get_genrich.params.dir
    shell:
        '''
        mkdir -p {params.dir}/pileups {params.dir}/beds {params.dir}/logs

        echo "Calling peaks for {input.bams}" >> {log}
        {params.genrich}/Genrich/Genrich -t '{input.bams}'  -o {output.genos} -f {output.pval} -k {output.pileup} -b  {output.bed} -m 30 -r -R {output.dups} -q 0.05 -y -j -e chrM -E {input.blacklist} -z -v &>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''