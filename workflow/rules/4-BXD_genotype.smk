'''
This Snakefile contains rules to genotype the BXD lines used:

get_genotypes: get BXD genotype annotation 
call_variants: Call relevant variants from available ATAC seq data
check_variants: Compare called variant to annotate one to confirm samples genotype
correct_samples: Correcte the two mislabbeled samples
'''

localrules: get_genotypes, check_variants


rule get_genotypes:
    '''
    Download BXD mm10 genotypes from qtl2data and cooreect mm10 issues

    Performs: 
    - Chromossome 10 correction using local mm9 annotation
    '''
    input:
        mm9_genotypes=ancient('rawdata/atac/BXD_mm9.geno.txt')
    output:
        genos='results/4-BXD_genotype/genotypes/genotypes.Rdata',  
        mm9='results/4-BXD_genotype/genotypes/BXD_mm9.geno.txt'       
    log:
        'logs/4-BXD_genotype/genotype.log'
    benchmark:
        'benchmarks/4-BXD_genotype/genotype.txt'
    resources:
        mem_mb = 100,
        time = '00:05:00'
    threads: 1
    params:
        dir='results/4-BXD_genotype/genotypes',
        metadata = config['atac_metadata']
    shell:
        '''
        mkdir -p {params.dir}
        module load r-light/4.5.2

        echo "Downloading genotypes for BXD lines" >> {log}
        wget https://raw.githubusercontent.com/rqtl/qtl2data/main/BXD/bxd.zip -P {params.dir}
        unzip {params.dir}/bxd.zip -d {params.dir}
        rm {params.dir}/bxd.zip

        cp {input.mm9_genotypes} {output.mm9}

        echo "Validating genotypes" >> {log}
        Rscript workflow/scripts/4.1-Genotype_get.R -m {params.metadata} -d {params.dir} 
        
        echo "Logs saved in <{log}>" >> {log}
        '''

rule call_variants:
    '''
    bfctools: Call genomic variants from ATAC reads
    Reference: https://samtools.github.io/bcftools/howtos/variant-calling.html

    Performs: 
    - Filter bam file for region 10Kb around the variant being tested
    - Calls BXD variants in the filtered bam file
    '''
    input:
        bam=ancient(rules.sort_aligment.output.full_co),
        geno=ancient(rules.get_genotypes.output.genos)
    output:
        idx  = protected('results/4-BXD_genotype/variants/{sample}_full.csv'),
    log:
        'logs/4-BXD_genotype/vars_{sample}.log'
    benchmark:
        'benchmarks/4-BXD_genotype/vars_{sample}.txt'
    resources:
        mem_mb = 3000,
        time = '00:20:00'
    threads: 20
    params:
        dir='results/4-BXD_genotype/variants'
    shell:
        '''
        module load r-light/4.5.2 bcftools/1.22 samtools/1.21
        mkdir -p {params.dir}
        
        echo "Calling genotypes" >> {log}
        Rscript workflow/scripts/4.2-Genotype_call.R -b {input.bam} -g {input.geno} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}

        '''

rule check_variants:
    '''
    Compare called variants with BXD genotypes
    bfctools: Call genomic variants from ATAC reads
    Reference: https://samtools.github.io/bcftools/howtos/variant-calling.html

    Performs: 
    - Filter bam file for region 10Kb around the variant being tested
    - Calls BXD variants in the filtered bam file
    '''
    input:
        ancient(get_variants)
    output:
        meta  = 'results/4-BXD_genotype/check/metadata_corrected.csv',
    log:
        'logs/4-BXD_genotype/check.log'
    benchmark:
        'benchmarks/4-BXD_genotype/check.txt'
    resources:
        mem_mb = 1000,
        time = '00:30:00'
    threads: 20
    params:
        dir='results/4-BXD_genotype/check',
        metadata = config['atac_metadata'],
        vars=rules.call_variants.params.dir,
        geno=rules.get_genotypes.output.genos,  
        bxd=rules.get_genotypes.params.dir
    shell:
        '''
        module load r-light/4.5.2
        mkdir -p {params.dir}
        
        echo "Commparing genotypes" >> {log}
        Rscript workflow/scripts/4.3-Genotype_check.R -m {params.metadata} -d {params.vars} -g {params.geno} -b {params.bxd} -o {params.dir}

        echo "Logs saved in <{log}>" >> {log}

        '''

rule correct_samples:
    '''
    Correct BXD line in mislabbeled samples identified in check_variants
    See: cat ./logs/4-BXD_genotype/check_4.3.log | grep -A4 mm10
    '''
    input:
        sample = ancient(correct_sample),
        check=ancient(rules.check_variants.output.meta),
        qc_check=ancient(rules.alignQC_summary.output.plot)
    output:
        full_co = protected('results/4-BXD_genotype/corrected_bams/coord/{sample}_full.co.bam'),
        sub1_co = protected('results/4-BXD_genotype/corrected_bams/coord/{sample}_sub1.co.bam'),
        sub2_co = protected('results/4-BXD_genotype/corrected_bams/coord/{sample}_sub2.co.bam'),
        full_nm = protected('results/4-BXD_genotype/corrected_bams/name/{sample}_full.nm.bam'),
        sub1_nm = protected('results/4-BXD_genotype/corrected_bams/name/{sample}_sub1.nm.bam'),
        sub2_nm = protected('results/4-BXD_genotype/corrected_bams/name/{sample}_sub2.nm.bam')
    log:
        'logs/4-BXD_genotype/correct_{sample}.log'
    benchmark:
        'benchmarks/4-BXD_genotype/correct_{sample}.txt'
    resources:
        mem_mb = 5000,
        time = '00:30:00'
    threads: 1
    params:
        dir='results/4-BXD_genotype/corrected_bams',
        source=rules.sort_aligment.params.dir
    shell:
        '''
        mkdir -p {params.dir}/coord {params.dir}/name

        echo "Copying <{input.sample}> files from <{params.source}> to <{params.dir}> " > {log}
        echo "cp {input.sample} {output.full_co}" >> {log}
        cp {input.sample} {output.full_co}

        echo "cp {input.sample}.bai {output.full_co}.bai" >> {log}
        cp {input.sample}.bai {output.full_co}.bai

        echo "cp $(echo {input.sample} | sed 's/full/sub1/g') {output.sub1_co}" >> {log}
        cp $(echo {input.sample} | sed 's/full/sub1/g') {output.sub1_co}
        
        echo "cp $(echo {input.sample} | sed 's/full/sub1/g').bai {output.sub1_co}.bai" >> {log}
        cp $(echo {input.sample} | sed 's/full/sub1/g').bai {output.sub1_co}.bai

        echo "cp $(echo {input.sample} | sed 's/full/sub2/g') {output.sub2_co}" >> {log}
        cp $(echo {input.sample} | sed 's/full/sub2/g') {output.sub2_co}
        
        echo "cp $(echo {input.sample} | sed 's/full/sub2/g').bai {output.sub2_co}.bai" >> {log}
        cp $(echo {input.sample} | sed 's/full/sub2/g').bai {output.sub2_co}.bai
        
        echo "cp $(echo {input.sample} | sed 's/coord/name/g') {output.full_nm}" >> {log}
        cp $(echo {input.sample} | sed 's/coord/name/g' | sed 's/\\.co\\./\\.nm\\./g') {output.full_nm}
        
        echo "cp $(echo {input.sample} | sed 's/full/sub1/g' | sed 's/coord/name/g') {output.sub1_nm}" >> {log}
        cp $(echo {input.sample} | sed 's/full/sub1/g' | sed 's/coord/name/g'| sed 's/\\.co\\./\\.nm\\./g') {output.sub1_nm}
        
        echo "cp $(echo {input.sample} | sed 's/full/sub2/g'| sed 's/coord/name/g') {output.sub2_nm}" >> {log}
        cp $(echo {input.sample} | sed 's/full/sub2/g'| sed 's/coord/name/g'| sed 's/\\.co\\./\\.nm\\./g') {output.sub2_nm}

        echo "Logs saved in <{log}>" >> {log}
    
        '''

 
