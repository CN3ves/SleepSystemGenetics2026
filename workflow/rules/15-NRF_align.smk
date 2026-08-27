'''
This Snakefile contains rules to align processed reads to reference genome:

star_index: Build index for aligment
align_rna: Align reads to genome
'''
rule star_index:
    '''
    Build index from genome and gene models
    '''

    output:
        idx  = 'results/15-NRF_align/index/Genome',
        gtf  = 'results/15-NRF_align/index/gencode.vM10.primary_assembly.annotation.gtf',
    log:
        'logs/15-NRF_align/index.log'
    benchmark:
        'benchmarks/15-NRF_align/index.txt'
    resources:
        mem_mb = 35000,
        time = '00:15:00'
    threads: 20
    params:
        dir='results/15-NRF_align/index',
        genome = config['atac_genome']

    shell:
        '''
        module load star/2.7.11b
        
        mkdir -p {params.dir}

        echo "Downloading index for {params.genome}" > {log}
        wget -P {params.dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/GRCm38.primary_assembly.genome.fa.gz
        wget -P {params.dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.primary_assembly.annotation.gff3.gz
        wget -P {params.dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.primary_assembly.annotation.gtf.gz
        
        gzip -d {params.dir}/*.gz

        echo "Build index" >> {log}
        STAR --runThreadN 20 --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir {params.dir} --genomeFastaFiles {params.dir}/GRCm38.primary_assembly.genome.fa --sjdbOverhang 150 --sjdbGTFfile {params.dir}/gencode.vM10.primary_assembly.annotation.gff3 >> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

rule align_rna:
    '''
    STAR: Spliced Transcripts Alignment to a Reference 
    Reference: 10.1093/bioinformatics/bts635; https://github.com/alexdobin/STAR

    Performs: 
    - Splice-aware read aligment to reference genome
    '''

    input:
        reads = ancient(rules.readsQC_nrf.output.trimmed),
        idx =  ancient(rules.star_index.output.idx),
        check = ancient(rules.readsQC_summary_nrf.output.plot)
    output:
        full = protected('results/15-NRF_align/aligment/{sample}_Aligned.out.sam'),
    log:
        'logs/15-NRF_align/align_{sample}.log'
    benchmark:
        'benchmarks/15-NRF_align/align_{sample}.txt'
    resources:
        mem_mb = 30000,
        time = '00:10:00'
    threads: 20
    params:
        dir='results/15-NRF_align/aligment',
        index=rules.star_index.params.dir,
        genome = rules.star_index.params.genome
    shell:
        '''
        module load star/2.7.11b

        mkdir -p {params.dir}

        echo "Aligning sample <{input.reads}> to <{params.genome}> transcripts annotation" >> {log}
        STAR --runMode alignReads --runThreadN 20 --genomeDir {params.index} --readFilesIn {input.reads} --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --winAnchorMultimapNmax 50 --outFileNamePrefix {params.dir}/{wildcards.sample}_ --outSAMunmapped Within  --quantTranscriptomeSAMoutput BanSingleEnd --quantMode GeneCounts &>> {log}
        
        echo "Logs saved in <{log}>" >> {log}
        '''

