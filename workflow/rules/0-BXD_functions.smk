'''
This Snakefile contains the fucntion used to agreagate intputs and connect all rules
'''

# 1-BXD_fastq
#Input function used in rule readsQC_summary
def get_QCjson(wildcards):
    '''
    This function lists outputs from readsQC to connect the DAG
    '''
    return [f'results/1-BXD_fastq/QC/json/{sample}.json'
            for sample in config['atac_samples']]

# 3-BXD_bam
# Input function used in rule alignQC_summary
def get_QCbams(wildcards):
    '''
    This function lists outputs from bamQC and make_wig to connect the DAG
    '''
    return [f'results/3-BXD_bam/bamQC/{sample}_full.RData'
            for sample in config['atac_samples']] + \
           [f'results/3-BXD_bam/wigs/{sample}.bw'
            for sample in config['atac_samples']]     

# 4-BXD_genotype               
# Input function used in rule check_variants
def get_variants(wildcards):
    '''
    This function lists outputs from call_variants to connect the DAG
    '''
    return [f'results/4-BXD_genotype/variants/{sample}_full.csv'
            for sample in config['atac_samples']]

# Input function used in rule correct_samples
def correct_sample(wildcards):
    # Change the string value of the wildcard for the input path
    sample=wildcards.sample 
    if '038' in sample: 
        sample=sample.replace('038','039')
    elif '039' in sample: 
        sample=sample.replace('039','038')

    return [f'results/3-BXD_bam/sort/coord/{sample}_full.co.bam']

# 5-BXD_peaks

# Input function used in rule call_peaks/call_mitochondria
def get_groups(wildcards):
    '''
    This function lists outputs from correct_samples that match the call_peaks replicate group
    '''
    corrected_samples = [sample.replace('tmp','039') for sample in
     [sample.replace('039','038') for sample in
     [sample.replace('038','tmp') for sample in config['atac_samples']]
     ]]

    return [f'results/4-BXD_genotype/corrected_bams/name/{sample}_full.nm.bam'
            for sample in corrected_samples if wildcards.group in sample]       

# Input function used in rule call_parents
def subset_parents(wildcards):
    '''
    This function lists outputs from correct_samples that match a parent line and then randomly pick a triplicate from them
    '''
    import random
    # enforce parent lines
    samples = [sample for sample in config['atac_samples'] if 'DBA' in sample] + [sample for sample in config['atac_samples'] if 'C57Bl6' in sample]
    files = [f'results/4-BXD_genotype/corrected_bams/name/{sample}_full.nm.bam' for sample in samples if wildcards.group in sample]   
    
    seed=random.Random(wildcards.group+wildcards.i) # sees is tied to wildcard for reproducibility/prevent input from changing

    subset = seed.sample(files,k=3)

    return subset    

# Input function used in rule call_subsampled
def subset_subsample(wildcards):
    '''
    This function lists outputs from correct_samples that match the call_subsampled replicate group
    '''
    corrected_samples = [sample.replace('tmp','039') for sample in
     [sample.replace('039','038') for sample in
     [sample.replace('038','tmp') for sample in config['atac_samples']]
     ]]

    return [f'results/4-BXD_genotype/corrected_bams/name/{sample}_sub{wildcards.i}.nm.bam'
            for sample in corrected_samples if wildcards.group in sample]       

# 6-BXD_features
# Input function used in rule make_features
def get_peaks(wildcards):
    '''
    This function lists outputs from call_peaks to aggregate all narrowpeaks files 
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    return [f'results/5-BXD_peaks/peaks/{group}.narrowPeak.gz'
            for group in sample_groups]       

# Input function used in rule make_mito
def get_mito(wildcards):
    '''
    This function lists outputs from call_peaks to aggregate all narrowpeaks files 
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    return [f'results/5-BXD_peaks/mito/{group}.narrowPeak.gz'
            for group in sample_groups]       

# Input function used in rule make_parents
def get_parents(wildcards):
    '''
    This function lists outputs from call_parents to aggregate all narrowpeaks files 
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    parent_samples = [sample for sample in sample_groups if 'DBA' in sample] + [sample for sample in sample_groups if 'C57Bl6' in sample]
    return [f'results/5-BXD_peaks/parents/{group}_rand{wildcards.i}.narrowPeak.gz'
            for group in parent_samples]       

# Input function used in rule make_subsample
def get_subsample(wildcards):
    '''
    This function lists outputs from call_subsampled to aggregate all narrowpeaks files 
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    return [f'results/5-BXD_peaks/subset/{group}_sub{wildcards.i}.narrowPeak.gz'
            for group in sample_groups]       

# Input function used in rule count_QC
def get_counts(wildcards):
    '''
    This function lists outputs from count_features 
    '''
    corrected_samples = [sample.replace('tmp','039') for sample in
        [sample.replace('039','038') for sample in
        [sample.replace('038','tmp') for sample in config['atac_samples']]
        ]]

    return [f'results/6-BXD_features/counts/{sample}.count'
            for sample in corrected_samples]    

def get_frips(wildcards):
    '''
    This function lists outputs from frip to pull the connections and force execution
    '''
    peaks = ['results/6-BXD_features/frip/frip_list.csv']
    mito = ['results/6-BXD_features/frip/frip_mito.csv']
    parents = [f'results/6-BXD_features/frip/frip_parent{i}.csv' 
        for i in range(config['atac_parent_triplicate_permutation'])]
    subs =  [f'results/6-BXD_features/frip/frip_subsample{i}.csv'
        for i in [1,2]]

    return peaks+mito+parents+subs

# 9-BXD_QTL
# Input function used in rule qtl_peaks
def get_permutations(wildcards):
    '''
    This function lists outputs from qtl_perms  
    '''

    perms=int(config['qtl_permutation_split'])+1
    
    return [f'results/9-BXD_qtl/QTL/{wildcards.exp}{wildcards.chr}{wildcards.fc}_QTL_perm{i}.RData' for i in range(1, perms)]

# Input function used in rule qtl_aggregate
def get_qtls(wildcards):
    '''
    This function lists outputs from qtl_peaks  
    '''

    atac = [[f'results/9-BXD_qtl/sigs/atacchr{chr}{fc}_QTL_sigs.csv' for chr in list(range(1,20))+['X','Y']] for fc in ['', 'FC']]
    rna = [f'results/9-BXD_qtl/sigs/rna{fc}_QTL_sigs.csv' for fc in ['', 'FC']]
    sleep = ['results/9-BXD_qtl/sigs/sleep_QTL_sigs.csv']
    
    return sum(atac, rna) + sleep # sum flatten sthe 2D atac list

# 12-BXD_footprint
# Input function used in rule merge_bams
def get_bams(wildcards):
    '''
    This function merges all replicate bam files into a single sample
    '''
    samples = [sample for sample in config['atac_samples'] if wildcards.group in sample]
    
    return [f'results/3-BXD_bam/sort/coord/{sample}_full.co.bam'
            for sample in samples] 

# Input function used in rule scan_motifs
def get_footprints(wildcards):
    '''
    This function gathers all detected footprint files
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    bed=[f'results/12-BXD_footprints/footprints/{group}.bed' for group in sample_groups] 
    sub1=[f'results/12-BXD_footprints/footprints/{group}_sub1.bed' for group in sample_groups] 
    sub2=[f'results/12-BXD_footprints/footprints/{group}_sub2.bed' for group in sample_groups] 
    return bed+sub1+sub2

# Input function used in rule differential_heatmap
def all_bed(wildcards):
    '''
    This function gathers all bed required 
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    lines = list(set([line.split('_')[0] for line in sample_groups]))
        
    return [f'results/12-BXD_footprints/motifs/{line}_CTRL_mpbs.bed' for line in lines] + [f'results/12-BXD_footprints/motifs/{line}_SD_mpbs.bed' for line in lines]

def all_bam(wildcards):
    '''
    This function gathers all bed required 
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    lines = list(set([line.split('_')[0] for line in sample_groups]))

    return [f'results/12-BXD_footprints/bams/{line}_CTRL.bam' for line in lines] + [f'results/12-BXD_footprints/bams/{line}_SD.bam' for line in lines]

# Input function used in rule footprint_plots
def diffprints(wildcards):
    '''
    This function gathers all differential footprint analyses
    '''
    sample_groups = list(set(['_'.join(sample.split('_')[1:]) for sample in config['atac_samples']]))
    lines = list(set([line.split('_')[0] for line in sample_groups]))

    return [f'results/12-BXD_footprints/diff/{line}{wildcards.sub}/differential_statistics.txt' for line in lines] 

# 14-NFR_fastq
#Input function used in rule readsQC_summary_nrf
def get_QCjson_nrf(wildcards):
    '''
    This function lists outputs from readsQC to connect the DAG
    '''
    return [f'results/14-NRF_fastq/QC/json/{sample}.json'
            for sample in config['nrf_samples']]

# 16-NRF_bam
# Input function used in rule alignQC_summary_nrf
def get_QCbams_nrf(wildcards):
    '''
    This function lists outputs from bamQC and make_wig to connect the DAG
    '''
    return [f'results/16-NRF_bam/bamQC/{sample}.RData'
            for sample in config['nrf_samples']]