# python3
# Henry North April 2025

# Required Libraries
import gzip
import pandas as pd
import os
import time
import argparse

# Read in arguments 
parser = argparse.ArgumentParser(description='Filter genotypes by individual depth thresholds.')

parser.add_argument('-OUTDIR', '--output_directory', required=True, type=str, help='Full path to output directory')
parser.add_argument('-c', '--chromosome', required=True, type=str, help='Chromosome number')
parser.add_argument('-mDP', '--minDP_constant', type=int, help='minimum depth of coevrage for all genotypes')
parser.add_argument('-DPr', '--DP_ratio_threshold', type=int, help='DP will be filtered to be mean_depth/DP_ratio_threshold \
    < DP < mean_depth*DP_ratio_threshold, e.g. if the mean autosomal depth of coverage for an individual is 100, and \
    the DP_ratio_threshold is 5, then the depth of coverage for that individual genotype will be constrained to 20 < DP <500.')

args = parser.parse_args()

OUTDIR = args.output_directory
chromosome = args.chromosome
minDP_constant = args.minDP_constant
DP_ratio_threshold = args.DP_ratio_threshold

# set start time
start_time = time.time()

vcf_path = f'{OUTDIR}/c{chromosome}_SNPs_baseFilter1.vcf.gz'

header_subsetter = f"zgrep '^#CHROM' {vcf_path}" # get the header line

# Use the header lines to define dataframe columns
columns = os.popen(header_subsetter).read().lstrip('#').strip().split('\t')

# Read in the VCF information, adding header column names
vcf=pd.read_csv(vcf_path, sep='\t', compression='gzip', header=None, names=columns, comment='#')

# Define the context columns
context_columns=['CHROM','POS','ID','REF', 'ALT','QUAL','FILTER','INFO','FORMAT']

# Define the IND columns as all those that aren't context columns
ind_columns = [col for col in vcf.columns if col not in context_columns]

# Read in mean individual depth data and map to 
ind_depth_path='/home/hln33/rds/jiggins-rds-fT31urweTx0/projects/project_helicoverpa/publication_analyses/00_filtering/biallelic_site_VCFs/autosomal_biallelic_snps_complete.idepth'
ind_depth=pd.read_csv(ind_depth_path, sep='\t')

# Confirm that all the elements of ind_columns are in ind_depth['INDV']
for col in ind_columns:
    if col not in ind_depth['INDV'].values:
        print(f"Error: {col} is not in the INDV column of the depth file.")
        exit(1) 

# Create a dictionary for ind_depth lookup
mean_depth_dict = dict(zip(ind_depth['INDV'], ind_depth['MEAN_DEPTH']))

# Define a function to filter a genotype field based on depth thresholds
def genotype_filter(GENOTYPE_ij, FORMAT_ij, mean_depth):
    
    # If genotype is missing, do not apply the filter
    if GENOTYPE_ij.startswith('./.'):
        return GENOTYPE_ij
    
    # Otherwise, extract the genotypic depth of coverage (DP) using the FORMAT field
    format_keys = FORMAT_ij.split(':')
    geno_values = GENOTYPE_ij.split(':')

    
    # If this is missing or misspecified, exit
    if 'DP' not in format_keys:
        print(f"Error: DP not found in FORMAT field '{FORMAT_ij}'")
        exit(1)
        
    try:
        dp_index = format_keys.index('DP')
        depth = int(geno_values[dp_index])
    except (IndexError, ValueError):
        print(f"Error: Unexpected format in the following:'{FORMAT_ij}'")
        exit(1)

    ######################## the filter ########################
    if depth < minDP_constant or depth < mean_depth / DP_ratio_threshold or depth > mean_depth * DP_ratio_threshold:
        geno_values[0] = './.'
    ############################################################
    
    return ':'.join(geno_values)


# Apply filter column-wise
for col in ind_columns:
    mean_depth = mean_depth_dict[col] # get mean depth for the individual
    
    vcf[col] = vcf.apply(lambda row: genotype_filter(row[col], row['FORMAT'], mean_depth), axis=1) #apply the filter function
    
    #print(f"Filter complete for individual {col}") # notify completion


# Now output the modified VCF

# Extract the original VCF header lines (all ## and #CHROM lines)
vcf_header_lines = []

with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            vcf_header_lines.append(line)
        elif line.startswith('#CHROM'):
            vcf_header_lines.append(line)
            break  # stop after the column header line

# Define output path
output = f'{OUTDIR}/c{chromosome}_SNPs_baseFilter1_4.vcf.gz'

# Write output VCF with updated genotype calls
with gzip.open(output, 'wt') as out_f:
    # Write header lines
    for line in vcf_header_lines:
        out_f.write(line)
    
    # Write the modified DataFrame
    vcf.to_csv(out_f, sep='\t', index=False, header=False)

# Record completion time
end_time = time.time()
elapsed = end_time - start_time
print(f"\nCompleted in {elapsed:.2f} seconds.")
