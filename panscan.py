# %%
import os 
import pandas as pd
from pathlib import Path
import json
import subprocess
from multiprocessing import Pool
import random
import re
import pickle
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import numpy as np


import argparse
import subprocess

def call_novel_seq(vcf1, vcf2):
    os.system(f"perl findNovelSeq.pl {vcf1} {vcf2}")
        
        
def call_variant_analysis(vcf):
    os.system(f"perl preprocessVCF.pl {vcf}")

def process_and_plot_data(result_file):
    # Load data
    result_df = pd.read_csv(result_file, index_col=0)
    hprc = pd.read_csv('hprc-matrix.csv', skipfooter=1, engine='python', index_col=0)
    cpc = pd.read_csv('cpc-matrix.csv', index_col=0)
    cpc = cpc.drop(columns=['Unnamed: 2', 'Unnamed: 3', 'Unnamed: 4'])
    
    # Plotting duplications
    duplicates_count = (result_df != 0).sum()
    duplicates_count = duplicates_count.sort_values()
    plt.figure(figsize=(40, 24))
    duplicates_count.plot(kind='bar', color='#1F1F6A', edgecolor='black')
    plt.xlabel('Genome by number of duplications')
    plt.ylabel('Duplicated genes per genome')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xticks(rotation=0, color='white', ha='right')
    plt.tight_layout()
    plt.savefig('genome-by-dup.png')
    plt.show()

    # Venn diagram
    #common_indices = hprc.index.intersection(result_df.index).intersection(cpc.index)
    #venn_data = {
    #    '100': len(hprc.index.difference(result_df.index).difference(cpc.index)),
    #    '010': len(result_df.index.difference(hprc.index).difference(cpc.index)),
    #    '001': len(cpc.index.difference(hprc.index).difference(result_df.index)),
    #    '110': len(hprc.index.intersection(result_df.index).difference(cpc.index)),
    #    '101': len(hprc.index.difference(result_df.index).intersection(cpc.index)),
    #    '011': len(result_df.index.difference(hprc.index).intersection(cpc.index)),
    #    '111': len(common_indices)
    #}
    #plt.figure(figsize=(40, 40))
    #venn = venn3(subsets=venn_data, set_labels=('HPRC', 'APR', 'CPC'))
    #for text in venn.set_labels:
    #    text.set_fontsize(80)
    #for text in venn.subset_labels:
    #    text.set_fontsize(60)
    #plt.savefig('venn-comparison.png')
    #plt.show()
    
    # Venn diagram
    common_indices = hprc.index.intersection(result_df.index).intersection(cpc.index)
    venn_data = {
        '100': len(hprc.index.difference(result_df.index).difference(cpc.index)),
        '010': len(result_df.index.difference(hprc.index).difference(cpc.index)),
        '001': len(cpc.index.difference(hprc.index).difference(result_df.index)),
        '110': len(hprc.index.intersection(result_df.index).difference(cpc.index)),
        '101': len(hprc.index.difference(result_df.index).intersection(cpc.index)),
        '011': len(result_df.index.difference(hprc.index).intersection(cpc.index)),
        '111': len(common_indices)
    }
    plt.figure(figsize=(40, 40))
    venn = venn3(subsets=venn_data, set_labels=('HPRC', 'APR', 'CPC'))

    # Update font size for the Venn diagram labels
    for text in venn.set_labels:
        if text:  # Check if the label exists
            text.set_fontsize(80)
    if venn.subset_labels:  # Check if subset labels exist
        for text in venn.subset_labels:
            if text:  # Check if the subset label exists
                text.set_fontsize(60)
    plt.savefig('venn-comparison.png')
    plt.show()

    # Frequency comparison for each category
    def frequency_comparison(df1, df2, output_file, label1, label2):
        df1_freq = df1.apply(lambda row: (row != 0).sum() / len(row), axis=1).to_frame(name='Frequency_x')
        df2_freq = df2.apply(lambda row: (row != 0).sum() / len(row), axis=1).to_frame(name='Frequency_y')
        
        comparison_df = df1_freq.merge(df2_freq, left_index=True, right_index=True)
        comparison_df['Absolute_Difference'] = abs(comparison_df['Frequency_x'] - comparison_df['Frequency_y'])
        comparison_df['Percentage_Difference'] = (comparison_df['Absolute_Difference'] / comparison_df[['Frequency_x', 'Frequency_y']].max(axis=1)) * 100
        
        filtered_df = comparison_df[comparison_df['Percentage_Difference'] >= 5]
        sorted_df = filtered_df.sort_values(by='Percentage_Difference', ascending=False)
        top_5_diff = sorted_df.head(5)
        
        plt.figure(figsize=(40, 36))
        plt.barh(top_5_diff.index, top_5_diff['Frequency_x'], color='#1F1F6A', label=label1)
        plt.barh(top_5_diff.index, -top_5_diff['Frequency_y'], color='#EA4E15', label=label2)
        plt.xlabel('CNV frequency', labelpad=50)
        plt.legend()
        plt.rcParams['font.size'] = 70
        plt.tight_layout(pad=4)
        plt.savefig(output_file)
        plt.show()

    # Plot frequency comparison with CPC
    frequency_comparison(result_df, cpc, 'frequency-comparsion-with-cpc.png', 'Result', 'CPC')

    # Plot frequency comparison with HPRC
    frequency_comparison(result_df, hprc, 'frequency-comparsion-with-hprc.png', 'Result', 'HPRC')



def run_variant_analysis(vcf_file):
    """
    Runs the panscan variant_analysis command.
    """
    try:
        # variant_analysis_cmd = f"panscan variant_analysis {vcf_file}"
        # subprocess.run(variant_analysis_cmd, shell=True, check=True)
        call_variant_analysis(vcf_file)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running panscan variant_analysis: {e}")

def run_gene_dup(csv_file):
    """
    Runs the panscan gene_dup command.
    """
    try:
        process_and_plot_data(csv_file)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running panscan gene_dup: {e}")

def run_panscan_novel_seq(vcf1, vcf2):
    """
    Runs the panscan novel_seq command.
    """
    try:
        call_novel_seq(vcf1, vcf2)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running panscan novel_seq: {e}")

def main():
    parser = argparse.ArgumentParser(description="Panscan tool to analyze complex loci, perform variant analysis, and gene duplication detection.")
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Complex command
    complex_parser = subparsers.add_parser("complex", help="Analyze complex loci from a VCF file.")
    complex_parser.add_argument("vcf_file", help="Path to the VCF file generated from Minigraph-cactus pipeline.")
    complex_parser.add_argument("-a", type=int, default=5, help="Number of alleles to define a complex site (default: 5).")
    complex_parser.add_argument("-n", type=int, default=1, help="Number of 10kb variants to define a complex site (default: 1).")
    complex_parser.add_argument("-s", type=int, default=10000, help="Minimum size of a variant to define a complex site (default: 10000).")
    complex_parser.add_argument("--regions", action="store_true", help="List complex regions.")
    complex_parser.add_argument("-l", type=int, default=100000, help="Length of the region to define a complex region (default: 100000).")
    complex_parser.add_argument("--sites", type=int, default=1, help="Number of complex sites in a region to define it as complex (default: 1).")
    complex_parser.add_argument("--sv", type=int, default=1, help="Number of secondary SVs in a region to define it as complex (default: 1).")
    
    # Novel seq command
    novel_seq_parser = subparsers.add_parser("novel_seq", help="Detect novel sequences from two VCF files.")
    novel_seq_parser.add_argument("vcf1", help="Path to the first VCF file.")
    novel_seq_parser.add_argument("vcf2", help="Path to the second VCF file.")
    
    # Variant analysis command
    variant_analysis_parser = subparsers.add_parser("variant_analysis", help="Perform variant analysis on a VCF file.")
    variant_analysis_parser.add_argument("vcf_file", help="Path to the VCF file for variant analysis.")
    variant_analysis_parser.add_argument("--output", help="Output directory to save the results.")

    
    # Gene duplication command
    gene_dup_parser = subparsers.add_parser("gene_dup", help="Detect gene duplications from a CSV file.")
    gene_dup_parser.add_argument("csv_file", help="Path to the CSV file for gene duplication detection.")
    
    args = parser.parse_args()
    
    if args.command == "complex":
        run_panscan_complex(args.vcf_file, args.a, args.n, args.s, args.l, args.sites, args.sv)
    elif args.command == "variant_analysis":
        run_variant_analysis(args.vcf_file)
    elif args.command == "gene_dup":
        run_gene_dup(args.csv_file)
    elif args.command == "novel_seq":
        run_panscan_novel_seq(args.vcf1, args.vcf2)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

# %% [markdown]
# Find complex regions
# print('running panscan ')
# %%
# vcf_file = 'apr_review_v1_2902_chm13.vcf'

# %%
def extract_info(info_str, key):
    for field in info_str.split(';'):
        if field.startswith(key + '='):
            return field[len(key)+1:]
    return None

def allele_lengths(ref, alt):
    alts = alt.split(',')
    return [len(a) for a in alts]

# %% [markdown]
# Find sites that fit a certain criteria

# %%
def get_all_sv_sites(vcf_file):
    all_sites = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                ref = columns[3]
                alt = columns[4]
                lengths = allele_lengths(ref, alt)
                
                if any((l-len(ref)) > 50 for l in lengths):
                    all_sites.append(line.strip().split('\t'))

        return pd.DataFrame(all_sites)
    

# %%
def get_interesting_sites(vcf_file):
    interesting_sites = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):  # Exclude header lines
                columns = line.strip().split('\t')
                ref = columns[3]
                alt = columns[4]
                
                lengths = allele_lengths(ref, alt)

                if len(lengths) >= 5 and any((l-len(ref)) > 10000 for l in lengths):
                    interesting_sites.append(line.strip())

    return interesting_sites

# %%
# interesting_sites = get_interesting_sites(vcf_file)

# %%
def get_df(sites):
    data = []
    for site in sites:
        data.append((site.split()[0], int(site.split()[1])))
    
    return pd.DataFrame(data, columns=['chrom', 'pos'])

# site_df = get_df(interesting_sites)

# # %%
# all_sites_df = get_all_sv_sites(vcf_file)

# %%
# all_sites_df[1] = pd.to_numeric(all_sites_df[1])

# %%

# window_size = 50000
# regions = []
# for i,site in site_df.iterrows():

#     chrom = site['chrom']
#     pos = site['pos']

#     sites_in_window = all_sites_df.loc[(all_sites_df[0] == chrom ) & (all_sites_df[1]>pos-window_size) & (all_sites_df[1]<pos+window_size)]
    
#     if sites_in_window.shape[0] >= 2:
#         regions.append([chrom, pos-window_size, pos+window_size])
        

# %%
# window_size = 100000

def get_complex_regions(window_size): 
    regions = []
    for chrom in site_df['chrom'].unique():
        temp_df = site_df.loc[site_df['chrom'] == chrom]
        max_position = temp_df['pos'].max()
        min_position = temp_df['pos'].min()
        
        ptr = min_position
        end = max_position
        
        
        
        while ptr < end:
            sites_in_window = temp_df.loc[(temp_df['pos']>ptr) & (temp_df['pos']<ptr+window_size)]
            if sites_in_window.shape[0] >= 2:
                regions.append((chrom, ptr, ptr+window_size))
                
            ptr += window_size

    return regions

# %%
# regions = get_complex_regions(window_size)

# %% [markdown]
# Find all the genes in these regions

# %%
# gff3_file = 'chm13_genomic_chr.gff3'

def extract_genes_from_gff3(filename):
    genes = []
    with open(filename, 'r') as gff:
        for line in gff:
            if not line.startswith("#"):  # Ignore comment lines
                parts = line.strip().split("\t")
                if parts[2] == "gene" or parts[2] == "pseudogene":
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    attributes = parts[8].split(";")
                    gene_name = None
                    is_protein_coding = False
                    for attribute in attributes:
                        if attribute.startswith("ID="):
                            gene_name = attribute.split("=")[1]
                        if attribute == "gene_biotype=protein_coding":
                                is_protein_coding = True
                            
                    if gene_name:
                        genes.append((chrom,start, end, '-'.join(gene_name.split('-')[1:]), is_protein_coding))
    return genes


# %%
# all_genes_df = pd.DataFrame(extract_genes_from_gff3(gff3_file))

# %%
def get_region_genes(regions, gene_df):
    region_genes = []
    for chr, start,end in regions:
        chr_genes = all_genes_df.loc[(all_genes_df[0] == chr) & (all_genes_df[1] >= start) & (all_genes_df[2] <= end)]
        pc_genes = list(chr_genes.loc[chr_genes[4] == True, 3])
        region_genes.append((chr, start, end, ','.join(chr_genes[3]), ','.join(list(pc_genes))))
     
        # for i, row in chr_genes.iterrows():
        #     gene = row[3]
        #     pc = row[4]
        #     region_genes.append((chr, start, end, gene, pc))

    return pd.DataFrame(region_genes, columns=['chr','start','end', 'all_genes', 'pc_genes'])


# %%
# region_genes = get_region_genes(regions, all_genes_df)

# # %%
# region_genes['pc_gene_list'] = region_genes['pc_genes'].apply(lambda x: x.split(','))

# # %%
# filtered_region_genes = region_genes[region_genes['pc_gene_list'].apply(len) >= 2]
# filtered_region_genes = filtered_region_genes['chr	start	end	all_genes	pc_genes'.split('\t')]

# # %% [markdown]
# # find regions that have more than one haplotypes

# # %%
# regions = [tuple(row) for row in filtered_region_genes[['chr', 'start', 'end']].itertuples(index=False, name=None)]

# %%
def produce_plottable(gfab, gff3_file, region, cutpoints, graph_base, viz_output, connected_output, ref, chrom, start, end, query_region, workdir, gene_alignments, df_all):
    Path(workdir).mkdir(parents=True, exist_ok=True)
    connected_command = f'gfabase sub {gfab} -o {connected_output} {query_region} --range --cutpoints {cutpoints} --view --connected'
    viz_command = f'gfabase sub {gfab} -o {viz_output} {query_region} --range --view --cutpoints {cutpoints}'

    print(connected_command)
    def run_command(command_string):
        cmd = command_string.split(' ')
        try:
            result = subprocess.run(cmd, check=True, text=True)
            return (cmd, None)
        except subprocess.CalledProcessError as e:
            return (cmd, None, e)
    if not os.path.isfile(viz_output):
        os.system(viz_command)
    
    if not os.path.isfile(connected_output):
        os.system(connected_command)
  
        
        
    def extract_genes_from_region(filename, chrom_filter, start_filter, end_filter):
        genes = []
        with open(filename, 'r') as gff:
            for line in gff:
                if not line.startswith("#"):  # Ignore comment lines
                    parts = line.strip().split("\t")
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    
                    # Filter for the specific region
                    if chrom == chrom_filter and parts[2] == "gene" and start_filter <= start and end_filter >= end:
                        attributes = parts[8].split(";")
                        gene_name = None
                        for attribute in attributes:
                            if attribute.startswith("gene_name="):
                                gene_name = attribute.split("=")[1]
                                break
                        if gene_name:
                            genes.append((chrom,start, end, gene_name))
        return genes

    def extract_protein_coding_genes_from_region(filename, chrom_filter, start_filter, end_filter):
        genes = []
        with open(filename, 'r') as gff:
            for line in gff:
                if not line.startswith("#"):  # Ignore comment lines
                    parts = line.strip().split("\t")
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    
                    # Filter for the specific region
                    if chrom == chrom_filter and parts[2] == "gene" and start_filter <= start and end_filter >= end:
                        attributes = parts[8].split(";")
                        gene_name = None
                        is_protein_coding = False
                        for attribute in attributes:
                            if attribute.startswith("gene_name="):
                                gene_name = attribute.split("=")[1]
                            if attribute == "gene_biotype=protein_coding":
                                is_protein_coding = True

                        if gene_name and is_protein_coding:
                            
                            genes.append((chrom,start, end, gene_name))

        return genes



    print('got_genes')

    all_genes_in_region = extract_genes_from_region(gff3_file, chrom, start, end) 
    pc_genes_in_region= extract_protein_coding_genes_from_region(gff3_file, chrom, start, end)
    all_genes = [items[-1] for items in all_genes_in_region]
    pc_genes = [items[-1] for items in pc_genes_in_region]
    
    print('getting genes')
    def hex_to_rgb(value):
        value = value.lstrip('#')
        length = len(value)
        return tuple(int(value[i:i + length // 3], 16) for i in range(0, length, length // 3))

    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])

    def gradient(start_hex, end_hex, num):
        start_rgb = hex_to_rgb(start_hex)
        end_rgb = hex_to_rgb(end_hex)
        colors = []
        for i in range(num):
            r = int(start_rgb[0] + (i / (num - 1)) * (end_rgb[0] - start_rgb[0]))
            g = int(start_rgb[1] + (i / (num - 1)) * (end_rgb[1] - start_rgb[1]))
            b = int(start_rgb[2] + (i / (num - 1)) * (end_rgb[2] - start_rgb[2]))
            colors.append(rgb_to_hex((r,g,b)))
        return colors
    
    
    def parse_string(s):
        return set(map(int, ''.join(c if c.isdigit() else ' ' for c in s).split()))

    
    
    def get_segement_colors_from_alignment(df_all, section_genes, color_files_dir):

        df_all = df_all[(df_all[0].isin(section_genes)) & (df_all[11] > 0)]

        df_all['index'] = df_all.groupby(0).cumcount() + 1
        
        df_all['Numbers'] = df_all[5].apply(parse_string)

        to_drop = []
        for i in range(len(df_all)):
            for j in range(i+1, len(df_all)):
                if df_all.iloc[i]['Numbers'].intersection(df_all.iloc[j]['Numbers']):
                    if df_all.iloc[i][1] < df_all.iloc[j][1]:
                        to_drop.append(df_all.iloc[i].name)
                    else:
                        to_drop.append(df_all.iloc[j].name)
                        
        df_all.drop(to_drop, inplace=True)
        df_all.drop(columns=['Numbers'], inplace=True)
        
        df = df_all[[0, 5]]
        df = df.groupby(0).agg(lambda x: ''.join(x)).reset_index()
        
        def random_color():
            return "#%06x" % random.randint(0, 0xFFFFFF)

        df['color'] = df.apply(lambda x: random_color(), axis=1)
        
        start_hex = '#073f40'
        end_hex = '#0e0d2b'
        
        


        items=[]
        for i, row in df.iterrows():
            gene = str(row[0])
            nodes = re.split(r'[<>]',row[5])
            colors = gradient(start_hex, end_hex, len(nodes))   
            with open(f'{color_files_dir}/gene_colors/{gene}.csv', 'w') as f:
                
            
                for node, color in zip(nodes, colors):
                    if node != '':
                        f.write(f'{node}, {row["color"]}\n')
                        items.append((node, row['color']))
        return pd.DataFrame(items), df_all
    
    
    color_files_dir = f'{workdir}/color_csvs'
    Path(color_files_dir).mkdir(parents=True, exist_ok=True)
    Path(f'{color_files_dir}/gene_colors').mkdir(parents=True, exist_ok=True)
    all_genes, all_genes_all = get_segement_colors_from_alignment(df_all,all_genes, color_files_dir)
    all_genes.to_csv(f'{color_files_dir}/all_colors.csv', header=None, index=None)
    pc_genes, pc_genes_all = get_segement_colors_from_alignment(df_all, pc_genes, color_files_dir)
    pc_genes.to_csv(f'{color_files_dir}/pc_colors.csv', header = None, index=None)
    viz_list = []
    sample_walks = {} 
    with open(viz_output, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line[0] == 'S':
                segment = line.split('\t')[1]
                viz_list.append(segment)
                
    with open(connected_output, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line[0] == 'W':
                items = line.split('\t')
                nodes = re.split(r'[><]', items[6].strip())
                if items[1] not in sample_walks:
                    sample_walks[items[1]] = []
                
                sample_walks[items[1]].extend([node for node in nodes if node != '' ])

    genes_walks = {}

    for i, row in pc_genes_all.iterrows():
        gene = row[0]
        gene_index = row['index']
        path = re.split(r'[<>]',row[5])
        s_i = f'{gene}_{gene_index}'
        if gene not in genes_walks:
            genes_walks[s_i] = []
            
        genes_walks[s_i].extend(path)
        
    for k,v in sample_walks.items():
        sample_walks[k] = set(viz_list).intersection(set(v))
    for k, v in sample_walks.items():
        colors = gradient('#073f40', '#2b184d', len(v))
        with open(f'{color_files_dir}/{k}.walks.csv', 'w') as f:
            for item, color in zip(v, colors):
                f.write(f'{item},{color}\n')
    sample_haps = {}
    for k, v in sample_walks.items():
        if k not in sample_haps:
            sample_haps[k] = []
            
        for kj, vj in genes_walks.items():
            if len(set(v).intersection(set(vj))) > 1:
                sample_haps[k].append(kj)
                
    hap_directory = f'{workdir}/sample_haps/'
    Path(hap_directory).mkdir(parents=True, exist_ok=True)
    with open(f'{hap_directory}/sample_haps.csv', 'w') as f:
        for i,haps in sample_haps.items():
            f.write(f'{i} : {",".join(haps)}\n')
            
    haps = open(f'{hap_directory}/sample_haps.csv').read().splitlines()
    haps_dict = {hap.split(':')[0] : tuple(hap.split(':')[-1].split(',')) for hap in haps}
            
    def count_unique_lists(d):
        count_dict = {}
    
    
    
        for key, value in d.items():
            # Convert list to tuple
            t = tuple(value)
            
            if t in count_dict:
                count_dict[t]['count'] += 1
                count_dict[t]['keys'].append(key)
            else:
                count_dict[t] = {'count': 1, 'keys': [key]}
        
        return count_dict
    
    result = count_unique_lists(haps_dict)
    json_ready_result = {str(list(key)): value for key, value in result.items()}

    # Save the result to a JSON file
    with open(f'{hap_directory}/haps_counts.json', 'w') as f:
        json.dump(json_ready_result, f, indent=4)

    return result


def call_novel_seq(vcf1, vcf2):
    subprocess.call(f"perl findNovelSeq.pl {vcf1} {vcf2}")
        
        
def call_variant_analysis(vcf):
    subprocess.run(f"perl preprocessVCF.pl {vcf}")

