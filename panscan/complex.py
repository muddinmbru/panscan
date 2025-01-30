import argparse
import pandas as pd  # type: ignore
from pathlib import Path
import subprocess
import os
import random
import re
import json

def extract_info(info_str, key):
    for field in info_str.split(';'):
        if field.startswith(key + '='):
            return field[len(key)+1:]
    return None

def allele_lengths(ref, alt):
    alts = alt.split(',')
    return [len(a) for a in alts]

def get_interesting_alleles(vcf_file, n_10kb, site_size):
    if n_10kb == None:
        n_10kb = 5
    if site_size == None:
        site_size=10000

    interesting_alleles = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):  # Exclude header lines
                columns = line.strip().split('\t')
                ref = columns[3]
                alt = columns[4]
                
                lengths = allele_lengths(ref, alt)

                if len(lengths) >= n_10kb and any((l-len(ref)) > site_size for l in lengths):
                    interesting_alleles.append(line.strip())

    return interesting_alleles

def get_df(sites):
    data = []
    for site in sites:
        data.append((site.split()[0], int(site.split()[1])))
    
    return pd.DataFrame(data, columns=['chrom', 'pos'])



def find_complex_sites(vcf_file, n_10kb, site_size):
    interesting_alleles = get_interesting_alleles(vcf_file, n_10kb, site_size)
    locations = [(allele.split()[0], allele.split()[1]) for allele in interesting_alleles]
    regions = [(region[0], int(region[1] )- 10000, int(region[1]) + 10000) for region in locations]

    site_df = get_df(interesting_alleles)

    return site_df

def find_interesting_regions(site_df, window_size):
    if window_size == None:
        window_size = 100000

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




def produce_plottable(*task):
    gfab, gff3_file, region, cutpoints, graph_base, viz_output, connected_output, ref, chrom, start, end, query_region, workdir, gene_alignments, df_all, sep_pattern = task
    
    Path(workdir).mkdir(parents=True, exist_ok=True)
    
    # Modify query region to match gfabase requirements
    modified_query_region = f'{ref}{sep_pattern}{chrom}:{start}-{end}'
    
    connected_command = f'gfabase sub {gfab} -o {connected_output} {modified_query_region} --range --cutpoints {cutpoints} --view --connected'
    viz_command = f'gfabase sub {gfab} -o {viz_output} {modified_query_region} --range --view --cutpoints {cutpoints}'
    
    print(connected_command)
    print(viz_command)
    
    os.system(connected_command)
    os.system(viz_command)
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
        #####debug zome
        print(df_all)
        print("Section genes:", section_genes)
        print("Initial df_all shape:", df_all.shape)
        #####
        df_all = df_all[(df_all[0].isin(section_genes)) & (df_all[11] > 0)]
        print("Filtered df_all shape:", df_all.shape)    
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

def run_panscan_complex(vcf_file, a, n, s, l, sites, sv, gfab, ref_fasta, gaf, sep_pattern, gff3_file, ref_name=None):
    if ref_name is None:
        ref_name = os.path.splitext(os.path.basename(ref_fasta))[0]


    print(f"Running complex with vcf_file={vcf_file}, a={a}, n={n}, s={s}, l={l}, sites={sites}, sv={sv}")
    site_df = find_complex_sites(vcf_file, n, s)
    complex_regions = find_interesting_regions(site_df, l)
 
    pd.DataFrame(complex_regions).to_csv("complex_regions.csv")
    graphs = [gfab]
    print('reading alignments')
    
        
    if not os.path.exists(gaf):
        os.system(f'GraphAligner -g {gfab} -f {gene_seq_fa} -t 2 -a {gaf} -x vg --multimap-score-fraction 0.1')
    
    # %%
    df_arp = pd.read_csv(gaf, sep = '\t', header=None)
    print(df_arp[0].head())

    df_arp[0] = df_arp[0].str.split('_', expand = True)[0].str.split('-', expand=True).iloc[:, 1:].fillna('').apply('-'.join, axis = 1).str.strip('-')
    print(df_arp[0].head())
    # %%
    # df_cpc = pd.read_csv('chm13.genes.cpc.gaf', sep = '\t', header=None)
    # df_cpc[0] = df_cpc[0].str.split('_', expand = True)[0].str.split('-', expand=True).iloc[:, 1:].fillna('').apply('-'.join, axis = 1).str.strip('-')


    # %%
    tasks = []
    for gfab in graphs:
        for chrom,start,end in complex_regions:
        
            # gfab = 'hprc-v2.1-mc-chm13.gfab'
            #gff3_file = "chm13v2.0_RefSeq_Liftoff_v5.1.gff3"
            fasta_file = ref_fasta
 
            region = 'complex'
            cutpoints = 1
            graph_base = gfab.split('.')[0]
            
            
            gene_alignments = gaf
          
            df_all = df_arp
            
            # if gfab.split('.')[0] == 'CPC': 
            #     ref = 'CHM13v2'
            #     gene_alignments = 'chm13.genes.cpc.gaf'
            #     sep = '.'
            #     df_all = df_cpc
            
            
            
            query_region = f'{ref_name}{sep_pattern}{chrom}:{start}-{end}'

            region = f'{chrom}{sep_pattern}{start}_{end}'

            workdir = f'all_plottables/{graph_base}.{region}.wd'
            viz_output = f'{workdir}/{region}.{cutpoints}.{graph_base}.gfa'
            connected_output = f'all_walks_gfas/{region}.{cutpoints}.{graph_base}.walks.gfa'
        
            os.system('mkdir -p all_walks_gfas')
        
            tasks.append((gfab, gff3_file, region, cutpoints, graph_base, viz_output, connected_output, ref_name, chrom, int(start), int(end), query_region, workdir, gene_alignments, df_all, sep_pattern))

    
    for task in tasks:
        produce_plottable(*task)
            




    

def main(args):
    run_panscan_complex(
        args.vcf_file, 
        args.a, 
        args.n, 
        args.s, 
        args.l, 
        args.sites, 
        args.sv, 
        args.gfab_file, 
        args.ref_fasta, 
        args.gaf_file, 
        args.sep_pattern, 
        args.gff3,
        ref_name=args.ref_name
            )

def add_subparser(subparsers):
    parser = subparsers.add_parser("complex", help="Analyze complex loci from a VCF file.")
    parser.add_argument("vcf_file", help="Path to the VCF file generated from Minigraph-cactus pipeline.")
    parser.add_argument("gfab_file", help="Path to the gfab file generated from the GFA file")
    parser.add_argument("--ref_fasta", help="Path to the reference file" )
    parser.add_argument("--gaf_file", help="Path to the gene alignments to the graph in GAF format")
    parser.add_argument("--sep_pattern", help="Separator for the sample and the sequence as present in GAF file (eg. '#')")
    parser.add_argument('--gff3', help='Path to the gff3 file')
    parser.add_argument("-a", type=int, default=5, help="Number of alleles to define a complex site (default: 5).")
    parser.add_argument("-n", type=int, default=1, help="Number of 10kb variants to define a complex site (default: 1).")
    parser.add_argument("-s", type=int, default=10000, help="Minimum size of a variant to define a complex site (default: 10000).")
    parser.add_argument("--regions", action="store_true", help="List complex regions.")
    parser.add_argument("-l", type=int, default=100000, help="Length of the region to define a complex region (default: 100000).")
    parser.add_argument("--sites", type=int, default=1, help="Number of complex sites in a region to define it as complex (default: 1).")
    parser.add_argument("--sv", type=int, default=1, help="Number of secondary SVs in a region to define it as complex (default: 1).")
    parser.add_argument("--ref_name", default=None, help="Reference name to use in region queries as present in GAF file (eg. CHM13)(default: None)")
    parser.set_defaults(func=main)





# %%
