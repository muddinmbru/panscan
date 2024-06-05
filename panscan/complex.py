import argparse
import pandas as pd  # type: ignore


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

def run_panscan_complex(vcf_file, a, n, s, l, sites, sv):
    print(f"Running complex with vcf_file={vcf_file}, a={a}, n={n}, s={s}, l={l}, sites={sites}, sv={sv}")
    site_df = find_complex_sites(vcf_file, n, s)
    complex_regions = find_interesting_regions(site_df, l)
 
    pd.DataFrame(complex_regions).to_csv("complex_regions.csv")




def main(args):
    run_panscan_complex(args.vcf_file, args.a, args.n, args.s, args.l, args.sites, args.sv)

def add_subparser(subparsers):
    parser = subparsers.add_parser("complex", help="Analyze complex loci from a VCF file.")
    parser.add_argument("vcf_file", help="Path to the VCF file generated from Minigraph-cactus pipeline.")
    parser.add_argument("-a", type=int, default=5, help="Number of alleles to define a complex site (default: 5).")
    parser.add_argument("-n", type=int, default=1, help="Number of 10kb variants to define a complex site (default: 1).")
    parser.add_argument("-s", type=int, default=10000, help="Minimum size of a variant to define a complex site (default: 10000).")
    parser.add_argument("--regions", action="store_true", help="List complex regions.")
    parser.add_argument("-l", type=int, default=100000, help="Length of the region to define a complex region (default: 100000).")
    parser.add_argument("--sites", type=int, default=1, help="Number of complex sites in a region to define it as complex (default: 1).")
    parser.add_argument("--sv", type=int, default=1, help="Number of secondary SVs in a region to define it as complex (default: 1).")
    parser.set_defaults(func=main)


if __name__=="__main__":


    vcf_file = "apr_review_v1_2902_chm13.vcf"
    a = 5
    n = 1
    s = 10000
    l = 100000
    sites = 1
    sv = 1

    args = argparse.Namespace(vcf_file=vcf_file, a=a, n=n, s=s, l=l, sites=sites, sv=sv)

    # Call main with the namespace object
    main(args)


