import subprocess
import sys

def run_variant_analysis(vcf):
    subprocess.run(['perl', 'scripts/variant_analysis.pl', vcf])

def main(args):
    run_variant_analysis(args.vcf_file)

def add_subparser(subparsers):
    parser = subparsers.add_parser("variant_analysis", help="Perform variant analysis on a VCF file.")
    parser.add_argument("vcf_file", help="Path to the VCF file for variant analysis.")
    parser.add_argument("--output", help="Output directory to save the results.")
    parser.set_defaults(func=main)



