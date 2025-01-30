import subprocess
import sys
import argparse

def run_preprocess_vcf(vcf_file):
    subprocess.run(['perl', 'scripts/preprocessVCF.pl', '-i', vcf_file])

def main(args):
    vcf_file = args.vcf
    print("Running VCF preprocessing on:", vcf_file)
    run_preprocess_vcf(vcf_file)

def add_subparser(subparsers):
    # Add the subcommand 'preprocess_vcf'
    subparser = subparsers.add_parser('preprocess_vcf', help="Preprocess VCF file by normalizing, decomposing, and merging genotypes.")
    subparser.add_argument("vcf", help="Path to the input VCF file.")
    subparser.set_defaults(func=main)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess VCF file by normalizing, decomposing, and merging genotypes.")
    parser.add_argument("vcf", help="Path to the input VCF file.")
    args = parser.parse_args()
    main(args)
