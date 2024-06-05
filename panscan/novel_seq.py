import subprocess
import sys
import argparse

def run_novel_seq(vcf1, vcf2):
    subprocess.run(['perl', 'scripts/findNovelSeq.pl', vcf1, vcf2])

def main(args):


    vcf1 = args.vcf1
    vcf2 = args.vcf2

    print("Running novel sequence analysis with: ", vcf1, vcf2)
    run_novel_seq(vcf1, vcf2)

def add_subparser(subparsers):
    parser = subparsers.add_parser("novel_seq", help="Detect novel sequences from two VCF files.")
    parser.add_argument("vcf1", help="Path to the first VCF file.")
    parser.add_argument("vcf2", help="Path to the second VCF file.")
    parser.set_defaults(func=main)

