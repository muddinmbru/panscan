import subprocess
import os
import argparse

def run_novel_seq(vcf1, vcf2):
    # Get the absolute path of the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to findNovelSeq.pl
    perl_script = os.path.join(script_dir, "scripts", "findNovelSeq.pl")
    perl_script = os.path.abspath(perl_script)

    if not os.path.isfile(perl_script):
        raise FileNotFoundError(f"Perl script not found: {perl_script}")

    # Run the Perl script with positional arguments
    cmd = ['perl', perl_script, vcf1, vcf2]
    
    # Run the command
    subprocess.run(cmd, check=True)

def main(args):
    vcf1 = args.vcf1
    vcf2 = args.vcf2

    print("Running novel sequence analysis with:", vcf1, vcf2)
    run_novel_seq(vcf1, vcf2)

def add_subparser(subparsers):
    parser = subparsers.add_parser("novel_seq", help="Detect novel sequences from two VCF files.")
    parser.add_argument("vcf1", help="Path to the first VCF file.")
    parser.add_argument("vcf2", help="Path to the second VCF file.")
    parser.set_defaults(func=main)
