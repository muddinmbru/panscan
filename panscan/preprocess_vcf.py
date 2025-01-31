import subprocess
import os

def run_preprocess_vcf(vcf_file):
    # Get the absolute path of the current script (preprocess_vcf.py)
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Find preprocessVCF.pl (assumes it's in ../scripts/)
    scripts_dir = os.path.abspath(os.path.join(script_dir, "..", "scripts"))
    perl_script = os.path.join(scripts_dir, "preprocessVCF.pl")

    # Locate perlModules dynamically (assumes it's in the same dir as preprocessVCF.pl)
    perl_modules = os.path.join(scripts_dir, "perlModules")

    # Debugging prints (optional)
    print(f"Running VCF preprocessing on: {vcf_file}")
    print(f"Running Perl script at: {perl_script}")
    print(f"Using Perl module path: {perl_modules}")

    # Run the Perl script
    subprocess.run(['perl', perl_script, '--i', vcf_file], check=True, env=env)

def main(args):
    run_preprocess_vcf(args.vcf)

def add_subparser(subparsers):
    # Add the subcommand 'preprocess_vcf'
    subparser = subparsers.add_parser('preprocess_vcf', help="Preprocess VCF file by normalizing, decomposing, and merging genotypes.")
    subparser.add_argument("--i", dest="vcf", required=True, help="Path to the input VCF file.")
    subparser.set_defaults(func=main)
