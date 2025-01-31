import subprocess
import os
import argparse

def run_preprocess_vcf(vcf_file):
    # Get the absolute path of the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct the absolute path to preprocessVCF.pl
    perl_script = os.path.join(script_dir, "..", "scripts", "perlModules", "preprocessVCF.pl")
    
    # Normalize the path
    perl_script = os.path.abspath(perl_script)

    # Debugging print (optional)
    print(f"Running Perl script at: {perl_script}")

    # Run the Perl script
    subprocess.run(['perl', perl_script, '-i', vcf_file], check=True)

def main(args):
    vcf_file = args.vcf
    print("Running VCF preprocessing on:", vcf_file)
    run_preprocess_vcf(vcf_file)

def add_subparser(subparsers):
    subparser = subparsers.add_parser('preprocess_vcf', help="Preprocess VCF file by normalizing, decomposing, and merging genotypes.")
    subparser.add_argument("vcf", help="Path to the input VCF file.")
    subparser.set_defaults(func=main)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess VCF file by normalizing, decomposing, and merging genotypes.")
    parser.add_argument("vcf", help="Path to the input VCF file.")
    args = parser.parse_args()
    main(args)
