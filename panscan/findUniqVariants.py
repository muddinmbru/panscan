import argparse
import subprocess
import os

def run_find_uniq_variants(vcf, var_type, db, overlap, output):
    # Get the absolute path of the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Change the current working directory to the script directory
    #os.chdir(script_dir)  # Ensures that Perl can find config.yaml in the same dir
    
    # Construct the absolute path to findUniqVariants.pl
    perl_script = os.path.join(script_dir, "scripts", "findUniqVariants.pl")
    perl_script = os.path.abspath(perl_script)

    # Locate perlModules dynamically (assumes it's in the same dir as findUniqVariants.pl)
    perl_modules = os.path.join(script_dir, "scripts")

    # Dynamically find the config.yaml in the same directory as the Perl script
    config_path = os.path.join(script_dir, "scripts", "config.yaml")
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file 'config.yaml' not found in {os.path.dirname(perl_script)}")

    # Debugging prints (optional)
    # print(f"Running Perl script at: {perl_script}")
    # print(f"Using Perl module path: {perl_modules}")
    # print(f"Using config file: {config_path}")

    # Set PERL5LIB to locate Perl modules
    env = os.environ.copy()
    env["PERL5LIB"] = f"{perl_modules}:{env.get('PERL5LIB', '')}"

    # Constructing the command based on input arguments
    cmd = [
        'perl', perl_script, 
        '--i', vcf,           # Input VCF file
        '--t', var_type,      # Variant type (SNP/INDEL/SV)
        '--db', db,           # Databases (comma-separated)
        '--op', str(overlap),  # Overlap value
        '--config', config_path  # Pass the dynamically found config file
    ]
    if output:
        cmd += ['--output', output]  # Optional output directory

    # Run the Perl script
    subprocess.run(cmd, check=True, env=env)

def main(args):
    run_find_uniq_variants(args.vcf_file, args.var_type, args.db, args.overlap, args.output)

def add_subparser(subparsers):
    # Parser for find_uniq_variants command
    parser = subparsers.add_parser("find_uniq_variants", help="Find unique variants in a VCF file.")
    
    # Arguments for the VCF file, variant type, database selection, overlap, and output
    parser.add_argument("--i", dest="vcf_file", required=True, help="Specify the input VCF file. (For SV comparison, input should be a pbsv VCF file)")
    parser.add_argument("--t", dest="var_type", choices=["SNP", "INDEL", "SV"], required=True, help="Specify the variant type to compare (SNP/INDEL/SV).")
    
    parser.add_argument("--db", dest="db", default="ALL", help="""

        Specify the databases (comma-separated for multiple databases).
        Available databases for:
            SNP: dbSNP, gnomAD, 1000Genomes, GME
            InDel: gnomAD, 1000Genomes, GME
            SV Insertion: DGV
            SV Deletion: DGV, 1000Genomes
        Default: ALL
    """)
    
    parser.add_argument("--op", dest="overlap", type=int, default=80, help="Specify the percentage of overlap for SV comparison. Default: 80")
    parser.add_argument("--output", help="Optional output directory to save results.")
    
    # Set the function to call when this subparser is used
    parser.set_defaults(func=main)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find unique variants in a VCF file.")
    subparsers = parser.add_subparsers()
    add_subparser(subparsers)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
