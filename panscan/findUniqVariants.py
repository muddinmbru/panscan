import subprocess
import os

def run_find_uniq_variants(vcf, var_type, db, overlap, output):
    # Get the absolute path of the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct the absolute path to findUniqVariants.pl
    perl_script = os.path.join(script_dir, "..", "scripts", "perlModules", "findUniqVariants.pl")
    
    # Normalize the path
    perl_script = os.path.abspath(perl_script)

    # Debugging print (optional)
    print(f"Running Perl script at: {perl_script}")

    # Constructing the command based on input arguments
    cmd = [
        'perl', perl_script, 
        '-i', vcf,           # Input VCF file
        '-t', var_type,      # Variant type (SNP/INDEL/SV)
        '-db', db,           # Databases (ALL or specific ones)
        '-o', str(overlap),  # Overlap value
    ]
    if output:
        cmd += ['--output', output]  # Optional output directory

    # Run the Perl script
    subprocess.run(cmd, check=True)

def main(args):
    run_find_uniq_variants(args.vcf_file, args.var_type, args.db, args.overlap, args.output)

def add_subparser(subparsers):
    # Parser for find_uniq_variants command
    parser = subparsers.add_parser("find_uniq_variants", help="Find unique variants in a VCF file.")
    
    # Arguments for the VCF file, variant type, database selection, overlap, and output
    parser.add_argument("vcf_file", help="Path to the VCF file.")
    parser.add_argument("-t", "--var_type", choices=["SNP", "INDEL", "SV"], required=True, help="Variant type to compare (SNP/INDEL/SV).")
    parser.add_argument("-db", "--db", default="ALL", help="Databases to compare against (default: 'ALL').")
    parser.add_argument("-o", "--overlap", type=int, default=80, help="Overlap percentage (default: 80).")
    parser.add_argument("--output", help="Optional output directory to save results.")
    
    # Set the function to call when this subparser is used
    parser.set_defaults(func=main)
