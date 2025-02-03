import subprocess
import os
import sys

def verify_perl_environment(perl_script, perl_modules):
    # Get the current directory of the Python script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to preprocessVCF.pl dynamically
    perl_script_path = os.path.join(script_dir , 'scripts', 'preprocessVCF.pl')

    if not os.path.isfile(perl_script_path):
        raise FileNotFoundError(f"Perl script not found: {perl_script_path}")
    
    # Now, use this path for subprocess call
    subprocess.run(['perl', perl_script_path, '--i', 'Sample.vcf'])

def run_preprocess_vcf(vcf_file):
    # Get the absolute path of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Find preprocessVCF.pl and perlModules
    scripts_dir = os.path.abspath(os.path.join(script_dir, "scripts"))
    perl_script = os.path.join(scripts_dir, "preprocessVCF.pl")
    perl_modules = scripts_dir
    
    try:
        # Verify environment before proceeding
        verify_perl_environment(perl_script, perl_modules)
        
        #print(f"Running VCF preprocessing on: {vcf_file}")
        #print(f"Running Perl script at: {perl_script}")
        #print(f"Using Perl module path: {perl_modules}")
        
        # Set PERL5LIB - make sure it's an absolute path
        env = os.environ.copy()
        env["PERL5LIB"] = f"{perl_modules}:{env.get('PERL5LIB', '')}"
        
        # Run with more detailed error output
        result = subprocess.run(
            ['perl', perl_script, '--i', vcf_file],
            check=True,
            env=env,
            capture_output=True,
            text=True
        )
        
        # Print stdout if any
        if result.stdout:
            print(result.stdout)
            
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        raise
    except subprocess.CalledProcessError as e:
        print(f"Error running Perl script:", file=sys.stderr)
        print(f"Exit code: {e.returncode}", file=sys.stderr)
        print(f"Error output:", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        raise

def main(args):
    run_preprocess_vcf(args.vcf)

def add_subparser(subparsers):
    subparser = subparsers.add_parser('preprocess_vcf', 
        help="Preprocess VCF file by normalizing, decomposing, and merging genotypes.")
    subparser.add_argument("--i", dest="vcf", required=True, 
        help="Path to the input VCF file.")
    subparser.set_defaults(func=main)
