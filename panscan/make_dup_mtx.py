import os
import glob
import subprocess
import pandas as pd
import sys

def process_assemblies(gencode_gff3, fofn_path, ref_fa, threads):
    # Step 1: Set up the main folder and directories
    os.makedirs('copy-num', exist_ok=True)

    # Read assembly paths from the file of files (fofn_path)
    with open(fofn_path, 'r') as f:
        assembly_paths = [line.strip() for line in f if line.strip()]

    # Step 2: Process each assembly sequentially
    for assembly_path in assembly_paths:
        assembly_name = os.path.splitext(os.path.basename(assembly_path))[0]
        assembly_dir = os.path.join('copy-num', assembly_name)
        os.makedirs(assembly_dir, exist_ok=True)
        gencode_symlink = os.path.join(assembly_dir, 'gencode.gff3')
        if not os.path.islink(gencode_symlink):
            os.symlink(gencode_gff3, gencode_symlink)
        liftoff_cmd = (
            f"liftoff -p {threads} -sc 0.90 -copies -g {gencode_symlink} "
            f"-u {os.path.join(assembly_dir, assembly_name)}.unmapped.txt "
            f"-o {os.path.join(assembly_dir, assembly_name)}.genedup.gff3 "
            f"-polish {assembly_path} {ref_fa} "
            f"> {os.path.join(assembly_dir, 'liftoff.log')} 2>&1"
        )
        subprocess.run(liftoff_cmd, shell=True, check=True)

    # Modified file reading and merging section
    gff3_files = glob.glob('copy-num/**/*.genedup.gff3', recursive=True)
    gff_data_list = []
    
    for file_path in gff3_files:
        sample_name = os.path.basename(file_path).rsplit('.', 1)[0]
        # Read only lines that don't start with #
        with open(file_path, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]

        # Convert lines to DataFrame
        gff_data = pd.read_csv(
            pd.io.common.StringIO('\n'.join(lines)),
            sep='\t',
            header=None,
            dtype=str  # Read all columns as strings initially
        )
        # Filter for gene entries
        gff_data = gff_data[gff_data[2] == 'gene']
        gff_data['sample'] = sample_name
        gff_data_list.append(gff_data)
        print(f"Processed {file_path}")

    # Combine all the DataFrames
    combined_gff_data = pd.concat(gff_data_list, ignore_index=True)
    combined_gff_data.to_csv('all_dup_data_gene.csv', index=False)
    
    # Rest of your original code
    combined_gff_data['cn'] = combined_gff_data[8].str.extract(r'extra_copy_number=([^;]+)').astype(float)
    combined_gff_data['gene_id'] = combined_gff_data[8].str.extract(r'gene_id=([^;]+)')
    combined_gff_data['gene_id'] = combined_gff_data['gene_id'].str.split('_').str[0]
    combined_gff_data['gene_name'] = combined_gff_data[8].str.extract(r'gene_name=([^;]+)')
    combined_gff_data['gene_type'] = combined_gff_data[8].str.extract(r'gene_type=([^;]+)')
    
    filtered_data = combined_gff_data[combined_gff_data['cn'] >= 1]
    matrix = filtered_data.pivot_table(index='gene_name', columns='sample', values='cn', aggfunc='max').fillna(0)
    matrix.to_csv('gene-dup-matrix.csv')

def add_subparser(subparsers):
    """Add subparser for the make_dup_mtx command."""
    parser = subparsers.add_parser(
        "make_dup_mtx",
        help="Generate a duplication matrix based on gene duplication data."
    )
    parser.add_argument("gencode_gff3", help="Path to the GENCODE GFF3 file.")
    parser.add_argument("fofn_path", help="Path to the file of filenames (FOFN).")
    parser.add_argument("ref_fa", help="Path to the reference FASTA file.")
    parser.add_argument("threads", type=int, help="Number of threads to use.")
    parser.set_defaults(func=run_make_dup_mtx)

def run_make_dup_mtx(args):
    """Run the make-dup-mtx process with the provided arguments."""
    process_assemblies(args.gencode_gff3, args.fofn_path, args.ref_fa, args.threads)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python make_dup_mtx.py <gencode_gff3> <fofn_path> <ref_fa> <threads>")
        sys.exit(1)
    
    gencode_gff3 = sys.argv[1]
    fofn_path = sys.argv[2]
    ref_fa = sys.argv[3]
    threads = int(sys.argv[4])
    process_assemblies(gencode_gff3, fofn_path, ref_fa, threads)
