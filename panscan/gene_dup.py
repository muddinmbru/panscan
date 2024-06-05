import argparse
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
from matplotlib_venn import venn3 # type: ignore


def process_and_plot_data(result_file):
    # Load data
    result_df = pd.read_csv(result_file, index_col=0)
    hprc = pd.read_csv('hprc-matrix.csv', skipfooter=1, engine='python', index_col=0)
    cpc = pd.read_csv('cpc-matrix.csv', index_col=0)
    cpc = cpc.drop(columns=['Unnamed: 2', 'Unnamed: 3', 'Unnamed: 4'])
    
    # Plotting duplications
    duplicates_count = (result_df != 0).sum()
    duplicates_count = duplicates_count.sort_values()
    plt.figure(figsize=(40, 24))
    duplicates_count.plot(kind='bar', color='#1F1F6A', edgecolor='black')
    plt.xlabel('Genome by number of duplications')
    plt.ylabel('Duplicated genes per genome')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xticks(rotation=0, color='white', ha='right')
    plt.tight_layout()
    plt.savefig('genome-by-dup.png')
    plt.show()

    
    # Venn diagram
    common_indices = hprc.index.intersection(result_df.index).intersection(cpc.index)
    venn_data = {
        '100': len(hprc.index.difference(result_df.index).difference(cpc.index)),
        '010': len(result_df.index.difference(hprc.index).difference(cpc.index)),
        '001': len(cpc.index.difference(hprc.index).difference(result_df.index)),
        '110': len(hprc.index.intersection(result_df.index).difference(cpc.index)),
        '101': len(hprc.index.difference(result_df.index).intersection(cpc.index)),
        '011': len(result_df.index.difference(hprc.index).intersection(cpc.index)),
        '111': len(common_indices)
    }
    plt.figure(figsize=(40, 40))
    venn = venn3(subsets=venn_data, set_labels=('HPRC', 'APR', 'CPC'))

    # Update font size for the Venn diagram labels
    for text in venn.set_labels:
        if text:  # Check if the label exists
            text.set_fontsize(80)
    if venn.subset_labels:  # Check if subset labels exist
        for text in venn.subset_labels:
            if text:  # Check if the subset label exists
                text.set_fontsize(60)
    plt.savefig('venn-comparison.png')
    plt.show()

    # Frequency comparison for each category
    def frequency_comparison(df1, df2, output_file, label1, label2):
        df1_freq = df1.apply(lambda row: (row != 0).sum() / len(row), axis=1).to_frame(name='Frequency_x')
        df2_freq = df2.apply(lambda row: (row != 0).sum() / len(row), axis=1).to_frame(name='Frequency_y')
        
        comparison_df = df1_freq.merge(df2_freq, left_index=True, right_index=True)
        comparison_df['Absolute_Difference'] = abs(comparison_df['Frequency_x'] - comparison_df['Frequency_y'])
        comparison_df['Percentage_Difference'] = (comparison_df['Absolute_Difference'] / comparison_df[['Frequency_x', 'Frequency_y']].max(axis=1)) * 100
        
        filtered_df = comparison_df[comparison_df['Percentage_Difference'] >= 5]
        sorted_df = filtered_df.sort_values(by='Percentage_Difference', ascending=False)
        top_5_diff = sorted_df.head(5)
        
        plt.figure(figsize=(40, 36))
        plt.barh(top_5_diff.index, top_5_diff['Frequency_x'], color='#1F1F6A', label=label1)
        plt.barh(top_5_diff.index, -top_5_diff['Frequency_y'], color='#EA4E15', label=label2)
        plt.xlabel('CNV frequency', labelpad=50)
        plt.legend()
        plt.rcParams['font.size'] = 70

        plt.tight_layout(pad=4)
        plt.savefig(output_file)
        plt.show()

    # Plot frequency comparison with CPC
    frequency_comparison(result_df, cpc, 'frequency-comparsion-with-cpc.png', 'Result', 'CPC')

    # Plot frequency comparison with HPRC
    frequency_comparison(result_df, hprc, 'frequency-comparsion-with-hprc.png', 'Result', 'HPRC')

def run_gene_dup(csv_file):
    process_and_plot_data(csv_file)

def main(args):
    run_gene_dup(args.csv_file)

def add_subparser(subparsers):
    parser = subparsers.add_parser("gene_dup", help="Detect gene duplications from a CSV file.")
    parser.add_argument("csv_file", help="Path to the CSV file for gene duplication detection.")
    parser.set_defaults(func=main)