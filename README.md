
# Panscan
Pangenome analyses toolkit

## Requirements
Requires Python >= 3.7

Matplotlib
perl
pandas
[GFABase](https://github.com/mlin/gfabase) must be installed and be accessible in the shell.
[Bandage](https://github.com/rrwick/Bandage) should be installed and added to the path for plotting

## Installation
To use Panscan, clone the respository and add it to your path.

## Complex loci analyses 
The complex loci can be anlayzed from the vcf file generated from [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 
The complex sites, as defined in the Arab Pangenome paper, are the sites with more than 5 alleles and atleast 1 10kb variant. To list all of the complex sites with the preceding definition run
```
panscan complex --list -a 5 -n1 -s 10000
```

use ```-a```, ```-n```, ```-s``` to redefine the parameters for a complex site. 

### Complex regions
Complex regions are regions of 100Kb with atleast one complex site and another SV. To list complex regions use the --regions flag with
```panscan complex``` subcommand
```-l``` can be used to define the length of the region
```--sites``` can be used to define how many sites should be present
```--sv``` can be used to define how many secondary SVs should be present

for a region to be considered complex.

### End-to-end
Run the ```panscan complex``` to run the command in full to produce complex regions and haplotype walks for each sample in all the regions. 

## Gene-duplication analyses


It takes in your gene duplication matrix, and visualizes the duplications in your data and compares them with the hprc and cpc duplications as well. The plots made are :
Duplications per assembly
Venn diagram of duplications w.r.t hprc and cpc
Frequency comparison of your duplications w.r.t. hprc and cpc . Plots the most distinct ones

To run the gene duplication pipeline, use the ```panscan gene_dup``` subcommand.

## Variant analyses

The program will convert multi-allelic VCF records to single-allelic ones. Next, complex indels will be decomposed into SNPs and indels using the RTG tools "decompose" program. Finally, the genotypes of variants at the same locus will be merged to produce the final pre-processed VCF file.

```panscan variant_analysis``` should be used to process the variants.

## Novel seq
The tool identifies novel sequences present in VCF file1 by comparing SV insertions with those in VCF file2 and reports them in FASTA format. Initially, the input VCF files undergo pre-processing, which involves splitting multi-allelic variants into single-allelic ones and decomposing complex variants into indels and SNPs using the "decompose" program from RTG Tools. After pre-processing, the SV insertions in the VCF files are compared using the "truvari bench" command, identifying novel SV insertions in VCF file1. These novel insertions at the same locus are clustered using the CD-HIT program, and the final novel sequence FASTA file is generated.

```panscan novel_seq``` should be used to process the novel sequences.
