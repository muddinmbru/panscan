
# Panscan
Pangenome analyses toolkit.

## Requirements
- Python >= 3.7 
- perl
- Matplotlib
- pandas
- liftoff 
- [cd-hit](https://github.com/weizhongli/cdhit)
- [rtg-tools](https://github.com/RealTimeGenomics/rtg-tools)
- [GFABase](https://github.com/mlin/gfabase)



## Installation
To use Panscan, clone the respository and add it to your path.

After cloning, enter the repository and run the comman 

```pip install .```

After successful installation run the tool with

``panscan```


## Complex loci analyses 
The complex loci can be anlayzed from the vcf file generated from [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 
The complex sites, as defined in the Arab Pangenome paper, are the sites with more than 5 alleles and atleast one 10kb variant. To list all of the complex sites with the preceding definition run
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

Example command: ```panscan complex --ref_fasta chm13v2.0.fa --gaf_file chm13_mapped_genes.gaf --sep_pattern '#0#' --gff3 chm13v2.0_RefSeq_Liftoff_v5.1.gff3 -a 5 -n 1 -s 10000 --regions -l 100000 --sites 1 --sv 1 --ref_name CHM13 panscan.vcf panscan.gfab ```

The gaf files needed for the complex command can be produced by aligning the gene sequences file to your pangenome. The gene sequence files and a script to produce these files are present in the 'complex' directory.

### End-to-end
Run the ```panscan complex``` to run the command in full to produce complex regions and haplotype walks for each sample in all the regions. 


## Gene-duplication analyses

There are 2 parts to this:
You have to first run ```panscan make-dup-mtx``` to produce the gene-duplication matrix from all your assmeblies.

Then the second comman ```panscan gene-dup``` takes in your gene duplication matrix, and visualizes the duplications in your data and compares them with the hprc and cpc duplications as well. The plots made are :
Duplications per assembly
Venn diagram of duplications w.r.t hprc and cpc
Frequency comparison of your duplications w.r.t. hprc and cpc . Plots the most distinct ones


## Preprocess Vcf

The program will convert multi-allelic VCF records to single-allelic ones. Next, complex indels will be decomposed into SNPs and indels using the RTG tools "decompose" program. Finally, the genotypes of variants at the same locus will be merged to produce the final pre-processed VCF file.

```panscan preprocess_vcf``` should be used to process the variants.

## Novel seq
The tool identifies novel sequences present in VCF file1 by comparing SV insertions with those in VCF file2 and reports them in FASTA format. Initially, the input VCF files undergo pre-processing, which involves splitting multi-allelic variants into single-allelic ones and decomposing complex variants into indels and SNPs using the "decompose" program from RTG Tools. After pre-processing, the SV insertions in the VCF files are compared using the "truvari bench" command, identifying novel SV insertions in VCF file1. These novel insertions at the same locus are clustered using the CD-HIT program, and the final novel sequence FASTA file is generated.

```panscan novel_seq``` should be used to process the novel sequences.
