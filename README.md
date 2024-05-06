
# Panscan
Pangenome analyses toolkit

## Requirements
Requires Python >= 3.7

[GFABase](https://github.com/mlin/gfabase) must be installed and be accessible in the shell.
[Bandage](https://github.com/rrwick/Bandage) should be installed and added to the path for plotting.

## Installation
To use Panscan, clone the respository and add it to your path.

## Complex loci analyses 
The complex loci can be anlayzed from the vcf file generated from [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 
The complex sites, as defined in the Arab Pangenome paper, are the sites with more than 5 alleles and atleast 1 10kb variant. To list all of the complex sites with the preceding definition run
panscan complex --list -a 5 -n1 -s 10000

use -a, -n, -s to redefine the parameters for a complex site. 

### Complex regions
Complex regions are regions of 100Kb with atleast one complex site and another SV. To list complex regions use the --regions flag with
**panscan complex** subcommand
-l can be used to define the length of the region
--sites can be used to define how many sites should be present
--sv can be used to define how many secondary SVs should be present

for a region to be considered complex.

### End-to-end
Run the **panscan complex** to run the command in full to produce complex regions and haplotype walks for each sample in all the regions. 



## Variant analyses

