#!/bin/bash

# remove duplicate/low quality individuals
# cd /home/soleary/FLOUNDER/VCF_Filtering/
# vcftools --vcf TotalRawSNPs.vcf --out RawSNPs --remove removeInd.txt --recode --recode-INFO-all

# make sure filename for Filter 0 is correct depending on whether or not pre-filtering individuals

# FILTER 0: SEQUENCE QUALITY
vcftools --vcf TotalRawSNPs.vcf --out SFL.F0 --minQ 20 --recode --recode-INFO-all
vcftools --vcf TotalRawSNPs.vcf --out SFL.F0a --minQ 10 --recode --recode-INFO-all
vcftools --vcf TotalRawSNPs.vcf --out SFL.F0b --minQ 15 --recode --recode-INFO-all

# FILTER 1: GENOTYPE CALL RATE & MINIMUM ALLELE COUNT
vcftools --vcf SFL.F0.recode.vcf --out A --geno 0.3 --mac 3 --recode --recode-INFO-all
vcftools --vcf SFL.F0.recode.vcf --out B0 --geno 0.5 --mac 5 --recode --recode-INFO-all
vcftools --vcf SFL.F0.recode.vcf --out B --geno 0.5 --mac 3 --recode --recode-INFO-all

# FILTER 2: ELIMINATE SITES BY REQUIRED MINIMUM DEPTH
vcftools --vcf A.recode.vcf --out A.1 --minDP 3 --recode --recode-INFO-all
vcftools --vcf A.recode.vcf --out A.2 --minDP 5 --recode --recode-INFO-all
vcftools --vcf A.recode.vcf --out A.3 --minDP 10 --recode --recode-INFO-all

vcftools --vcf B.recode.vcf --out B.1 --minDP 3 --recode --recode-INFO-all
vcftools --vcf B.recode.vcf --out B.2 --minDP 5 --recode --recode-INFO-all
vcftools --vcf B.recode.vcf --out B.3 --minDP 10 --recode --recode-INFO-all

# FILTER 3: ELIMINATE INDIVIDUALS/LOCI WITH HIGH PROPORTION OF MISSING DATA I
# A: filter loci
vcftools --vcf A.1.recode.vcf --out A.1a --geno 0.5 --recode --recode-INFO-all
vcftools --vcf A.1a.recode.vcf --out A_1a --missing
vcftools --vcf A.2.recode.vcf --out A.2a --geno 0.5 --recode --recode-INFO-all
vcftools --vcf A.2a.recode.vcf --out A_2a --missing
vcftools --vcf A.3.recode.vcf --out A.3a --geno 0.5 --recode --recode-INFO-all

vcftools --vcf B.1.recode.vcf --out B.1a --geno 0.6 --recode --recode-INFO-all
vcftools --vcf B.1a.recode.vcf --out B_1a --missing
vcftools --vcf B.2.recode.vcf --out B.2a --geno 0.6 --recode --recode-INFO-all
vcftools --vcf B.2a.recode.vcf --out B_2a --missing

# B: filter individuals
mawk -v x=0.99 '$5 > x' A_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1
mawk -v x=0.95 '$5 > x' A_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.2
mawk -v x=0.9 '$5 > x' A_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.3

mawk -v x=0.9 '$5 > x' A_2a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.2a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.2.1

mawk -v x=0.99 '$5 > x' B_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.1
mawk -v x=0.95 '$5 > x' B_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.2
mawk -v x=0.9 '$5 > x' B_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3

# FILTER 4: MINOR ALLELE FREQUENCY AND MEAN DEPTH I
vcftools --vcf A.1.1.recode.vcf --out A.1.1.1 --maf 0.01 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf A.1.1.recode.vcf --out A.1.1.2 --maf 0.02 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf A.1.1.recode.vcf --out A.1.1.3 --maf 0.01 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf A.1.1.recode.vcf --out A.1.1.4 --maf 0.02 --min-meanDP 10 --recode --recode-INFO-all

vcftools --vcf B.1.3.recode.vcf --out B.1.3.1 --maf 0.01 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.2 --maf 0.01 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.3 --maf 0.02 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.4 --maf 0.02 --min-meanDP 10 --recode --recode-INFO-all

# FILTER 5: ELIMINATE INDIVIDUALS/LOCI WITH HIGH PROPORTION OF MISSING DATA II
# A: filter loci
vcftools --vcf A.1.1.1.recode.vcf --out A.1.1.1.1 --geno 0.75 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.1.recode.vcf --out A_1.1.1.1 --missing
vcftools --vcf A.1.1.1.recode.vcf --out A.1.1.1.2 --geno 0.8 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.recode.vcf --out A_1.1.1.2 --missing

vcftools --vcf B.1.3.3.recode.vcf --out B.1.3.3.1 --geno 0.75 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.1.recode.vcf --out B_1.3.3.1 --missing
vcftools --vcf B.1.3.3.recode.vcf --out B.1.3.3.2 --geno 0.8 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.recode.vcf --out B_1.3.3.2 --missing
vcftools --vcf B.1.3.3.recode.vcf --out B.1.3.3.3 --geno 0.85 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.3.recode.vcf --out B_1.3.3.3 --missing

# B: filter individuals
mawk -v x=0.95 '$5 > x' A_1.1.1.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.1.1

mawk -v x=0.9 '$5 > x' A_1.1.1.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.1
mawk -v x=0.95 '$5 > x' A_1.1.1.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2

mawk -v x=0.85 '$5 > x' B_1.3.3.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.1
mawk -v x=0.8 '$5 > x' B_1.3.3.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.2
mawk -v x=0.75 '$5 > x' B_1.3.3.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3

mawk -v x=0.85 '$5 > x' B_1.3.3.3.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.3.1
mawk -v x=0.8 '$5 > x' B_1.3.3.3.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.3.2
mawk -v x=0.75 '$5 > x' B_1.3.3.3.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.3.3

# FILTER 6: MINOR ALLELE FREQUENCY AND MEAN DEPTH II
vcftools --vcf A.1.1.1.1.1.recode.vcf --out A.1.1.1.1.1.1 --maf 0.03 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.1.1.recode.vcf --out A.1.1.1.1.1.2 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all

vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.1 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.2 --maf 0.05 --min-meanDP 15 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.3 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.4 --maf 0.03 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.5 --maf 0.03 --min-meanDP 15 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.6 --maf 0.03 --min-meanDP 20 --recode --recode-INFO-all

vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.1 --maf 0.03 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.2 --maf 0.03 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.3 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.4 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all

# FILTER 7: ELIMINATE INDIVIDUALS/LOCI WITH HIGH PROPORTION OF MISSING DATA III
# A: filter loci
vcftools --vcf A.1.1.1.2.2.3.recode.vcf --out A.1.1.1.2.2.3.1 --geno 0.85 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.1.recode.vcf --out A_1.1.1.2.2.3.1 --missing
vcftools --vcf A.1.1.1.2.2.3.recode.vcf --out A.1.1.1.2.2.3.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.2.recode.vcf --out A_1.1.1.2.2.3.2 --missing

vcftools --vcf B.1.3.3.2.3.2.recode.vcf --out B.1.3.3.2.3.2.1 --geno 0.85 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --out B_1.3.3.2.3.2.1 --missing
vcftools --vcf B.1.3.3.2.3.2.recode.vcf --out B.1.3.3.2.3.2.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --out B_1.3.3.2.3.2.2 --missing

vcftools --vcf B.1.3.3.2.3.4.recode.vcf --out B.1.3.3.2.3.4.1 --geno 0.85 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.1.recode.vcf --out B_1.3.3.2.3.4.1 --missing
vcftools --vcf B.1.3.3.2.3.4.recode.vcf --out B.1.3.3.2.3.4.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.2.recode.vcf --out B_1.3.3.2.3.4.2 --missing

# B: filter individuals
mawk -v x=0.85 '$5 > x' A_1.1.1.2.2.3.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.2.3.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.1
mawk -v x=0.8 '$5 > x' A_1.1.1.2.2.3.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.2.3.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.2

mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.1
mawk -v x=0.65 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.2
mawk -v x=0.6 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.3
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.4

mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.1
mawk -v x=0.65 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.2
mawk -v x=0.6 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.3
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.4

mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.4.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.1
mawk -v x=0.65 '$5 > x' B_1.3.3.2.3.4.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.2

mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.4.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.1
mawk -v x=0.65 '$5 > x' B.1.3.3.2.3.4.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.2

# dDocent filters (requires user input) - will stop script
dDocent_filters B.1.3.3.2.3.2.1.2.recode.vcf B.1.3.3.2.3.2.1.2
dDocent_filters B.1.3.3.2.3.2.2.1.recode.vcf B.1.3.3.2.3.2.2.2
dDocent_filters A.1.1.1.2.2.3.1.2.recode.vcf A.1.1.1.2.2.3.1.2
dDocent_filters B.1.3.3.2.3.4.1.2.recode.vcf B.1.3.3.2.3.4.1.2
dDocent_filters B.1.3.3.2.3.4.2.1.recode.vcf B.1.3.3.2.3.4.2.1

# Remove INDELS
vcfallelicprimitives B.1.3.3.2.3.2.1.2.recode.vcf --keep-info --keep-geno > SFL.prim.vcf
vcftools --vcf SFL.prim.vcf --out B.1.3.3.2.3.2.1.2.SNP --remove-indels --recode --recode-INFO-all

vcfallelicprimitives B.1.3.3.2.3.2.2.2.recode.vcf --keep-info --keep-geno > SFL.prim.vcf
vcftools --vcf SFL.prim.vcf --out B.1.3.3.2.3.2.2.2.SNP --remove-indels --recode --recode-INFO-all

vcfallelicprimitives A.1.1.1.2.2.3.1.2.recode.vcf --keep-info --keep-geno > SFL.prim.vcf
vcftools --vcf SFL.prim.vcf --out A.1.1.1.2.2.3.1.2.SNP --remove-indels --recode --recode-INFO-all

vcfallelicprimitives B.1.3.3.2.3.4.1.2.recode.vcf --keep-info --keep-geno > SFL.prim.vcf
vcftools --vcf SFL.prim.vcf --out B.1.3.3.2.3.4.1.2.SNP --remove-indels --recode --recode-INFO-all

vcfallelicprimitives B.1.3.3.2.3.4.2.1.recode.vcf --keep-info --keep-geno > SFL.prim.vcf
vcftools --vcf SFL.prim.vcf --out B.1.3.3.2.3.4.2.1.SNP --remove-indels --recode --recode-INFO-all

# Final Threshold values
# Filter MAF
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.recode.vcf --out B.1.3.3.2.3.2.1.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.2.2.SNP.recode.vcf --out B.1.3.3.2.3.2.2.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.1.2.SNP.recode.vcf --out A.1.1.1.2.2.3.1.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.1.2.SNP.recode.vcf --out B.1.3.3.2.3.4.1.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.2.1.SNP.recode.vcf --out B.1.3.3.2.3.4.2.1.SNP.final0 --maf 0.05 --recode --recode-INFO-all

# filter missing loci
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.1.2.finala.1 --geno 0.95 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.1.recode.vcf --out B_1.3.3.2.3.2.1.2.finala.1 --missing

vcftools --vcf B.1.3.3.2.3.2.2.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.2.2.finala.1 --geno 0.95 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.1.recode.vcf --out B_1.3.3.2.3.2.2.2.finala.1 --missing

vcftools --vcf A.1.1.1.2.2.3.1.2.SNP.final0.recode.vcf --out A.1.1.1.2.2.3.1.2.finala.1 --geno 0.95 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.1.recode.vcf --out A_1.1.1.2.2.3.1.2.finala.1 --missing

vcftools --vcf B.1.3.3.2.3.4.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.1.2.finala.1 --geno 0.95 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.1.recode.vcf --out B_1.3.3.2.3.4.1.2.finala.1 --missing

vcftools --vcf B.1.3.3.2.3.4.2.1.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.2.1.finala.1 --geno 0.95 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.1.recode.vcf --out B_1.3.3.2.3.4.2.1.finala.1 --missing

vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.1.2.finala.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.2.recode.vcf --out B_1.3.3.2.3.2.1.2.finala.2 --missing

vcftools --vcf B.1.3.3.2.3.2.2.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.2.2.finala.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.2.recode.vcf --out B_1.3.3.2.3.2.2.2.finala.2 --missing

vcftools --vcf A.1.1.1.2.2.3.1.2.SNP.final0.recode.vcf --out A.1.1.1.2.2.3.1.2.finala.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.2.recode.vcf --out A_1.1.1.2.2.3.1.2.finala.2 --missing

vcftools --vcf B.1.3.3.2.3.4.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.1.2.finala.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.2.recode.vcf --out B_1.3.3.2.3.4.1.2.finala.2 --missing

vcftools --vcf B.1.3.3.2.3.4.2.1.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.2.1.finala.2 --geno 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.2.recode.vcf --out B_1.3.3.2.3.4.2.1.finala.2 --missing

# filter missing individuals
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.1.2.finala.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.2.SNP.finalb.1

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.2.2.finala.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.2.SNP.finalb.1

mawk -v x=0.5 '$5 > x' A_1.1.1.2.2.3.1.2.finala.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.2.SNP.finalb.1

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.1.2.finala.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.2.SNP.finalb.1

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.2.1.finala.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.1.SNP.finalb.1

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.1.2.finala.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.2.SNP.finalb.2

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.2.2.finala.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.2.SNP.finalb.2

mawk -v x=0.5 '$5 > x' A_1.1.1.2.2.3.1.2.finala.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.2.SNP.finalb.2

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.1.2.finala.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.2.SNP.finalb.2

mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.2.1.finala.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.1.SNP.finalb.2

# remove unnecessary files
rm *.vcf.idx *.count *.DEPTH *.depth meandepthpersite *.frq *.imiss lowDP.indv *.ldepth *.lmiss *.lowQDloci *.qual A_*.log B_*.log

# create folder with all VCF-files
mkdir vcf
mv *.vcf vcf/

# create folder with all log-files
mkdir Filter.logs
mv *.log Filter.logs/

# delete rmaining files
rm B.* A.*

# generate filter stats (SNP, Contig, Indv)

cd vcf

echo "FILTER SNP CONTIG INDV" > Final_Filter.count

for i in *.vcf
do
  SNP=$(grep -cv '#' $i)
  CONTIG=$(grep -v '#' $i | cut -f 1 | sort | uniq | wc -l)
  INDV=$(vcfsamplenames $i | wc -l)
  echo "$i $SNP $CONTIG $INDV" >> Final_Filter.count
done

# SOL
