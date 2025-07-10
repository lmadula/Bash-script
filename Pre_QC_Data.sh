#Download plink versions on hemera
#Create the directory called bin
mkdir ~/bin

#Then export it to the path
export PATH=$PATH:~/bin/ 

#Download plink version 1 to the above directory (bin)
wget https://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip

#Then unzip it using unzip
unzip plink-1.07-x86_64.zip

#Download the plink version 2 to the above directory (bin)
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240625.zip

#Then unzip it using unzip
unzip plink2_linux_x86_64_20240625.zip

#Check if variants are present in files 
grep '40310434' flipped_GP2_merge_STELLENBOS.bim
grep '40340400' flipped_GP2_merge_STELLENBOS.bim
grep '161350208' flipped_GP2_merge_STELLENBOS.bim
grep '161785820' flipped_GP2_merge_STELLENBOS.bim
grep '161350203' flipped_GP2_merge_STELLENBOS.bim
grep '20639990' flipped_GP2_merge_STELLENBOS.bim

#Use masterkey to add phenotype and sex information for binary files (Upload the file tha Abigail sent to hemera)
awk -F',' '{print $1,$4}' STELLENBOS_key_aug7_2024.csv > pheno_info_mono.txt       

awk 'BEGIN {FS=OFS=" "}
    NR==FNR { phenotype[$1] = $2; next }
    { $6 = phenotype[$2]; print }' pheno_info_mono.txt flipped_GP2_merge_STELLENBOS.fam > updated.fam

mv updated.fam flipped_GP2_merge_STELLENBOS.fam 

awk -F',' '{print $1,$7}' STELLENBOS_key_aug7_2024.csv > sex_info_mono.txt       

awk 'BEGIN {FS=OFS=" "}
    NR==FNR { phenotype[$1] = $2; next }
    { $5 = phenotype[$2]; print }' sex_info_mono.txt flipped_GP2_merge_STELLENBOS.fam > updated.fam

mv updated.fam flipped_GP2_merge_STELLENBOS.fam

#Convert all positions to chr:pos:ref:alt (Test ) 
./plink2 --bfile flipped_GP2_merge_STELLENBOS --set-all-var-ids @:# --make-bed --out gp2_merge_STEL
###Not working - skip
./plink2 --bfile flipped_GP2_merge_STELLENBOS --set-all-var-ids 'chr'@:#":"\$r":"\$a --make-bed --out gp2_variant_recoded
###Not working - skip

# SNP Missingness 5%
./plink --bfile flipped_GP2_merge_STELLENBOS --geno 0.05 --make-bed --out gp2__geno_0.05
#Test at what threshold individual previously removed is included
plink --bfile flipped_GP2_merge_STELLENBOS --geno 0.05 --make-bed --out gp2_geno_0.045 

#Check if variants are present in files
grep '40310434' gp2__geno_0.05.bim
grep '40340400' gp2__geno_0.05.bim
grep '161350208' gp2__geno_0.05.bim
grep '161785820' gp2__geno_0.05.bim
grep '161350203' gp2__geno_0.05.bim
grep '20639990' gp2__geno_0.05.bim
 
#Individual missingness 10%
plink --bfile gp2__geno_0.05 --mind 0.1 --make-bed --out gp2_mind_0.1
 
#Check for sex discrepancies (F less than 0,2 = female and F larger than 0.8 = male)
#Need to check if their reported and observed sex is the same or different
plink --bfile gp2_mind_0.1 --check-sex
# Remove the list of individuals with the status PROBLEM. 
grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}'> sex_discrepancy.txt
plink --bfile gp2_mind_0.1 --remove sex_discrepancy.txt --make-bed --out gp2_sex_check_qc
 
#Delete SNPs not in HWE 
#Controls
plink --bfile gp2_sex_check_qc --hwe 1e-6 --make-bed --out hwe_controls
#Cases
plink --bfile hwe_controls --hwe 1e-10 --hwe-all --make-bed --out gp2_hwe

#Check if variants are present in files
grep '40310434' gp2_hwe.bim
grep '40340400' gp2_hwe.bim
grep '161350208' gp2_hwe.bim
grep '161785820' gp2_hwe.bim
grep '161350203' gp2_hwe.bim
grep '20639990' gp2_hwe.bim
 
#Runs of heterozygosity (LD)
plink --bfile gp2_hwe --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile gp2_hwe --extract indepSNP.prune.in --het --out R_check
 
#To plot the heterozygosity rate distribution (Start here). hao
het <- read.table("R_check.het", head=TRUE)
pdf("heterozygosity.pdf")
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
 
#The following code was run in R and produced a list of individuals who deviated more than 3 standard deviations from the mean
het <- read.table("R_check.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)
dev.off() 
 
#Output of the command above: fail-het-qc.txt .
#Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt 
 
#Remove heterozygosity rate outliers.
plink --bfile gp2_hwe --remove het_fail_ind.txt --make-bed --out gp2_het

#Check if variants are present in files
grep '40310434' gp2_het.bim
grep '40340400' gp2_het.bim
grep '161350208' gp2_het.bim
grep '161785820' gp2_het.bim
grep '161350203' gp2_het.bim
grep '20639990' gp2_het.bim 

#Fix the alleles switches made by PLINK 
#Code to generate the Homo_sapiens_assembly38.fasta file
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta
bwa index Homo_sapiens_assembly38.fasta

#Python script
import os
# Define input binary files
dataIn = "gp2_het"
# Convert binary files to VCF format with PLINK
os.system(f"plink --bfile ./{dataIn} --recode vcf-iid bgz --out ./gp2_fl --output-chr chrMT")
# Index the resulting VCF file with BCFtools
os.system(f"bcftools index ./gp2_fl.vcf.gz")
# Define the path to the reference FASTA file
fasta = f"/home/lmadula/Homo_sapiens_assembly38.fasta"
# Fix the reference alignment and create a new VCF file
os.system(f"bcftools +fixref ./gp2_fl.vcf.gz -Oz -o gp2_fl_1.vcf.gz -- -f {fasta} -m flip -d")
# Re-index the normalized VCF file
os.system(f"bcftools index ./gp2_fl_1.vcf.gz -f")
# Convert the normalized VCF file back to PLINK binary format (BED/BIM/FAM)
os.system(f"plink --vcf ./gp2_fl_1.vcf.gz --make-bed --out ./gp2_fl_binary --double-id")

#Extract 22 PD genes
plink --bfile gp2_het --extract range pd_gene_list.txt --make-bed --out extract_genes 

#Check if variants are present in files
grep '40310434' extract_genes.bim
grep '40340400' extract_genes.bim
grep '161350208' extract_genes.bim
grep '161785820' extract_genes.bim
grep '161350203' extract_genes.bim
grep '20639990' extract_genes.bim (done)

#Convert to .vcf
plink --bfile extract_genes --recode vcf --out extract_genes 

#Split samples in .vcf file using VCF_Split_#3.py 
#!/bin/bash
#PBS -N SplitGeneExtract
#PBS -l select=1:ncpus=24:mpiprocs=4
#PBS -P HEAL1381
#PBS -q serial
#PBS -l walltime=48:00:00
#PBS -e lmadula@hemera.mb.sun.ac.za:/home/lmadula/qc_output/split.err
#PBS -o lmadula@hemera.mb.sun.ac.za:/home/lmadula/qc_output/split.out
#PBS -m abe
#PBS -M lmadula@sun.ac.za
### start here and run it in exucutable. download bcftools of hemera
## How to download bcftools

cd   lmadula@hemera.mb.sun.ac.za:/home/lmadula/NBA_filesfix

## spack load bcftools

for sample in `bcftools query -l /home/lmadula/NBA_filesfix/extract_genes.vcf`
do
        bcftools view -c1 -s ${sample} /home/lmadula/NBA_filesfix/extract_genes.vcf -o /home/lmadula/NBA_filesfix/${sample}.vcf
done

#Check if sample.vcf files contain known variants (inspete the specific individuals)
grep '40310434' {sample}.vcf
#LRRK2 R1441C
grep '40340400' {sample}.vcf
#LRRK2 G2019S
grep '161350208' {sample}.vcf
#PRKN G430D
grep '161785820' {sample}.vcf
#PRKN R275W
grep '161350203' {sample}.vcf
#PRKN M432V
grep '20639990' {sample}.vcf
#PINK1 Y258X

grep '40310434' 0_m-STELLENBOS_000029_s1.vcf
grep '40340400' 0_m-STELLENBOS_000029_s1.vcf
grep '161350208' 0_m-STELLENBOS_000029_s1.vcf
grep '161785820' 0_m-STELLENBOS_000029_s1.vcf
grep '161350203' 0_m-STELLENBOS_000029_s1.vcf
grep '20639990' 0_m-STELLENBOS_000029_s1.vcf


#Check if variants are present in GRCH37 harmonized files (ignore)
grep '40704236' filename
#LRRK2 R1441C
grep '40734202' filename
#LRRK2 G2019S
grep '161771240' filename
#PRKN G430D
grep '162206852' filename
#PRKN R275W
grep '161771235' filename
#PRKN M432V
grep '20966483' filename
#PINK1 Y258X



## Essential information deals with the data manipulation on R programming language
https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html
## Explaining reeluar expressions in simple terms
https://www.freecodecamp.org/news/practical-regex-guide-with-real-life-examples/
### Bioinformatics links that deals with the fixing of the strand and alleles and reference alleles
https://shicheng-guo.github.io/blog/page8/
## Online couse
https://www.futurelearn.com/courses/interpreting-genomic-variation-overcoming-challenges-in-diverse-populations/1/steps/1919948


### code I used to copy the data from the remote server to the local server
rsync -avzhe ssh lmadula@lengau.chpc.ac.za:/mnt/lustre/groups/HEAL1381/gp2_NBA_fixed/flipped_GP2_merge_STELLENBOS.log .



### Installing the VEP on hemera , Initially is was designed in a command line
git clone https://github.com/Ensembl/ensembl-vep.git

### Insall the perl script inside the ensembl-vep
cd ensembl-vep
perl INSTALL.pl

###  VEP can use plugin modules written in Perl to add functionality to the software.Plugins are a powerful way to extend, filter and manipulate the VEP output.
perl INSTALL.pl -a p -g list
https://github.com/Ensembl/VEP_plugins/blob/release/112/CADD.pm
mv CADD.pm ~/.vep/Plugins/

Copy Number Variation (CNV)
Genetic alterations that involve the duplication or deletion of segments of DNA (typically 50 bases or more), resulting in changes in copy numbers of specific regions of the genome

Heterozygous
Refers to a gene or an individual having two different alleles at a genetic locus.

Homozygous
Refers to a gene or an individual having two copies of the same allele at a locus.


Pathogenic variant
A change in the DNA sequence that causes a person to have or be at risk of developing a certain genetic disorder or disease.

Single Nucleotide Variations (SNV)
Single nucleotide variations involve changes at a single nucleotide within the DNA sequence, encompassing substitutions, insertions, or deletions.

Single nucleotide insertion: Addition of a single nucleotide into a DNA or RNA sequence.
Single nucleotide deletion: Removal of a single nucleotide from a DNA or RNA sequence.
Substitution: Replacement of one nucleotide with another in a DNA or RNA sequence. These are the most common variations observed.

Copy Number Variations (CNV)
Copy number variations are genetic alterations that involve the duplication or deletion of segments of DNA (typically 50 bases or more), resulting in changes in copy numbers of specific regions of the genome.

Copy number deletion: Loss of a portion of a DNA segment, resulting in fewer copies of a specific genetic sequence.
Copy number amplification: Increase in the number of copies of a specific genetic sequence within the DNA. This includes duplications, triplications and higher amplifications.