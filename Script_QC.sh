
#Create the directory called bin
mkdir ~/bin

#Then export it to the path
export PATH=$PATH:~/bin/ 

#Downolad the plink version 1 on the above directory (bin)
wget https://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip


#Then unzip it using the below code
unzip plink-1.07-x86_64.zip

#Download the plink version 2 on the same directory (bin)
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240625.zip

#Then unzip it using the below code
unzip plink2_linux_x86_64_20240625.zip

#Start doing QC (First check  if variants are present)
#Check if variants are present in files 
grep '40310434' flipped_GP2_merge_STELLENBOS.bim
grep '40340400' flipped_GP2_merge_STELLENBOS.bim
grep '161350208' flipped_GP2_merge_STELLENBOS.bim
grep '161785820' flipped_GP2_merge_STELLENBOS.bim
grep '161350203' flipped_GP2_merge_STELLENBOS.bim
grep '20639990' flipped_GP2_merge_STELLENBOS.bim

#Use masterkey to add phenotype and sex information for binary files 
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

#SNP Missingness 5%
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
 
#Check for sex discrepancies (F less than 0,5 = female and F larger than 0.8 = male)
#Need to check if their reported and observed sex is the same or different
plink --bfile flipped_GP2_merge_STELLENBOS --check-sex 0.5 0.8
#Remove the list of individuals with the status PROBLEM. 
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
 
#To plot the heterozygosity rate distribution. 
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

#Fix the reference.
#It is essential to run the code to fix the allele switches made by plink (The script is written in python)
#First download the human reference
#Pre-processing
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta
bwa index Homo_sapiens_assembly38.fasta

#Python script
import os
#Define input binary files
dataIn = "gp2_het"
#Convert binary files to VCF format with PLINK
os.system(f"plink --bfile ./{dataIn} --recode vcf-iid bgz --out ./gp2_fl --output-chr chrMT")
#Index the resulting VCF file with BCFtools
os.system(f"bcftools index ./gp2_fl.vcf.gz")
#Define the path to the reference FASTA file
fasta = f"/home/lmadula/Homo_sapiens_assembly38.fasta"
#Fix the reference alignment and create a new VCF file
os.system(f"bcftools +fixref ./gp2_fl.vcf.gz -Oz -o gp2_fl_1.vcf.gz -- -f {fasta} -m flip -d")
#Re-index the normalized VCF file
os.system(f"bcftools index ./gp2_fl_1.vcf.gz -f")
#Convert the normalized VCF file back to PLINK binary format (BED/BIM/FAM)
os.system(f"plink --vcf ./gp2_fl_1.vcf.gz --make-bed --out ./gp2_fl_binary")

#Check if variants are present in files
grep '40310434' gp2_fl_binary.bim
grep '40340400' gp2_fl_binary.bim
grep '161350208' gp2_fl_binary.bim
grep '161785820' gp2_fl_binary.bim
grep '161350203' gp2_fl_binary.bim
grep '20639990' gp2_fl_binary.bim 

#Extract 22 PD genes
plink --bfile gp2_fl_binary --extract range pd_gene_list.txt --make-bed --out extract_genes 

#Check if variants are present in files
grep '40310434' extract_genes.bim
grep '40340400' extract_genes.bim
grep '161350208' extract_genes.bim
grep '161785820' extract_genes.bim
grep '161350203' extract_genes.bim
grep '20639990' extract_genes.bim 

#Convert to .vcf
plink --bfile extract_genes --recode vcf --out extract_genes

#Split samples in .vcf file using VCF_Split_#3.py
#!/bin/bash
#PBS -N SplitGeneExtract
#PBS -l select=1:ncpus=24:mpiprocs=4
#PBS -P HEAL1381
#PBS -q serial
#PBS -l walltime=48:00:00
#PBS -e lmadula@hemera.mb.sun.ac.za:/home/lmadula/NBA_filesfix/split.err
#PBS -o lmadula@hemera.mb.sun.ac.za:/home/lmadula/NBA_filesfix/split.out
#PBS -m abe
#PBS -M lmadula@sun.ac.za


cd   lmadula@hemera.mb.sun.ac.za:/home/lmadula/NBA_filesfix

#Load the bcftools
Spack load bcftools

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

#Install the annotation tool (ANNOVAR)
#First register and you will receive the link in shortly
curl https://eur03.safelinks.protection.outlook.com/?url=http%3A%2F%2Fwww.openbioinformatics.org%2Fannovar%2Fdownload%2F0wgxR2rIVP%2Fannovar.latest.tar.gz&data=05%7C02%7C%7C64a5dc472e7c4c79464b08dcc67ab954%7Ca6fa3b030a3c42588433a120dffcd348%7C0%7C0%7C638603478059776283%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C20000%7C%7C%7C&sdata=a2iaWON7xKH2AgILvLFmyvZR1o%2BSX%2FALwP3ZyHPr8yk%3D&reserved=0 annotate.latest.tar.gz
 
#Unpack the packages using the below command line
tar xvfz annovar.latest.tar.gz

#Download appropriate database files using annotate_variation.pl
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240611 humandb/


# Split samples in .vcf file

for sample in $(bcftools query -l "$OUTPUT_DIR/extract_genes.vcf"); do
    bcftools view -c1 -s ${sample} "$OUTPUT_DIR/extract_genes.vcf" -o "$OUTPUT_DIR/${sample}.vcf"
done


#Script for positive controls
#!/bin/bash
#Directory where the VCF files are located  
OUTPUT_DIR="/home/lmadula/fourth_qcoutput"

#Directory where selected VCF files will be copied to
SELECTED_VCF="/home/lmadula/Solved_cases"

#Combined selected file (merge)
COMBINED_FILE="$SELECTED_VCF/combined.vcf"

#Patterns to match specific vcf files
PATTERNS=("_001505_s1.vcf" "_001520_s1.vcf" "_001508_s1.vcf" "_001522_s1.vcf" "_001521_s1.vcf" "_001525_s1.vcf" "_001517_s1.vcf" "_001504_s1.vcf" "_001512_s1.vcf" "_001499_s1.vcf" "_001515_s1.vcf" "_001526_s1.vcf")

#Find all VCF files and filter using grep for the specified patterns
find "$OUTPUT_DIR" -type f -name "*.vcf" | grep -E "$(IFS=\|; echo "${PATTERNS[*]}")" | while read -r file; do
    echo "Found and copying file: $file"
    cp "$file" "$SELECTED_VCF /"
done

#Combine the selected VCF files into one file
echo "Combining files into $COMBINED_FILE"
cat "$SELECTED_VCF"/*.vcf > "$COMBINED_FILE"

#Optionally, you might want to remove the individual VCF files after combining (It's not essential)
rm "$SELECTED_VCF"/*.vcf

echo "Combination complete. Check the '$COMBINED_FILE' file for the combined VCF data."

#Change the combined positive control vcf file into a file.avinput
perl convert2annovar.pl -format vcf4 /path/to/file/file.vcf > /path/to/file/file.avinput
 

#Then do the annotation using annovar
~/annovar/table_annovar.pl ~/selected_vcf/combined.avinput ~/annovar/humandb -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,clinvar -operation g,r,f,f,f,f -nastring NA -csvout -polish


#!/bin/bash

OUTPUT_DIR="$HOME/fourth_qcoutput"  # Use HOME for output directory
SELECTED_VCF="$HOME/Solved_cases"     # Use HOME for selected VCF directory
COMBINED_FILE="$SELECTED_VCF/combined.vcf"
ANNOVAR_DB="$HOME/annovar/humandb"    # Use HOME for ANNOVAR database directory

# Patterns to match specific VCF files
PATTERNS=("_001505_s1.vcf" "_001520_s1.vcf" "_001508_s1.vcf" "_001522_s1.vcf" "_001521_s1.vcf" "_001525_s1.vcf" "_001517_s1.vcf" "_001504_s1.vcf" "_001512_s1.vcf" "_001499_s1.vcf" "_001515_s1.vcf" "_001526_s1.vcf")

# Find and copy selected VCF files
find "$OUTPUT_DIR" -type f -name "*.vcf" | grep -E "$(IFS=\|; echo "${PATTERNS[*]}")" | while read -r file; do
    echo "Found and copying file: $file"
    cp "$file" "$SELECTED_VCF/"
done

# Combine the selected VCF files into one file
echo "Combining files into $COMBINED_FILE"
cat "$SELECTED_VCF"/*.vcf > "$COMBINED_FILE"

# Optionally, remove the individual VCF files after combining
rm "$SELECTED_VCF"/*.vcf

echo "Combination complete. Check the '$COMBINED_FILE' file for the combined VCF data."

# Convert the combined VCF file to .avinput format
perl "$HOME/annovar/convert2annovar.pl" -format vcf4 "$COMBINED_FILE" > "$SELECTED_VCF/combined.avinput"

# Annotate using ANNOVAR
perl "$HOME/annovar/table_annovar.pl" "$SELECTED_VCF/combined.avinput" "$ANNOVAR_DB" -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,clinvar -operation g,r,f,f,f,f -nastring NA -csvout -polish

# Move the final annotation output to the selected VCF directory
mv "myanno.hg38_multianno.csv" "$SELECTED_VCF/myanno.hg38_multianno.csv" 





# Assign the missingness threshold to a variable
MISSINGNESS_THRESHOLD=0.1

# Use the variable in the PLINK command
plink --bfile gp2__geno_0.05 --mind $MISSINGNESS_THRESHOLD --make-bed --out gp2_mind_$MISSINGNESS_THRESHOLD
 


 https://rpubs.com/lumumba99/1026665
 ## Genomic data
 https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.amp-pd.org/data/genomic-data&ved=2ahUKEwin3Kfi8qKJAxXM-LsIHR5rO-wQFnoECCYQAQ&usg=AOvVaw31_rprbJIyb3LtYSuZpc-s
