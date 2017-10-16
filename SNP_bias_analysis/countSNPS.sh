#!/bin/bash

#### INPUTS ####
bam="$1"
igv_dir="$2"
full_mpileup="$3"

#### SCRIPT BEGIN ####
count_script_path='/hpc/users/holtj02/SINAI_SCRIPTS/MISC/countPerBPPercentage.py'
hg19Fasta='/sc/orga/projects/GCF/REFERENCES/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
basename=${bam%.Aligned.out.sorted.bam}
append_file=${bam%.Aligned.out.sorted.bam}.append.sam
sam_file=${bam%.bam}.sam
#### LOAD MODULES ####

module purge
module load samtools
module load python
module load py_packages
module load bcftools

if [[ $full_mpileup == '' ]]; then
    samtools mpileup -f $hg19Fasta $bam > ${bam%.Aligned.out.sorted.bam}.mpileup
    full_mpileup=${bam%.Aligned.out.sorted.bam}.mpileup
fi

if [[ ! -f $append_file ]]; then
    ## extract the header from the provided bam file... this file will be appended reads which have high enough coverage and aren't garbage
    samtools view -H $bam > $append_file
    samtools view $bam > $sam_file
fi

## append append_file with reads with high enough coverage, and good quality pass distribution 
#  also write the typical output to the txt file provided (output_txt)
#  then convert the sam to bed format to provide regions of interest to the mpileup
output_txt=${bam%.Aligned.out.sorted.bam}.new.txt
#sam2bed_file=${append_file%.sam}.bed
mkdir -p sams
sample_name=$(basename `pwd`)
if [[ ! -f sams/$output_txt ]]; then
    # make it so the count_script also outputs the sam alignments to the provided $append_file
    # this also needs to be adjusted to take in a specified output file for the output of this script
    cd sams
    sam_dir=`pwd`
    python $count_script_path ../$sam_file ../$full_mpileup 10 $output_txt ../$append_file
    sams=($(ls sam_alignments_*.sam))
    echo "there were ${#sams[@]} sam_alignments_*.sam file generated by $count_script_path"
    mkdir bams beds first_vcf mpileup vcf snps
    for i in ${sams[@]}; do
        sam2bed_file=beds/${i%.sam}.bed
        append_file_bam=bams/${i%.sam}.bam
        append_file_first_vcf=first_vcf/${i%.sam}.1st.vcf
        append_file_mpileup_raw=mpileup/${i%.sam}.mpileup
        append_file_vcf=vcf/${i%.sam}.vcf
        snps_file_vcf=snps/"snps.${i%.sam}.vcf"
        sam2bed $i $sam2bed_file
        samtools view -b $i > $append_file_bam
        cd bams
        samtools index $append_file_bam
        cd $igv_dir
        ln -s $sam_dir/$append_file_bam $(basename $append_file_bam)
        ln -s $sam_dir/$append_file_bam.bai $(basename $append_file_bam).bai
        cd $sam_dir
        samtools mpileup -L $sam2bed_file -C 50 -Q 20 -uf $hg19Fasta $append_file_bam > $append_file_first_vcf
        samtools mpileup -L $sam2bed_file -C 50 -Q 20 -f $hg19Fasta $append_file_bam > $append_file_mpileup_raw
        bcftools call -m $append_file_first_vcf > $append_file_vcf
        bcftools view --no-header -v snps $append_file_vcf | vcfutils.pl varFilter -d 25 > $snps_file_vcf
    done
fi

exit


