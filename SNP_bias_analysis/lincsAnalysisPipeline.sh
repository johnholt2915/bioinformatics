#!/bin/bash

workdir="$1"
min_len="$2"
mismatches="$3"

filter_script="/hpc/users/holtj02/SINAI_SCRIPTS/NGS_Analysis_Templates/otherFilter.py"
bwa_genome="/sc/orga/projects/GCF/REFERENCES/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
count_BP_Percentage="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/countPerBPPercentage.py"
count_SNPS="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/countSNPS.sh"
bsubMitter="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/bsubMitter.py"
merge_new_txt_script="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/mergeNewTxtFiles.py"
igv_dir="/hpc/users/holtj02/www/bams/$(basename $workdir)"
echo $igv_dir
echo "making $workdir/bwa_out directory"
mkdir -p $workdir/bwa_out
echo "making $igv_dir directory"
mkdir -p $igv_dir

module purge
module load python py_packages
#module load bwa
#module load samtools

cd $workdir/fastq

fqs=($(find $workdir/fastq -maxdepth 1 -name '*.fastq*'))
if [[ ${#fqs} == 0 ]]; then
    echo "need to put fastq files into $workdir/fastq"
fi

echo "${fqs[@]}"

for fqstart in ${fqs[@]}; do
    python $bsubMitter $workdir $fqstart $min_len $mismatches
done

len=$(bjobs | wc -l)
count=0
while [[ $len -gt 0 ]]; do 
    count=$(($count+3))
    sleep 3m
    echo "Waited $count minutes for scripts to finish..."
    len=$(bjobs | wc -l)
done
echo "All bsub scripts complete!"

cd $workdir
echo "Merging results and generating urls for igv."
python $merge_new_txt_script $workdir
echo "$workdir/analysis/igv_locations.txt"

exit

# this for loop needs new stuff to iterate over ... 
for i in $(ls); do
    echo "fqstart: $fqstart"
    base=$(basename ${fqstart%.fastq})
    echo "base: $base"
    fq=$(basename $fqstart)
    if [[ $fqstart == *".fastq.gz" ]]; then
        gunzip $fqstart
        base=$(basename ${fqstart%.fastq.gz})
        fq=$base.fastq
    fi
    echo "fq: $fq"
    echo "mkdir -p $base"
    mkdir -p $base
    if [[ ! -f $base/$fq ]]; then
        echo "cp $fq $base/."
        cp $fq $base/.
    else
        echo "$base/$fq exists, do not need to copy $fq to $base/."
    fi
    cd $base
    echo "cd $base"
    if [[ ! -f $base.done.txt ]]; then
        python $filter_script $fq $min_len $mismatches 1> $base.Log.out 2> $base.Log.err
        touch $base.done.txt
    else
        echo "$base.done.txt exists... filtering already complete..."
    fi
    continue
    pf_fastqs=$(ls pass_filter*)
    mkdir -p $workdir/bwa_out/$base
    for pf_fq in $pf_fastqs; do
        if [[ $pf_fq == *"5prime"* ]]; then
            mkdir -p $workdir/bwa_out/$base/5prime_$base
            cd $workdir/bwa_out/$base/5prime_$base
            rg="@RG\tID:5prime_$base\tLB:5prime_$base\tPL:Illumina\tPU:HiSeq4000\tSM:5prime_$base\tCN:MSSM\tDS:5prime_$base"
            echo "bwa mem -t 12 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 5prime_$base.Aligned.out.sam"
            bwa mem -t 12 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 5prime_$base.Aligned.out.sam
            echo "samtools view -b 5prime_$base.Aligned.out.sam > 5prime_$base.Aligned.out.bam"
            samtools view -b 5prime_$base.Aligned.out.sam > 5prime_$base.Aligned.out.bam
            echo "samtools sort -@ 4 -T tmp -O bam -o 5prime_$base.Aligned.out.sorted.bam 5prime_$base.Aligned.out.bam"
            samtools sort -@ 4 -T tmp -O bam -o 5prime_$base.Aligned.out.sorted.bam 5prime_$base.Aligned.out.bam
            echo "samtools index 5prime_$base.Aligned.out.sorted.bam"
            samtools index 5prime_$base.Aligned.out.sorted.bam
            echo "bash $count_SNPS 5prime_$base.Aligned.out.sorted.bam"
            bash $count_SNPS 5prime_$base.Aligned.out.sorted.bam $igv_dir
        elif [[ $pf_fq == *"3prime"* ]]; then
            mkdir -p $workdir/bwa_out/$base/3prime_$base
            cd $workdir/bwa_out/$base/3prime_$base
            rg="@RG\tID:3prime_$base\tLB:3prime_$base\tPL:Illumina\tPU:HiSeq4000\tSM:3prime_$base\tCN:MSSM\tDS:3prime_$base"
            echo "bwa mem -t 12 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 3prime_$base.Aligned.out.sam"
            bwa mem -t 12 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 3prime_$base.Aligned.out.sam
            echo "samtools view -b 3prime_$base.Aligned.out.sam > 3prime_$base.Aligned.out.bam"
            samtools view -b 3prime_$base.Aligned.out.sam > 3prime_$base.Aligned.out.bam
            echo "samtools sort -@ 4 -T tmp -O bam -o 3prime_$base.Aligned.out.sorted.bam 3prime_$base.Aligned.out.bam"
            samtools sort -@ 4 -T tmp -O bam -o 3prime_$base.Aligned.out.sorted.bam 3prime_$base.Aligned.out.bam
            echo "samtools index 3prime_$base.Aligned.out.sorted.bam"
            samtools index 3prime_$base.Aligned.out.sorted.bam
            echo "bash $count_SNPS 3prime_$base.Aligned.out.sorted.bam"
            bash $count_SNPS 3prime_$base.Aligned.out.sorted.bam $igv_dir
        fi
    done
    cd $workdir/fastq
done
echo "remove this exit and rerun so that merging can happen."
# find out how many regions fall within 3 prime utr regions when comparing the *new.txt.
new_files=($(find $workdir/bwa_out -name *new.txt))
cd $workdir
# make a file containing the header from *new.txt files
merge_new_txt_script="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/mergeNewTxtFiles.py"
python $merge_new_txt_script $workdir
merge_new_txt_path="$workdir/analysis/merged_$(basename $workdir).new.txt"

exit

