import sys
import os
import subprocess

def runCommand(command):
    subprocess.call(command, shell=True)

def runBsub(workdir,fastq,min_len,mismatches):
    name = fastq.split('.fa')[0]
    lsf_file = name+".lsf"
    outfile = open(lsf_file,'w')
    outfile.write("#!/bin/bash\n");
    outfile.write("#BSUB -J "+name+"\n");
    outfile.write("#BSUB -W 04:00\n");
    outfile.write("#BSUB -n 8\n");
    outfile.write("#BSUB -m manda\n");
    outfile.write("#BSUB -q premium\n");
    outfile.write("#BSUB -P acc_GCF\n");
    outfile.write("#BSUB -o "+name+"_bsub.out\n");
    outfile.write("#BSUB -e "+name+"_bsub.err\n");
    outfile.write("#BSUB -u john.holt@mssm.edu\n");
    outfile.write("#BSUB -L /bin/bash\n");
    outfile.write("#BSUB -R rusage[mem=4200]\n");
    outfile.write("#BSUB -R span[hosts=1]\n");
    outfile.write('workdir="'+workdir+'"\n')
    outfile.write('fqstart="'+fastq+'"\n')
    outfile.write('min_len="'+min_len+'"\n')
    outfile.write('mismatches="'+mismatches+'"\n')
    outfile.write('filter_script="/hpc/users/holtj02/SINAI_SCRIPTS/NGS_Analysis_Templates/otherFilter.py"\n')
    outfile.write('bwa_genome="/sc/orga/projects/GCF/REFERENCES/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"\n')
    outfile.write('count_BP_Percentage="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/countPerBPPercentage.py"\n')
    outfile.write('count_SNPS="/hpc/users/holtj02/SINAI_SCRIPTS/MISC/countSNPS.sh"\n')
    outfile.write('igv_dir="/hpc/users/holtj02/www/bams/$(basename $workdir)"\n')
    outfile.write('module purge\n')
    outfile.write('module load python py_packages\n')
    outfile.write('module load bwa\n')
    outfile.write('module load samtools\n')
    outfile.write('cd $workdir/fastq\n')
    outfile.write('echo "fqstart: $fastq"\n')
    outfile.write('if [[ $fqstart == *".gz" ]]; then\n')
    outfile.write('base=$(basename ${fqstart%.fastq.gz})\n')
    outfile.write('fq=${fqstart%.fastq.gz}.fastq\n')
    outfile.write('echo "base: $base"\n')
    outfile.write('echo "mkdir -p $base"\n')
    outfile.write('mkdir -p $base\n')
    outfile.write('cd $base\n')
    outfile.write('echo "cd $base"\n')
    outfile.write('if [[ ! -f $base.done.txt ]]; then\n')
    outfile.write('gunzip $fqstart\n')
    outfile.write('python $filter_script $fq $min_len $mismatches 1> $base.Log.out 2> $base.Log.err\n')
    outfile.write('echo $(du -h ./*) > $base.done.txt\n')
    outfile.write('rm fail_filter*.fastq\n')
    outfile.write('gzip ${fqstart%.gz}\n')
    outfile.write('echo "gzipped ${fqstart%.gz}"\n')
    outfile.write('else\n')
    outfile.write('echo "$base.done.txt exists... filtering already complete..."\n')
    outfile.write('fi\n')
    outfile.write('else\n')
    outfile.write('base=$(basename ${fqstart%.fastq})\n')
    outfile.write('fq=$fqstart\n')
    outfile.write('echo "base: $base"\n')
    outfile.write('echo "mkdir -p $base"\n')
    outfile.write('mkdir -p $base\n')
    outfile.write('cd $base\n')
    outfile.write('echo "cd $base"\n')
    outfile.write('if [[ ! -f $base.done.txt ]]; then\n')
    outfile.write('python $filter_script $fq $min_len $mismatches 1> $base.Log.out 2> $base.Log.err\n')
    outfile.write('echo $(du -h ./*) > $base.done.txt\n')
    outfile.write('rm fail_filter*.fastq\n')
    outfile.write('gzip $fqstart\n')
    outfile.write('echo "gzipped $fqstart"\n')
    outfile.write('else\n')
    outfile.write('echo "$base.done.txt exists... filtering already complete..."\n')
    outfile.write('gzip $fqstart\n')
    outfile.write('echo "gzipped $fqstart"\n')
    outfile.write('fi\n')
    outfile.write('fi\n')
    outfile.write('pf_fastqs=$(ls pass_filter*)\n')
    outfile.write('mkdir -p $workdir/bwa_out/$base\n')
    outfile.write('for pf_fq in $pf_fastqs; do\n')
    outfile.write('if [[ $pf_fq == *"5prime"* ]]; then\n')
    outfile.write('mkdir -p $workdir/bwa_out/$base/5prime_$base\n')
    outfile.write('cd $workdir/bwa_out/$base/5prime_$base\n')
    outfile.write('rg="@RG\\tID:5prime_$base\\tLB:5prime_$base\\tPL:Illumina\\tPU:HiSeq4000\\tSM:5prime_$base\\tCN:MSSM\\tDS:5prime_$base"\n')
    outfile.write('echo "bwa mem -t 8 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 5prime_$base.Aligned.out.sam"\n')
    outfile.write('bwa mem -t 8 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 5prime_$base.Aligned.out.sam\n')
    outfile.write('echo "samtools view -b 5prime_$base.Aligned.out.sam > 5prime_$base.Aligned.out.bam"\n')
    outfile.write('samtools view -b 5prime_$base.Aligned.out.sam > 5prime_$base.Aligned.out.bam\n')
    outfile.write('echo "samtools sort -@ 4 -T tmp -O bam -o 5prime_$base.Aligned.out.sorted.bam 5prime_$base.Aligned.out.bam"\n')
    outfile.write('samtools sort -@ 4 -T tmp -O bam -o 5prime_$base.Aligned.out.sorted.bam 5prime_$base.Aligned.out.bam\n')
    outfile.write('echo "samtools index 5prime_$base.Aligned.out.sorted.bam"\n')
    outfile.write('samtools index 5prime_$base.Aligned.out.sorted.bam\n')
    outfile.write('echo "bash $count_SNPS 5prime_$base.Aligned.out.sorted.bam"\n')
    outfile.write('bash $count_SNPS 5prime_$base.Aligned.out.sorted.bam $igv_dir\n')
    outfile.write('elif [[ $pf_fq == *"3prime"* ]]; then\n')
    outfile.write('mkdir -p $workdir/bwa_out/$base/3prime_$base\n')
    outfile.write('cd $workdir/bwa_out/$base/3prime_$base\n')
    outfile.write('rg="@RG\\tID:3prime_$base\\tLB:3prime_$base\\tPL:Illumina\\tPU:HiSeq4000\\tSM:3prime_$base\\tCN:MSSM\\tDS:3prime_$base"\n')
    outfile.write('echo "bwa mem -t 8 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 3prime_$base.Aligned.out.sam"\n')
    outfile.write('bwa mem -t 8 -R $rg $bwa_genome $workdir/fastq/$base/$pf_fq > 3prime_$base.Aligned.out.sam\n')
    outfile.write('echo "samtools view -b 3prime_$base.Aligned.out.sam > 3prime_$base.Aligned.out.bam"\n')
    outfile.write('samtools view -b 3prime_$base.Aligned.out.sam > 3prime_$base.Aligned.out.bam\n')
    outfile.write('echo "samtools sort -@ 4 -T tmp -O bam -o 3prime_$base.Aligned.out.sorted.bam 3prime_$base.Aligned.out.bam"\n')
    outfile.write('samtools sort -@ 4 -T tmp -O bam -o 3prime_$base.Aligned.out.sorted.bam 3prime_$base.Aligned.out.bam\n')
    outfile.write('echo "samtools index 3prime_$base.Aligned.out.sorted.bam"\n')
    outfile.write('samtools index 3prime_$base.Aligned.out.sorted.bam\n')
    outfile.write('echo "bash $count_SNPS 3prime_$base.Aligned.out.sorted.bam"\n')
    outfile.write('bash $count_SNPS 3prime_$base.Aligned.out.sorted.bam $igv_dir\n')
    outfile.write('fi\n')
    outfile.write('done\n')
    outfile.write('exit\n')
    outfile.write('\n')
    outfile.close()
    runCommand("bsub < "+lsf_file)

def main():
    workdir = sys.argv[1]
    fastq = sys.argv[2]
    min_len = sys.argv[3]
    mismatches = sys.argv[4]
    runBsub(workdir,fastq,min_len,mismatches)

if __name__ == '__main__':
    main()

