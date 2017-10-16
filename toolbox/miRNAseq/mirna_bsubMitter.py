import os
import sys
import subprocess
import optparse

def runCommand(cmd):
    subprocess.call(cmd,shell=True)

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-w', '--workdir', dest='workdir')
    parser.add_option( '-1', '--mirna-fq', dest='mirna_fq')
    parser.add_option( '-2', '--mirna-fq-read-2', dest='read_2')
    parser.add_option( '-a', '--adapter-fasta', dest='adapter_fa')
    parser.add_option( '-i', '--bowtie2-index', dest='bwt2_idx')
    parser.add_option( '-t', '--threads', dest='threads')
    
    ( options, args ) = parser.parse_args()
    if not options.mirna_fq:
        print "Please provide a miRNAseq fastq file to run this analysis."
        sys.exit()
    if not options.bwt2_idx:
        options.bwt2_idx = "/sc/orga/projects/GCF/miRNA_QC/mirBase/Bowtie2Index/homo_sapiens_mature_mirna_cdna"
    if not options.threads:
        options.threads = '20'
    return options

def bsubMit(workdir,mirna_fq,mirna_fq_read_2,adapter_fa,bwt2_idx,nproc):
    filename= mirna_fq.split('/')[-1]
    mirna_name = filename.split('.')[0]
    # infer the file ending from the fastq
    if '.fq.gz' in filename:
        fq_tail = '.fq.gz'
    elif '.fastq.gz' in filename:
        fq_tail = '.fastq.gz'
    elif '.fq' in filename:
        fq_tail = '.fq'
    elif '.fastq' in filename:
        fq_tail = '.fastq'
    # if a second fastq if provided, it is assumed to be the paired end file complementary to mirna_fq
    # two different trimmomatic and bowtie2 commands are required in this case
    if not mirna_fq_read_2:
        trimm_cmd = "java -jar $TRIMMOMATIC_JAR SE -threads $threads -phred33 -trimlog $workdir/data/$mirna_name/analysis/trimmed_fastq/Log.trimmomatic.txt $mirna_fq $trimm_output ILLUMINACLIP:$adapter_fa:2:30:10 MINLEN:15"
        bwt2_cmd = "bowtie2 -k 100 -p %s --local --very-sensitive-local $rg_info --un-gz $workdir/data/$mirna_name/analysis/Unaligned.$mirna_name.$fq_tail -x $bwt2_idx -U $trimm_output -S $mirna_name.sam 1> Log.alignment_stats.out 2> Log.alignment_stats.err" % (nproc)
    else:
        trimm_cmd = "java -jar $TRIMMOMATIC_JAR PE -threads $threads -phred33 -trimlog $workdir/data/$mirna_name/analysis/trimmed_fastq/Log.trimmomatic.txt $mirna_fq $read_2_fq $trimm_output_unpaired_1 $trimm_output_paired_1 $trimm_output_unpaired_2 $trimm_output_paired_2 ILLUMINACLIP:$adapter_fa:2:30:10 MINLEN:15"
        bwt2_cmd = "bowtie2 -k 100 -p %s --local --very-sensitive-local $rg_info --un-gz $workdir/data/$mirna_name/analysis/Unaligned.$mirna_name.$fq_tail -x $bwt2_idx -1 $trimm_output_paired_1 -2 $trimm_output_paired_2 -S $mirna_name.sam 1> Log.alignment_stats.out 2> Log.alignment_stats.err" % (nproc)
    rg_info = "ID:$mirna_name\tPL:Illumina\tLB:$mirna_name\tID:$mirna_name"
    outfilename = os.path.join(workdir,mirna_name+".lsf")
    fp=open(outfilename,"w")
    fp.write("#!/bin/bash\n")
    fp.write("#BSUB -J "+mirna_name+"_mirna_align_qc\n")
    fp.write("#BSUB -W 3:00\n")
    fp.write("#BSUB -m manda\n")
    fp.write("#BSUB -q premium\n")
    fp.write("#BSUB -P acc_GCF\n")
    fp.write("#BSUB -o "+workdir+"/"+mirna_name+"_bsub.out\n")
    fp.write("#BSUB -o "+workdir+"/"+mirna_name+"_bsub.err\n")
    fp.write("#BSUB -n "+nproc+"\n")
    fp.write("#BSUB -u john.holt@mssm.edu\n")
    fp.write("#BSUB -L /bin/bash\n")
    fp.write("#BSUB -R rusage[mem=4200]\n")
    fp.write("#BSUB -R span[hosts=1]\n")
    fp.write("workdir="+workdir+"\n")
    fp.write("mirna_fq="+mirna_fq+"\n")
    fp.write("mirna_name="+mirna_name+"\n")
    if mirna_fq_read_2:
        fp.write('read_2_fq='+mirna_fq_read_2+'\n')
    fp.write('fq_tail='+fq_tail+'\n')
    fp.write('threads='+nproc+'\n')
    fp.write('adapter_fa='+adapter_fa+'\n')
    fp.write('mirna_qc_script=/hpc/users/holtj02/SINAI_SCRIPTS/MISC/bioinformatics_git_repos/miRNAseq/miRNA_qc.py\n')
    fp.write('bwt2_idx='+bwt2_idx+'\n')
    #fp.write('rg_info="@RG\tID:$mirna_name\tPL:Illumina\tLB:$mirna_name\tID:$mirna_name"\n')
    fp.write('module load bowtie2 samtools trimmomatic\n')
    fp.write('mkdir -p $workdir/data/$mirna_name/analysis/trimmed_fastq\n')
    #fh.write('mkdir -p $workdir/RData\n')
    fp.write("cd $workdir/data/$mirna_name/analysis\n")
    if mirna_fq_read_2:
        fp.write('trimm_output_unpaired_1="$workdir/data/$mirna_name/analysis/trimmed_fastq/$mirna_name.trimmed.unpaired.read1.$fq_tail"\n')
        fp.write('trimm_output_unpaired_2="$workdir/data/$mirna_name/analysis/trimmed_fastq/$mirna_name.trimmed.unpaired.read2.$fq_tail"\n')
        fp.write('trimm_output_paired_1="$workdir/data/$mirna_name/analysis/trimmed_fastq/$mirna_name.trimmed.paired.read1.$fq_tail"\n')
        fp.write('trimm_output_paired_2="$workdir/data/$mirna_name/analysis/trimmed_fastq/$mirna_name.trimmed.paired.read2.$fq_tail"\n')
        fp.write('if [[ ! -f $trimm_output_paired_1 ]]; then\n')
    else:
        fp.write('trimm_output="$workdir/data/$mirna_name/analysis/trimmed_fastq/$mirna_name.trimmed.$fq_tail"\n')
        fp.write('if [[ ! -f $trimm_output ]]; then\n')
    fp.write(trimm_cmd+'\n')
    fp.write('fi\n')
    fp.write('bwt2_mirna_out_file="$workdir/data/$mirna_name/analysis/$mirna_name.sorted.bam.bai"\n')
    fp.write('if [[ ! -f $bwt_mirna_out_file ]]; then\n')
    fp.write('echo "aligning $mirna_fq to the human genome"\n')
    fp.write('rg_info="--rg-id ID --rg $mirna_name --rg-id PL --rg Illumina --rg-id LB --rg $mirna_name --rg-id SN --rg $mirna_name"\n')
    fp.write(bwt2_cmd+'\n')
    fp.write('samtools view -H $mirna_name"_Aligned.out.sam" > $mirna_name"_Aligned.out.sam"\n')Aligned.sortedByCoord.out.bam
    fp.write('samtools view -F 4 $mirna_name.sam >> Aligned.$mirna_name.sam\n')
    fp.write('samtools view -b $mirna_name.sam > $mirna_name.bam\n')
    fp.write('samtools sort -@ $threads -T tmp -O bam -o $mirna_name.sorted.bam $mirna_name.bam\n')
    fp.write('samtools index $mirna_name.sorted.bam\n')
    fp.write('fi\n')
    fp.write('mirna_qc_output="$workdir/data/$mirna_name/analysis/metrics.tsv"\n')
    fp.write('if [[ ! -f $mirna_qc_output ]]; then\n')
    fp.write('python $mirna_qc_script -s $workdir/data/$mirna_name/analysis/$mirna_name.sam -o $workdir/data/$mirna_name/analysis/$mirna_name.counts.txt -a $workdir/data/$mirna_name/analysis/Log.alignment_stats.out\n')
    fp.write('fi\n')
    fp.write('exit\n')
    fp.close()
    #runCommand("bsub < "+outfilename)

def main():
    options = argParser()
    workdir = options.workdir
    mirna_fq = options.mirna_fq
    mirna_fq_read_2 = options.read_2
    adapter_fa = options.adapter_fa
    bwt2_idx = options.bwt2_idx
    threads = options.threads
    bsubMit(workdir,mirna_fq,mirna_fq_read_2,adapter_fa,bwt2_idx,threads)

if __name__ == '__main__':
    main()

