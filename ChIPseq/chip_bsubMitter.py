import os
import sys
import subprocess
import optparse

def runCommand(cmd):
    subprocess.call(cmd,shell=True)

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-w', '--workdir', dest='workdir')
    parser.add_option( '-c', '--chip-fq', dest='chip_fq')
    parser.add_option( '-i', '--input-fq', dest='input_fq')
    parser.add_option( '-t', '--threads', dest='threads')
    parser.add_option( '-r', '--read-len', dest='read_len')
    ( options, args ) = parser.parse_args()

    if not options.chip_fq:
        print "Please provide a ChIPseq fastq file to run this analysis."
        sys.exit(1)
    if not options.input_fq:
        options.input_fq = ''
    if not options.threads:
        options.threads = '4'
    if not options.read_len:
        options.read_len = '100'
    return options

def bsubMit(workdir,chip_fq,input_fq,nproc,read_len):
    chip_name = chip_fq.split('/')[-1].split('.')[0]
    if input_fq == '':
        input_name = ''
    else:
        input_name = input_fq.split('/')[-1].split('.')[0]
    rg_info = "@RG\tID:$chip_name\tPL:Illumina\tLB:$chip_name\tID:$chip_name"
    outfilename = os.path.join(workdir,chip_name+".lsf")
    fp=open(outfilename,"w")
    fp.write("#!/bin/bash\n")
    fp.write("#BSUB -J "+name+"_chip_align_peak_qc\n")
    fp.write("#BSUB -W 3:00\n")
    fp.write("#BSUB -m manda\n")
    fp.write("#BSUB -q premium\n")
    fp.write("#BSUB -P acc_GCF\n")
    fp.write("#BSUB -o "+workdir+"/"+name+"_bsub.out\n")
    fp.write("#BSUB -o "+workdir+"/"+name+"_bsub.err\n")
    fp.write("#BSUB -n "+nproc+"\n")
    fp.write("#BSUB -u john.holt@mssm.edu\n")
    fp.write("#BSUB -L /bin/bash\n")
    fp.write("#BSUB -R rusage[mem=4200]\n")
    fp.write("#BSUB -R span[hosts=1]\n")
    fp.write("workdir="+workdir+"\n")
    fp.write("chip_fq="+chip_fq+"\n")
    fp.write("chip_name="+chip_name+"\n")
    fp.write("input_fq="+input_fq+"\n")
    fp.write("input_name="+input_name+"\n")
    fp.write('threads='+nproc+'\n')
    fp.write('chip_qc_script=/hpc/users/holtj02/SINAI_SCRIPTS/NGS_Analysis_Templates/chip_sample_qc.R\n')
    fp.write('bwa_genome=/sc/orga/projects/GCF/REFERENCES/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa\n') 
    fp.write('rg_info="@RG\tID:$chip_name\tPL:Illumina\tLB:$chip_name\tID:$chip_name"\n')
    fp.write('module load bwa samtools macs/2.1.0 R\n')
    fp.write('mkdir -p $workdir/data/$chip_name/analysis\n')
    fh.write('mkdir -p $workdir/RData\n')
    fp.write("cd $workdir/data/$chip_name/analysis\n")
    fp.write('bwa_chip_out_file="$workdir/data/$chip_name/analysis/$chip_name.sorted.bam.bai"\n')
    fp.write('if [[ ! -f $bwa_chip_out_file ]]; then\n')
    fp.write('echo "aligning $chip_fq to the human genome"\n')
    fp.write('rg_info="@RG\tID:$chip_name\tPL:Illumina\tLB:$chip_name\tID:$chip_name"\n')
    fp.write('bwa mem -t $threads -R $rg_info $bwa_genome $chip_fq > $chip_name.sam\n')
    fp.write('samtools view -b $chip_name.sam > $chip_name.bam\n')
    fp.write('samtools sort -@ $threads -T tmp -O bam -o $chip_name.sorted.bam $chip_name.bam\n')
    fp.write('samtools index $chip_name.sorted.bam\n')
    fp.write('echo "filtering duplicated reads from $chip_name.sorted.bam"\n')
    fp.write('macs2 filterdup -i $chip_name.sorted.bam --keep-dup=1 -o $chip_name.filterdup.bed 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('post_filter_chip_reads=$(wc -l $chip_name.filterdup.bed | cut -d ' ' -f 1)\n')
    fp.write('fi\n')
    fp.write('bwa_input_out_file="$workdir/data/$chip_name/analysis/$input_name.sorted.bam.bai"\n')
    fp.write('if [[ ! -f $bwa_input_out_file ]]; then\n')
    fp.write('echo "aligning $input_fq to human genome"\n')
    fp.write('rg_info="@RG\tID:$input_name\tPL:Illumina\tLB:$input_name\tID:$input_name"\n')
    fp.write('bwa mem -t $threads -R $rg_info $bwa_genome $input_fq > $input_name.sam\n')
    fp.write('samtools view -b $input_name.sam > $input_name.bam\n')
    fp.write('samtools sort -@ $threads -T tmp -O bam -o $input_name.sorted.bam $input_name.bam\n')
    fp.write('samtools index $input_name.sorted.bam\n')
    fp.write('echo "filtering duplicated reads from $input_name.sorted.bam"\n')
    fp.write('macs2 filterdup -i $input_name.sorted.bam --keep-dup=1 -o $input_name.filterdup.bed 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('post_filter_input_reads=$(wc -l $input_name.filterdup.bed | cut -d ' ' -f 1)\n')
    fp.write('fi\n')
    fp.write("## peak calling ##\n")
    fp.write('echo "predicting footprint size of DNA bound proteins."\n')
    fp.write('macs2 predictd -i $chip_name.filterdup.bed -g hs -m 5 50 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('mv predictd predictd.$chip_name.R\n')
    fp.write('extsize=') # figure out how to collect this info...
    fp.write('macs2 pileup -i $chip_name.filterdup.bed -o $chip_name.filterdup.pileup.bdg --extsize $extsize 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('echo "normalizing ChIP sample to Input based on footprint size, input sequencing depth, and multiple local bias fields."\n')
    fp.write('macs2 pileup -i $input_name.filterdup.bed -B --extsize $(($extsize/2)) -o d_bg_$chip_name.norm.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 pileup -i $input_name.filterdup.bed -B --extsize 500 -o 1k_bg_$input_name.norm.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 pileup -i $input_name.filterdup.bed -B --extsize 5000 -o 10k_bg_$input_name.norm.bdg  1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgopt -i 1k_bg_$input_name.norm.bdg -m multiply -p $(bc -l <<< "$extsize/1000") -o 1k_bg_$chip_name.norm.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgopt -i 10k_bg_$input_name.norm.bdg -m multiply -p $(bc -l <<< "$extsize/1000") -o 10k_bg_$chip_name.norm.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgcmp -m max -t 1k_bg_$chip_name.norm.bdg -c 10k_bg_$chip_name.norm.bdg -o 1k_10k_bg_$chip_name.norm.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgcmp -m max -t 1k_10k_bg_$chip_name.norm.bdg -c d_bg_$chip_name.norm.bdg -o d_1k_10k_bg_$chip_name.norm.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgopt -i d_1k_10k_bg_$chip_name.norm.bdg -m max -p $(bc -l <<< "$(($post_filter_input_reads * $extsize))/2700000000") -o local_bias_raw_$chip_name.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgopt -i local_bias_raw_$chip_name.bdg -m multiply -p $(bc -l <<< "$post_filter_chip_reads/$post_filter_input_reads") -o local_lambda_$chip_name.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgcmp -t $chip_name.filterdup.pileup.bdg -c local_lambda_$chip_name.bdg -m qpois -o $chip_name.qvalue.bdg 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('macs2 bdgpeakcall -i $chip_name.qvalue.bdg -c 1.301 -l $extsize -g '+read_len+' -o $chip_name.peaks.bed 1>> Log.$chip_name.out 2>> Log.$chip_name.err\n')
    fp.write('Rscript $chip_qc_script $workdir/data/$chip_name/$chip_name.sorted.bam $workdir/data/$chip_name/$chip_name.peaks.bed $chip_name\n')
    fp.close()
    runCommand("bsub < "+outfilename)

def main():
    options = argParser()
    chip_fq = options.chip_fq
    input_fq = options.input_fq
    workdir = options.workdir
    threads = options.threads
    read_len = options.read_len
    bsubMit(workdir,chip_fq,input_fq,threads,read_len)

if __name__ == '__main__':
    main()

