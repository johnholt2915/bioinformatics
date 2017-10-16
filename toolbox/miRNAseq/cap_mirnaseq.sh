#!/bin/bash
qcdir="$1"
sampleName="$2"

cd $qcdir/Raw/RNA*/Sample_$sampleName
sampleDir=`pwd`
# find the fastq file for this sample.
fqs=($(ls *.fastq.gz))
mkdir -p cap_mirna
cd cap_mirna
outdir=`pwd`
cutadapt_out=$outdir/$sampleName.cutadapt.fastq.gz
if [[ ! -f $cutadapt_out ]]; then
    ca_params="-b AGATCGGAAGAGCACACGTCT -O 3 -m 17 -f fastq -o $cutadapt_out"
    module purge
    module load python
    module load py_packages
    cutadapt $ca_params ${fqs[@]}
fi
bowtie_out=$outdir/$sampleName'_Aligned.sortedByCoord.bam.bai'
if [[ ! -f $bowtie_out ]]; then
    bowtie_genome="/sc/orga/projects/GCF/REFERENCES/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome"
    bt_params="-p 12 -S --sam-RG ID:$sampleName\tLB:$sampleName\tPL:Illumina\tPU:MiSeq\tSM:$sampleName\tCN:MSSM\tDS:mirnaseq -q -n 1 -e 80 -l 30 -a -m 5 --best --strata"
    module purge
    module load bowtie
    module load samtools
    if [[ ${#fqs[@]} -gt 1 ]]; then
        bowtie_input_files="${fqs[@]}"
    else
        bowtie_input_files="-1 ${fqs[0]} -2 ${fqs[1]}"
    fi
    bowtie $bt_params $bowtie_genome $bowtie_input_files $outdir/$sampleName'_Aligned.out.sam'
    samtools view -b $outdir/$sampleName'_Aligned.out.sam'
    samtools sort -@ 12 -T tmp -O bam -o $outdir/$sampleName'_Aligned.sortedByCoord.bam' $outdir/$sampleName'_Aligned.out.bam'
    samtools index $outdir/$sampleName'_Aligned.sortedByCoord.bam'
fi


Usage: 
bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]

  <m1>    Comma-separated list of files containing upstream mates (or the
          sequences themselves, if -c is set) paired with mates in <m2>
  <m2>    Comma-separated list of files containing downstream mates (or the
          sequences themselves if -c is set) paired with mates in <m1>
  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
          a mixture of paired and unpaired.  Specify "-" for stdin.
  <s>     Comma-separated list of files containing unpaired reads, or the
          sequences themselves, if -c is set.  Specify "-" for stdin.
  <hit>   File to write hits to (default: stdout)
Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
  -f                 query input files are (multi-)FASTA .fa/.mfa
  -r                 query input files are raw one-sequence-per-line
  -c                 query sequences given on cmd line (as <mates>, <singles>)
  -C                 reads and index are in colorspace
  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
  -s/--skip <int>    skip the first <int> reads/pairs in the input
  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
  --phred33-quals    input quals are Phred+33 (default)
  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
  --solexa-quals     input quals are from GA Pipeline ver. < 1.3
  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
  --integer-quals    qualities are given as space-separated integers (not ASCII)
  --large-index      force usage of a 'large' index, even if a small one is present
Alignment:
  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
    or
  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
  -l/--seedlen <int> seed length for -n (default: 28)
  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)
  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)
  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
  --nofw/--norc      do not align to forward/reverse-complement reference strand
  --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)
  -y/--tryhard       try hard to find valid alignments, at the expense of speed
  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
Reporting:
  -k <int>           report up to <int> good alignments per read (default: 1)
  -a/--all           report all alignments per read (much slower than low -k)
  -m <int>           suppress all alignments if > <int> exist (def: no limit)
  -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
  --best             hits guaranteed best stratum; ties broken by quality
  --strata           hits in sub-optimal strata aren't reported (requires --best)
Output:
  -t/--time          print wall-clock time taken by search phases
  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
  --quiet            print nothing but the alignments
  --refout           write alignments to files refXXXXX.map, 1 map per reference
  --refidx           refer to ref. seqs by 0-based index rather than name
  --al <fname>       write aligned reads/pairs to file(s) <fname>
  --un <fname>       write unaligned reads/pairs to file(s) <fname>
  --max <fname>      write reads/pairs over -m limit to file(s) <fname>
  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output
  --fullref          write entire ref name (default: only up to 1st space)
Colorspace:
  --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)
     or
  --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
  --col-cseq         print aligned colorspace seqs as colors, not decoded bases
  --col-cqual        print original colorspace quals, not decoded quals
  --col-keepends     keep nucleotides at extreme ends of decoded alignment
SAM:
  -S/--sam           write hits in SAM format
  --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments
  --sam-nohead       supppress header lines (starting with @) for SAM output
  --sam-nosq         supppress @SQ header lines for SAM output
  --sam-RG <text>    add <text> (usually "lab=value") to @RG line of SAM header
Performance:
  -o/--offrate <int> override offrate of index; must be >= index's offrate
  -p/--threads <int> number of alignment threads to launch (default: 1)
  --mm               use memory-mapped I/O for index; many 'bowtie's can share
  --shmem            use shared mem for index; many 'bowtie's can share
Other:
  --seed <int>       seed for random number generator
  --verbose          verbose output (for debugging)
  --version          print version information and quit
  -h/--help          print this usage message
