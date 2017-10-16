import os
import sys
import optparse

"""
Name: John Holt
Date Created: 8/8/17
Purpose: This script is for counting the number of alignments to each known miRNA from the mirBase mature.fa post alignment by bowtie2.
"""

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-s', '--samfile', dest='sam',help="The filename and/or path of the sam file containing miRNA alignments.")
    parser.add_option( '-o', '--output-filename', dest='outfilename',help="The name of the output file containing the counts for each mature miRNA.")
    parser.add_option( '-a', '--alignment-stats-file', dest='align_stats',help="A file containing the stdout from running bowtie2.")
    ( options, args ) = parser.parse_args()

    if not options.sam:
        print "Please provide a sam file with the alignments from the miRNA sequencing."
        sys.exit(1)
    if not options.align_stats:
        print "Please provide the filename containing the stdout from bowtie2 alignment for running this QC."
        sys.exit(1)
    if not options.outfilename:
        options.outfilename = options.sam.split('.sam')[0]+'.miRNA.counts.txt'
        print "No output filename detected, writing counts to",options.outfilename
    return options

def qcMetrics(samfile,outfilename):
    counts = {}
    # For unaligned reads
    counts['*'] = 0
    with open(samfile,'r') as fh:
        line = fh.readline()
        while "@" in line[0]:
            if "@SQ" in line[:3]:
                mat_mirna = line.split()[1].split('SN:')[1]
                counts[mat_mirna] = 0
            line = fh.readline()
        print "Done with Header."
        while line:
            mat_mirna = line.split()[2]
            counts[mat_mirna] += 1
            line = fh.readline()
    unaligned = counts['*']
    del counts['*']
    outfile = open(outfilename,'w')
    outfile.write('Mature_miRNA,Count\n')
    total_mature = 0
    mapped_mature = 0
    for key in counts:
        if counts[key] > 0:
            mapped_mature += 1
        total_mature += 1
        outfile.write(key+','+str(counts[key])+'\n')
    outfile.write("Unaligned,"+str(unaligned)+'\n')
    outfile.close()
    print "Wrote miRNA counts to",outfilename
    return total_mature,mapped_mature

STATS = {0:'Total Reads',2:'Unaligned',3:'Aligned 1 Time',4:'Aligned > 1 Time',1:None,5:None}

def overallQCMetrics(total,found,alignfilename,stats=STATS):
    path = os.path.dirname(alignfilename)
    outfilename = os.path.join(path,'metrics.tsv')
    outfile = open(outfilename,'w')
    header = "Total Reads\tUnaligned\tAligned 1 Time\tAligned > 1 Time\tTotal Mature miRNAs\tNumber Mapped to Mature miRNAs\tPercent Mapped to mature miRNAs\n"
    outfile.write(header)
    out_stats = ''
    with open(alignfilename,'r') as fh:
        contents = fh.readlines()
        for i in range(len(contents)):
            sep = contents[i].split()
            if stats[i] and len(sep) > 3:
                out_stats += sep[0]+'\t'
    out_stats += str(total)+'\t'+str(found)+'\t'+str((float(found)/int(total))*100)+'\n'
    outfile.write(out_stats)
    outfile.close()
    print "Wrote QC metrics to",outfilename

def main():
    options = argParser()
    samfile = options.sam
    alignment_stats = options.align_stats
    outfilename = options.outfilename
    total,mapped = qcMetrics(samfile,outfilename)
    overallQCMetrics(total,mapped,alignment_stats)

if __name__ == '__main__':
    main()

