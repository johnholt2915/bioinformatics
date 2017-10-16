import gzip
import optparse
import sys
import os
import random


def get_random_lines(filename,outfilename,n_lines,gzipped=False):
    filesize = os.path.getsize(filename)
    if gzipped:
        fastq = gzip.open(filename,'rb')
    else:
        fastq = open(filename,'r')
    outfile = open(outfilename,'w')
    count = 0
    found_reads = {}
    while count < n_lines:
        offset = random.randrange(filesize)
        fastq.seek(offset)
        fastq.readline()
        found = False
        while not found:
            found = True
            line = fastq.readline()
            split = line.split(":")
            if len(split) == 10 and split[0][0] == '@': # if the line is the 1st of a fastq read ... 
                try: # if we've seen it before, then continue searching ... NO REPEATS!
                    assert found_reads[line[:-1]]
                except: # otherwise write the read to the outfile, increment the count
                    found_reads[line[:-1]] = line[1:][:-1]
                    outfile.write(line)
                    for i in range(3):
                        outfile.write(fastq.readline())
                    count += 1

def extract_fq_read(fastq,readname,outfilename,gzipped):
    if os.path.exists(readname): # if readname is a file...
        rn_fh = open(readname,'r')
        readname = [rn.split('\n')[0] for rn in rn_fh]
    else: # otherwise this is just a string, so find this read in the fastq file
        readname = [readname]
    if gzipped:
        fh = gzip.open(fastq,'rb')
    else:
        fh = open(fastq,'r')
    outfile = open(outfilename,'w')
    line = fh.readline()
    len_readname = len(readname)
    found = 0
    print ""
    print "readname",readname
    print "fastq_format",[line[1:][:-1]]
    cont = raw_input('Does the readname format match the format in the fastq file?\nPress enter if they do, and type any character if they do not.\n>>> ')
    if cont != '':
        print "exiting..."
        sys.exit()
    while line:
        if line[1:][:-1] in readname:
            outfile.write(line)
            for i in range(3):
                outfile.write(fh.readline())
            line = fh.readline()
            found += 1
            if found == len_readname:
                break
        else:
            for i in range(4):
                line = fh.readline()
    fh.close()
    outfile.close()
    if found == len_readname:
        print "found all reads!"
    else:
        print "found",found,"of",len_readname,"reads"

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--fastq-in', dest='fastq')
    parser.add_option( '-r', '--read-name',dest='readname',help="This can be a read name for a single read, or a file with each line containing a specific readname.")
    parser.add_option( '-o', '--fasta-out', dest='outfilename')
    parser.add_option( '-z', '--gzipped', dest='gzip', action='store_true', default=False)
    parser.add_option( '-a', '--rand_reads',dest='rand_reads', action='store_true',default=False)
    parser.add_option( '-n', '--n-reads',dest='n_reads')
    ( options, args ) = parser.parse_args()

    if not options.fastq:
        print "Please provide a fastq files as input."
        sys.exit(2)
    if not options.readname and not options.rand_reads:
        print "Please provide a read name to extract from input fastq."
        sys.exit(2)
    if not options.outfilename:
        if 'fastq' in options.fastq:
            options.outfilename = os.path.join(os.path.dirname(options.fastq),options.fastq.split('.fastq')[0]+'.fa')
        elif 'fq' in options.fastq:
            options.outfilename = os.path.join(os.path.dirname(options.fastq),options.fastq.split('.fq')[0]+'.fa')
    if options.gzip:
        if '.gz' not in options.fastq and '.gzip' not in options.fastq:
            print "file not recognized as gzipped, please gzip or change the file name to end in <fastq_filename>.gz"
            sys.exit(2)
    if options.rand_reads and not options.n_reads:
        print "please provide the number of random reads you would like from the fastq file."
        sys.exit(2)
    return options

def main():
    args = argParser()
    fastq = args.fastq
    readname = args.readname
    outfilename = args.outfilename
    gzipped = args.gzip
    rand_reads = args.rand_reads
    n_reads = args.n_reads
    if rand_reads:
        get_random_lines(fastq,outfilename,n_reads,gzipped)
    else:
        extract_fq_read(fastq,readname,outfilename,gzipped)

if __name__ == '__main__':
    main()
