import gzip
import optparse

def fq2fa(fastq,outfilename,gzipped):
    if gzipped:
        fh = gzip.open(fastq,'rb')
    else:
        fh = open(fastq,'r')
    outfile = open(outfilename,'w')
    count = 0
    for line in fh:
        if count % 4 == 0:
            outfile.write('>'+line[1:])
        elif count % 4 == 1:
            outfile.write(line)
        count += 1
    fh.close()
    outfile.close()
    print "converted",fastq,"to",outfilename

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--fastq-in', dest='fastq')
    parser.add_option( '-o', '--fasta-out', dest='outfilename')
    parser.add_option( '-z', '--gzipped', dest='gzip', action='store_true', default=False)
    ( options, args ) = parser.parse_args()

    if not options.fastq:
        print "Please provide a fastq files as input"
    if not options.outfilename:
        if 'fastq' in options.fastq:
            options.outfilename = os.path.join(os.path.dirname(options.fastq),options.fastq.split('.fastq')[0]+'.fa')
        elif 'fq' in options.fastq:
            options.outfilename = os.path.join(os.path.dirname(options.fastq),options.fastq.split('.fq')[0]+'.fa')
    if options.gzip:
        if '.gz' not in options.fastq and '.gzip' not in options.fastq:
            print "file not recognized as gzipped, please gzip or change the file name to end in <fastq_filename>.gz"
            sys.exit(2)
    return options

def main():
    args = argParser()
    fastq = args.fastq
    outfilename = args.outfilename
    gzip = args.gzip
    fq2fa(fastq,outfilename,gzip)

if __name__ == '__main__':
    main()
