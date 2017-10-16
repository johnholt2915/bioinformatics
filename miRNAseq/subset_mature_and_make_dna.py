import os
import sys
import optparse
import re

RNA_DICT = {'A':'U','C':'G','G':'C','T':'A','U':'A','N':'N'}
DNA_DICT = {'A':'T','C':'G','G':'C','T':'A','U':'A','N':'N'}
NT_DICTS = [DNA_DICT,RNA_DICT]

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-f', '--fasta', dest='fasta',help="Fasta must be a fasta file from mirBase.")
    parser.add_option( '-s', '--subset', dest='subset',help="Filter must be in the format 'Genus species'.")
    parser.add_option( '-o', '--output-filename', dest='outfilename',help="Output filename can be full path or filename. If a filename, make sure you want data in local directory.")
    ( options, args ) = parser.parse_args()
    if not options.fasta:
        print help
        sys.exit()
    if not options.subset:
        print "No subset filter detected. Returning all reads."
        options.subset = '>'
    if not options.outfilename:
        print "No output filename detected. Using default subset.fa as output."
        options.outfilename = 'subset.fa'
    return options

def readConverter(read,RT=False,to_RNA=False,BP_DICT=NT_DICTS):
    """This function returns the nucleotide sequence passed as input, as either the reverse transcript of the read, or the DNA version of the read."""
    read = read.split('\n')[0]
    if RT:
        if to_RNA:
            new_read = [BP_DICT[1][nt.upper()] for nt in read][::-1]
        else:
            new_read = [BP_DICT[0][nt.upper()] for nt in read][::-1]
        return new_read
    elif to_RNA:
        return re.sub('T','U',read)
    else:
        return re.sub('U','T',read)

def isolateReads(mat_fa,filter,filename):
    outfile = open(filename,'w')
    with open(mat_fa,'r') as fh:
        line = fh.readline()
        while line:
            if '>' in line and filter in line:
                outfile.write(line)
                line = fh.readline()
                while line and '>' not in line:
                    read = readConverter(line.split('\n')[0]) + '\n'
                    outfile.write(read)
                    line = fh.readline()
            else:
                line = fh.readline()
    outfile.close()
    print "Wrote fasta reads to",filename,"using",filter,"to filter fasta lines."

def main():
    options = argParser()
    fasta = options.fasta
    subset_filter = options.subset
    output_filename = options.outfilename
    isolateReads(fasta,subset_filter,output_filename)

if __name__ == '__main__':
    main()

