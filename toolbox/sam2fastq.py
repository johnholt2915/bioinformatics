import os 
import sys

def main(argv):
    sam_input_file = argv[1]
    outfilename = sam_input_file.split('.sam/')[0]+'.fastq'
    outfile = open(outfilename,'w')
    with open(sam_input_file,'r') as fh:
        for line in fh:
            if line[0] != '@':
                sep = line.split('\t')
                outfile.write('@'+sep[0]+'\n')
                outfile.write(sep[9]+'\n')
                outfile.write('+\n')
                outfile.write(sep[10]+'\n')
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
