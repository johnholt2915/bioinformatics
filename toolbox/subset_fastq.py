import math
import random
import re
import sys
import os

def getRandomLines(filename,n_lines):
    filesize = os.path.getsize(filename)
    #print "filesize:",filesize
    fastq = open(filename,'r')
    fasta_lines = {}
    size_per_read = 0
    count = 0
    for line in fastq:
        if count == 4:
            break
        size_per_read += len(line)
        count += 1
    #print "size_per_read:",size_per_read
    fasta_lines = {}
    count = 0
    output = []
    while count < n_lines:
        offset = random.randrange(filesize)
        #print "initial offset:",offset
        if offset >= size_per_read:
            offset = offset - size_per_read
        else:
            offset = 0
        #print "offset:",offset
        fastq.seek(offset)
        found = False
        while not found:
            line = fastq.readline()
            #print line.strip()
            if line[0] != '@':
                continue
            split = line.split("cluster")
            if len(split) == 2:
                #print "FOUND LINE!!!",count
                found = True
                try:
                    assert fastq_lines[line]
                    print "JUSTKIDDING!!!",count
                except:
                    output.append(line)
                    for i in range(3):
                        output.append(fastq.readline())
                    fasta_lines[line] = None
                    #output.extend([seq_id,fasta_lines[seq_id]])
                    count += 1
    fastq.close()
    fasta_file_out = outputFileName(filename,"fa")
    print fasta_file_out
    #for line in output:
    #    print line.strip()
    fh = open(fasta_file_out,'w')
    fh.write("".join(output))
    fh.close()
    #return "".join(output)

def main(argv):
    filename = argv[1]
    n_lines = int(argv[2])
    getRandomLines(filename,n_lines)

if __name__ == '__main__':
    main(sys.argv)

