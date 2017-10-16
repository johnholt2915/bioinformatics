import sys

def filterReads(fq1,min_len,mismatches):
    fh1 = open(fq1,"r")
    fq1 = fq1.split('/')[-1]
    pass_filter_1_5prime = open("pass_filter_5prime_"+fq1,"w")
    pass_filter_no_trim_5prime = open("no_trim_pass_filter_5prime_"+fq1,"w")
    pass_filter_1_3prime = open("pass_filter_3prime_"+fq1,"w")
    pass_filter_no_trim_3prime = open("no_trim_pass_filter_3prime_"+fq1,"w")
    fail_filter_1 = open("fail_filter_"+fq1,"w")
    read1 = []
    no_trim_read1 = []
    # read a fastq entry (4 lines) at a time and save one copy for trimming and one for keeping the poly A/T stretch.
    for i in range(4):
        line = fh1.readline()
        read1.append(line)
        # used for keeping the original un-trimmed read.
        no_trim_read1.append(line)
    # read through the entire file and trim reads as appropriate
    while read1[0] != '':
        N_len = 0
        mismatch = 0
        nt_pos = 0
        # check for length of poly A/T regions through the recursive method 'passOrFailRead' 
        last_out_written,read1 = passOrFailRead(read1,N_len,min_len,mismatch,mismatches,nt_pos,rev=False)
        # if the read has a poly A/T stretch on its 5' end that is longer than 15 bp before 1 mismatch
        # and if the read has a relatively uniform distribution of AGCT after trimming
        # the read will pass.
        if last_out_written:
            pass_filter_1_5prime.write("".join(read1))
            pass_filter_no_trim_5prime.write("".join(no_trim_read1))
        else:
            # if the read failed from 5' side, then the 3' end will be tested in the same way.
            nt_pos = -1
            last_out_written,read1 = passOrFailRead(read1,N_len,min_len,mismatch,mismatches,nt_pos,rev=True)
            if last_out_written:
                pass_filter_1_3prime.write("".join(read1))
                pass_filter_no_trim_3prime.write("".join(no_trim_read1))
            else:
                # if the read fails both ends, then write the failed reads to failed_filter_1 file handle.
                fail_filter_1.write("".join(read1))
        for i in range(4):
            line = fh1.readline()
            read1[i] = line
            no_trim_read1[i] = line
    pass_filter_1_5prime.close()
    pass_filter_no_trim_5prime.close()
    pass_filter_1_3prime.close()
    pass_filter_no_trim_3prime.close()
    fail_filter_1.close()
    fh1.close()

def passOrFailRead(read1,n_len,min_len,mismatch,mismatches,nt_pos,rev=False):
    if rev:
        if read1[1][:-1][nt_pos] == 'A':
            #print "reverse, starting with A"
            nt_pos -= 1
            n_len += 1
            return polyStretchRev('A',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read1[1][:-1][nt_pos] == 'T':
            #print "reverse, starting with T"
            nt_pos -= 1
            n_len += 1
            return polyStretchRev('T',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read1[1][:-1][nt_pos] == 'G':
            nt_pos -= 1
            n_len += 1
            return polyStretchRev('G',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read1[1][:-1][nt_pos] ==  'C':
            nt_pos -= 1
            n_len += 1
            return polyStretchRev('C',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        else:
            #print "reverse, starting with N, C or G"
            if mismatch > mismatches:
                #print "read doesn't start with poly T/A"
                return False,read1
            mismatch += 1
            nt_pos -= 1
            return passOrFailRead(read1,n_len,min_len,mismatch,mismatches,nt_pos,rev=True)
    else:
        if read1[1][:-1][nt_pos] == 'A':
            #print "forward, starting with A"
            nt_pos += 1
            n_len += 1
            return polyStretch('A',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read1[1][:-1][nt_pos] == 'T':
            #print "forward, starting with T"
            nt_pos += 1
            n_len += 1
            return polyStretch('T',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read[1][:-1][nt_pos] == 'G':
            nt_pos += 1
            n_len += 1
            return polyStretch('G',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read[1][:-1][nt_pos] == 'C':
            nt_pos += 1
            n_len += 1
            return polyStretch('C',read1,min_len,n_len,nt_pos,mismatch,mismatches)
        else:
            #print "forward, starting with N, C or G"
            if mismatch > mismatches:
                #print "read doesn't start with polyT/A"
                return False,read1
            mismatch += 1
            nt_pos += 1
            return passOrFailRead(read1,n_len,min_len,mismatch,mismatches,nt_pos,rev=False)

NT_MAP= {'A':['G','C','T'],'T':['G','C','A'],'G':['C','T','A'],'C':['G','T','A']}

def polyStretchRev(poly_n,read1,min_len,N_len,nt_pos,mismatch,mismatches,nt_map=NT_MAP):
    #nt = poly_n
    #other_nt = ''
    other_nt = nt_map[poly_n]
    #if nt == 'A':
    #    other_nt = 'T'
    #else:
    #    other_nt = 'A'
    read = read1[1][:-1]
    possibilities = [poly_n]
    while mismatch <= mismatches and (len(read)-(-nt_pos-1)) >= min_len:
        #if read[nt_pos] == nt:
        if read[nt_pos] == poly_n:
            N_len += 1
            nt_pos -= 1
        else:
            mismatch += 1
            other_nt 
        #elif read[nt_pos] == other_nt and (-nt_pos-1) <= mismatches:
            N_len = 1
            nt_pos -= 1
            mismatch += 1
            nt = other_nt
            other_nt = poly_n
        elif read[nt_pos] == poly_n and mismatch <= mismatches:
            if N_len >= 2:
                N_len += 1
            else:
                N_len = 2
            nt_pos -= 1
            other_nt = nt
            nt = poly_n
        else:
            mismatch += 1
            nt_pos -= 1
    if (len(read) - (-nt_pos-1)) < min_len:
        #print 'not enough read left for alignment: min length after'
        return False,read1
    if N_len >= 15:
        #print "read passes!"
        nt_pos += 2
        read1[1] = read1[1][:nt_pos]+"\n"
        read1[3] = read1[3][:nt_pos]+"\n"
        return True,read1
    elif N_len < 15: # poly_n is too short...
        #print "poly_n stretch is less than 10"
        return False,read1

def polyStretch(poly_n,read1,min_len,N_len,nt_pos,mismatch,mismatches):
    nt = poly_n
    other_nt = ''
    if nt == 'A':
        other_nt = 'T'
    else:
        other_nt = 'A'
    read = read1[1][:-1]
    while mismatch <= mismatches and (len(read)-nt_pos) >= min_len:
        #print "___init___"
        #print "N_len:",N_len,"nt_pos:",nt_pos,"read[nt_pos]:",read[nt_pos],"nt:",nt,'other_nt:',other_nt,"mismatch:",mismatch,"mismatches:",mismatches,"poly_n",poly_n
        if read[nt_pos] == nt:
            #print "read[nt_pos] == nt"
            N_len += 1
            nt_pos += 1
        elif read[nt_pos] == other_nt and nt_pos <= mismatches:
            #print "read[nt_pos] == other_nt and nt_pos < mismatches"
            N_len = 1
            nt_pos += 1
            mismatch += 1
            nt = other_nt
            other_nt = poly_n
        elif read[nt_pos] == poly_n and mismatch <= mismatches:
            #print "read[nt_pos] == poly_n and nt_pos <= mismatches"
            if N_len >= 2:
                N_len += 1
            else:
                N_len = 2
            nt_pos += 1
            other_nt = nt
            nt = poly_n
        else:
            #print "mismatch"
            mismatch += 1
            nt_pos += 1
    if (len(read) - nt_pos) < min_len:
        #print 'not enough read left for alignment: min length after '
        return False,read1
    if N_len >= 15:
        #print "read passes!"
        read1[1] = read1[1][nt_pos:]
        read1[3] = read1[3][nt_pos:]
        return True,read1
    elif N_len < 15: # poly_n is too short...
        #print "poly_n stretch is less than 10"
        return False,read1

def main():
    fq1 = sys.argv[1]
    min_len = int(sys.argv[2])
    mismatches = 0
    try:
        mismatches = int(sys.argv[3])
    except:
        mismatches = 1
    filterReads(fq1,min_len,mismatches)

if __name__ == '__main__':
    main()

