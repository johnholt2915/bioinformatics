import csv
import sys

def removePolyN(filename,poly,mismatches,rev=False):
    nt_pos = 0
    end_T = 0
    if rev:
        nt_pos = -1
        end_T = -1
    len_hist = dict()
    outfile = open("no_poly_T_"+".".join(filename.split('.')[:-1])+".fastq","a")
    fh = open(filename,'r')
    read = []
    for i in range(4):
        read.append(fh.readline())
    while read[0] != '':
        mismatch = 0
        read = passOrFailRead(read,mismatch,mismatches,nt_pos,rev=rev)

def passOrFailRead(read,mismatch,mismatches,nt_pos,rev=False):
    #print line
    #n_len = 0
    if rev:
        if read[1][:-1][nt_pos] == 'A':
            #print "reverse, starting with A"
            nt_pos -= 1
            #n_len += 1
            return polyStretchRev('A',line,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read[1][:-1][nt_pos] == 'T':
            #print "reverse, starting with T"
            nt_pos -= 1
            #n_len += 1
            return polyStretchRev('T',line,min_len,n_len,nt_pos,mismatch,mismatches)
        else:
            #print "reverse, starting with N, C or G"
            if mismatch > mismatches:
                #print "read doesn't start with poly T/A"
                return False
            mismatch += 1
            nt_pos -= 1
            return passOrFailRead(read,mismatch,mismatches,nt_pos,rev=rev)
    else:
        if read[1][:-1][nt_pos] == 'A':
            #print "forward, starting with A"
            nt_pos += 1
            #n_len += 1
            return polyStretch('A',line,min_len,n_len,nt_pos,mismatch,mismatches)
        elif read[1][:-1][nt_pos] == 'T':
            #print "forward, starting with T"
            nt_pos += 1
            #n_len += 1
            return polyStretch('T',line,min_len,n_len,nt_pos,mismatch,mismatches)
        else:
            #print "forward, starting with N, C or G"
            if mismatch > mismatches:
                #print "read doesn't start with polyT/A"
                return False
            mismatch += 1
            nt_pos += 1
            return passOrFailRead(read,mismatch,mismatches,nt_pos,rev=rev)

def polyStretchRev(poly_n,read,min_len,N_len,nt_pos,mismatch,mismatches):
    #N_len = init_n_len
    nt = poly_n
    other_nt = ''
    if nt == 'A':
        other_nt = 'T'
    else:
        other_nt = 'A'
    #nt_pos = init_pos
    #mismatch = init_mismatch
    while mismatch <= mismatches and (len(read)-(-nt_pos-1)) >= min_len:
        if read[nt_pos] == nt:
            N_len += 1
            nt_pos -= 1
        elif read[nt_pos] == other_nt and (-nt_pos-1) <= mismatches:
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
        return False
    if N_len >= 15:
        #print "read passes!"
        return True
    elif N_len < 15: # poly_n is too short...
        #print "poly_n stretch is less than 10"
        return False

def polyStretch(poly_n,read,min_len,N_len,nt_pos,mismatch,mismatches):
    #N_len = init_n_len
    nt = poly_n
    other_nt = ''
    if nt == 'A':
        other_nt = 'T'
    else:
        other_nt = 'A'
    #nt_pos = init_pos
    #mismatch = init_mismatch
    #print poly_n
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
        return False
    if N_len >= 15:
        #print "read passes!"
        return True
    elif N_len < 15: # poly_n is too short...
        #print "poly_n stretch is less than 10"
        return False

 
"""
    with open(filename,"r") as fh:
        count = 0
        end_T = 0
        for line in fh:
            if count % 4 == 1:
                for i in range(len(line)):
                    if line[i] == poly:
                        end_T += 1
                        continue
                    else:
                        outfile.write(line[i:])
                        end_T = i
                        break
            elif count % 4 == 3:
                outfile.write(line[end_T:])
                end_len = 83 - (end_T+1)
                try:
                    len_hist[str(end_len)] += 1
                except KeyError:
                    len_hist[str(end_len)] = 1
                end_T = 0
            else:
                outfile.write(line)
            count += 1
    dict2csv(len_hist,"lengths_after_poly"+poly+"_rm","lengths")

def removePolyN(filename,poly,basename):
    len_hist = dict()
    outfile = open("no_poly_T_"+basename+".fastq","a")
    with open(filename, "r") as fh:
        count = 0
        ind = 0
        for line in fh:
            #print line
            if count % 4 == 1:
                missed_one = False
                end_T = 0
                for i in range(len(line)):
                    #print line[i]
                    if line[i] == poly:
                        if missed_one:
                            end_T += 1
                        continue
                    elif not missed_one:
                        missed_one = True
                        continue
                    else:
                        missed_one = False
                        if end_T > 0:
                            #print line[i:]
                            outfile.write(line[i:])
                            ind = i
                        else:
                            #print line[i-1:]
                            outfile.write(line[i-1:])
                            ind = i-1
                        end_T = 0
                        break
            elif count % 4 == 3:
                #print line[ind:]
                outfile.write(line[ind:])
                end_len = 83 - (ind+1)
                try:
                    len_hist[str(end_len)] += 1
                except KeyError:
                    len_hist[str(end_len)] = 1
                ind = 0
            else:
                outfile.write(line)
            count += 1
    dict2csv(len_hist,"lengths_after_poly"+poly+"_rm","lengths")
"""

def dict2csv(mydict,name,x):
    with open(name+'.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([x,"frequency"])
        for key, value in mydict.items():
            writer.writerow([key, value])

def main():
    try:
        filename = sys.argv[1]
        poly_n = sys.argv[2]
        mismatches = sys.argv[3]
        rev = sys.argv[4]
        if rev == 'true':
            rev = True
        elif rev == 'false':
            rev = False
    except:
        usage()
        sys.exit()
    removePolyN(filename,poly_n,mismatches,rev)

if __name__ == '__main__':
    main()

