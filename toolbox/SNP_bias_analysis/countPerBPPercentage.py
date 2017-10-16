import sys
import json


CIGAR_SET = set(["I","D","N","S","H","P","=","X"])

def makeIndex():
    chrs_dict = {}
    chrs = ['X','Y','M']
    chrs.extend(range(1,23,1))
    for i in chrs:
        chrs_dict['chr'+str(i)] = dict()
    return chrs_dict

def makeOutfile(fileIn):
    split = fileIn.split(".mpileup")[0]
    return split+".txt"

def readHasPoorNTsDistribution(read):
    nts = {'A':0,'G':0,'N':0,'C':0,'T':0}
    for nt in read:
        nts[nt] += 1
    if (len(read) - nts['A']) < (len(read)/2):
        return True
    elif (len(read) - nts['T']) < (len(read)/2):
        return True
    elif (len(read) - nts['G']) < (len(read)/2):
        return True
    elif (len(read) - nts['C']) < (len(read)/2):
        return True
    else:
        return False

def cigarPasses(cigar,set=CIGAR_SET):
    return any(elem in set for elem in cigar)

def bitwiseReverseCalculator(flag):
    if flag == 0:
        return False
        #print "bitwise flag is 0."
    else:
        bit = 0
        while flag > 0:
            flag_bit = flag%2
            if flag_bit == 1:
                if 2**bit == 16:
                    return True
            bit += 1
            flag = flag/2
        return False

def grabAlignmentsPerCoverageThreshold(sam,cov,headerfilename):
    with open(headerfilename,'r') as fh:
        header = fh.read()
    index = makeIndex()
    revcomp_index = makeIndex()
    with open(sam,"r") as fh:
        line = fh.readline()
        sep = []
        cigar = ''
        cig_pass = False
        lines_to_maybe_write = ''
        try:
            while line != None:
                #print line
                sep = line.split()
                cigar = sep[5]
                cig_fail = cigarPasses(cigar)
                # make sure line has all matches in cigar
                while cig_fail:
                #    print "cigar fail!",cigar,cig_fail,"\n"
                    line = fh.readline()
                #    print line
                    sep = line.split()
                    cigar = sep[5]
                    cig_fail = cigarPasses(cigar)
                #print "cigar pass!",cigar,cig_fail
                # all matched read
                # now count how many reads have all matched alignments starting in same position as
                # the read found with all matches
                if readHasPoorNTsDistribution(sep[9]):
                #    print "read with poor NTs"
                #    print sep[9]
                    line = fh.readline()
                    continue
                #print "read has good NTs"
                #print sep[9]
                if len(sep[9]) > 70:
                #    print "line too long"
                    line = fh.readline()
                    continue
                # now save info about this line, and see how many lines map to same location...
                #print "line length is good",len(sep[9])
                if bitwiseReverseCalculator(int(sep[1])):
                    pos = int(sep[3]) + len(sep[9])-1
                    try:
                        revcomp_index[sep[2]][str(pos)].append(line)
                    except:
                        revcomp_index[sep[2]][str(pos)] = [line]
                    line = fh.readline()
                    continue
                chr = sep[2]
                begin_pos = sep[3]
                #init_strand = sep[1]
                tmp_list = [sep[9]+";"+sep[10]]
                #plus_minus = {'+':0,'-':0}
                #print chr,begin_pos,"\n"
                #print "adding line to line_str"
                line_to_maybe_write = header + line
                #print line_to_maybe_write
                # begin to check for other lines which map to same location
                line = fh.readline()
                sep = line.split()
                start = sep[3]
                #print "start=",start
                while start == begin_pos: # while current line has same start position as the first line to pass cigar and NT filters ...
                    cigar = sep[5]
                    cig_fail = cigarPasses(cigar)
                    # check this cigar
                    if cig_fail:
                        # if fail, grab next line to see if it still has same start coordinates
                        line = fh.readline()
                        sep = line.split()
                        start = sep[3]
                #        print "cig_fail"
                        continue
                #    print "adding",sep[9],"to list"
                    tmp_list.append(sep[9]+";"+sep[10])
                #    print "length tmp_list=",len(tmp_list)
                #    print "adding line to line_str"
                    line_to_maybe_write = line_to_maybe_write + line
                    line = fh.readline()
                    sep = line.split()
                    start = sep[3]
                #    print "start:",start
                # check the tmp list for sufficient coverage (seeing if the first 20 bps will have more than cov reads aligned to it...)
                if len(tmp_list) < cov:
                    continue
                #print "pass"
                #print "tmp_list length:",len(tmp_list)
                # write the length of tmp_list to the first entry of the list,
                # this is to check that a reference was found for this start coordinate later on in the code
                tmp_list.insert(0,len(tmp_list))
                index[chr][begin_pos] = tmp_list
                next_sam = open("sam_alignments_"+chr+"_"+str(begin_pos)+".sam","w")
                next_sam.write(line_to_maybe_write)
                next_sam.close()
                #print "end"
        except:
            pass
    print "writing reverse complemented alignments..."
    for chr in revcomp_index:
        for key in revcomp_index[chr].keys():
            length = len(revcomp_index[chr][key])
            if length >= cov:
                next_sam = open("sam_alignments_"+chr+"_"+str(key)+"_revcomp.sam","w")
                next_sam.write(header+"".join(revcomp_index[chr][key]))
                next_sam.close()
                revcomp_index[chr][key].insert(0,length)
            else:
                del revcomp_index[chr][key]
    return index,revcomp_index

def dupCountedSam(counted_dict,chr,pos):
    overlap = False
    for i in range(int(pos)-70,int(pos)+1,1):
        try:
            assert counted_dict[chr][str(i)]
            overlap = True
            return overlap
        except:
            continue
    return overlap

def referenceFound(index):
    outfile = open("bad_lines.txt","a")
    outfile.write("All below regions have no reference")
    count = 0
    for i in index:
        for j in index[i]:
            if index[i][j][0] == len(index[i][j][1:]):
                outfile.write(str(i)+" "+str(j)+"\n")
                count += 1
    if count > 0:
        print "Some regions don't have a reference from the mpileup file provided."
        print "Look at the bad_lines.txt file to see which regions are bad"
        sys.exit()

def retrieveReference(mpileup,index,revcomp_index):
    size_ref = 20
    adding_ref = False
    ref_string = ''
    begin = ''
    pileup_dict = makeIndex()
    with open(mpileup,"r") as fh:
        for line in fh:
            sep = line.split()
            chr = sep[0]
            coord = sep[1]
            ref = sep[2]
            pileup_dict[chr][coord] = str.upper(ref)
    count = 0
    print "finding references for index"
    for chr in index:
        for pos in index[chr]:
            ref_str = ''.join(pileup_dict[chr][str(i)] for i in range(int(pos),int(pos)+20,1))
            index[chr][pos].append(ref_str)
            count  += 1
    count = 0
    print "finding references for rev_index"
    for chr in revcomp_index:
        for pos in revcomp_index[chr]:
            ref_str = ''.join(pileup_dict[chr][str(i)] for i in range(int(pos),int(pos)-20,-1))
            revcomp_index[chr][pos].append(ref_str)
            count += 1
    return index,revcomp_index

NT_DICT = {'A':'T','T':'A','G':'C','C':'G'}

def reverseComplement(read,nt_dict=NT_DICT):
    rev = ''
    for nt in read:
        rev = nt_dict[nt]+rev
    return rev

def countBasesAndWriteTSV(mpileup,index,revcomp_index,outfilename):
    outfile = open(outfilename,"a")
    outfile.write("Chromo\tPosition\tReference\tPercent_Match\tPercent_1st_Non_Reference_Mismatch\tA\tG\tC\tT\tTotal_Coverage\tGood_Quality\n")
    print "Chromo\tPosition\tReference\tPercent_Match\tPercent_1st_Non_Reference_Mismatch\tA\tG\tC\tT\tTotal_Coverage\tGood_Quality\n"
    nts = ['A','G','C','T']
    print "writing tsv for index"
    for chr in index:
        for pos in index[chr]:
            write_line = ntFrequencies(index[chr][pos],chr,pos,nts)
            #print write_line
            outfile.write(write_line)
    print "writing tsv for rev_index"
    for chr in revcomp_index:
        for pos in revcomp_index[chr]:
            write_line = revcompNTFrequencies(revcomp_index[chr][pos],chr,pos,nts)
            outfile.write(write_line)
    outfile.close()

def revcompNTFrequencies(read_list,chr,pos,nts):
    ref = read_list[-1]
    nt_dict = {'A':[],'G':[],'C':[],'T':[],'M':[],'mM':[],'TC':[],'GQ':[]}
    for i in range(1,len(read_list[:-1])):
        sep = read_list[i].split()
        read_list[i] = sep[9]+";"+sep[10]
    for i in range(1,len(ref)+1,1):
        for key in nt_dict:
            if 'M' not in key:
                nt_dict[key].append(0)
            else:
                nt_dict[key].append(0.0)
        total = 0
        unfiltered_total = 0
        for seq in read_list[1:-1]:
            sep = seq.split(';')
            unfiltered_total += 1
            if sep[0][-i] != 'N' and (ord(sep[1][-i])-33) >= 20:
                total += 1
                nt_dict[sep[0][-i]][i-1] += 1
        freq = 0
        for nt in nts:
            if nt == ref[i-1]:
                continue
            elif nt_dict[nt][i-1] >= freq:
                freq = nt_dict[nt][i-1]
        nt_dict['GQ'][i-1] = total
        nt_dict['TC'][i-1] = unfiltered_total
        if total > 0:
            nt_dict['M'][i-1] = float(nt_dict[ref[i-1]][i-1])/total
            nt_dict['mM'][i-1] = float(freq)/total
        else:
            nt_dict['M'][i-1] = float(nt_dict[ref[i-1]][i-1])/1
            nt_dict['mM'][i-1] = float(freq)/1
    A = ";".join(str(x) for x in nt_dict['A'])
    G = ";".join(str(x) for x in nt_dict['G'])
    C = ";".join(str(x) for x in nt_dict['C'])
    T = ";".join(str(x) for x in nt_dict['T'])
    M = ";".join(str(x) for x in nt_dict['M'])
    mM = ";".join(str(x) for x in nt_dict['mM'])
    TC = ";".join(str(x) for x in nt_dict['TC'])
    GQ = ";".join(str(x) for x in nt_dict['GQ'])
    write_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,M,mM,A,G,C,T,TC,GQ)
    return write_line

def ntFrequencies(read_list,chr,pos,nts):
    ref = read_list[-1]
    nt_dict = {'A':[],'G':[],'C':[],'T':[],'M':[],'mM':[],'TC':[],'GQ':[]}
    #print chr,pos
    #print "Reference!:",ref
    #print "read_list:",read_list
    for i in range(len(ref)):
        for key in nt_dict:
            if 'M' not in key:
                nt_dict[key].append(0)
            else:
                nt_dict[key].append(0.0)
        #print "nt_dict",nt_dict
        total = 0
        unfiltered_total = 0
        for seq in read_list[1:-1]:
            sep = seq.split(';')
            unfiltered_total += 1
            if sep[0][i] != 'N' and (ord(sep[1][i])-33) >= 20:
                total += 1
                nt_dict[sep[0][i]][i] += 1
        freq = 0
        for nt in nts:
            if nt == ref[i]:
                continue
            elif nt_dict[nt][i] >= freq:
                freq = nt_dict[nt][i]
        nt_dict['TC'][i] = unfiltered_total
        nt_dict['GQ'][i] = total
        if total > 0:
            nt_dict['M'][i] = float(nt_dict[ref[i]][i])/total
            nt_dict['mM'][i] = float(freq)/total
        else:
            nt_dict['M'][i] = float(nt_dict[ref[i]][i])/1
            nt_dict['mM'][i] = float(freq)/1
    A = ";".join(str(x) for x in nt_dict['A'])
    G = ";".join(str(x) for x in nt_dict['G'])
    C = ";".join(str(x) for x in nt_dict['C'])
    T = ";".join(str(x) for x in nt_dict['T'])
    M = ";".join(str(x) for x in nt_dict['M'])
    mM = ";".join(str(x) for x in nt_dict['mM'])
    TC = ";".join(str(x) for x in nt_dict['TC'])
    GQ = ";".join(str(x) for x in nt_dict['GQ'])
    write_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,M,mM,A,G,C,T,TC,GQ)
    return write_line

def main():
    sam = sys.argv[1]
    mpileup = sys.argv[2]
    cov = int(sys.argv[3])
    # filename to write the script output.
    outfile = sys.argv[4]
    # file to tack on sam file lines that pass ...
    append_file = sys.argv[5]
    json_file = outfile.split('.txt')[0] + '.json'
    print "sam",sam
    print "mpileup",mpileup
    print "cov",cov
    print "outfile",outfile
    print "append_file",append_file
    try:
        with open(json_file,"r") as fh:
            print "openned",json_file,"..."
            json_str = fh.read()
            index = json.loads(json_str)
            print "loaded json"
        with open("rev_"+json_file,"r") as fh:
            print "openned rev_"+json_file+"..."
            json_str = fh.read()
            rev_index = json.loads(json_str)
            print "loaded rev json"
    except Exception as e:
        print "could not find json, making sam index..."
        index,rev_index = grabAlignmentsPerCoverageThreshold(sam,cov,append_file)
        print "writing index and rev_index to json files..."
        with open(json_file,"w") as fh:
            json.dump(index,fh)
        with open("rev_"+json_file,"w") as fh:
            json.dump(rev_index,fh)
        #for i in index:
        #    print i,len(index[i])
        #    for pos in index[i]:
        #        print pos,len(index[i][pos])
        #for i in rev_index:
        #    print i,len(index[i])
        #    for pos in index[i]:
        #        print pos,len(index[i][pos])
        # loop through mpileup and get reference at bp which are in the index
        ## both 1 based coordinates
    print "retrieving reference for index and rev_index..."
    index,rev_index = retrieveReference(mpileup,index,rev_index)
    print "validating index references..."
    referenceFound(index)
    print "validating rev_index references..."
    referenceFound(rev_index)
    
    # count numbers of matches/mismatches per region and write to file
    print "counting bases and writing tsv"
    countBasesAndWriteTSV(mpileup,index,rev_index,outfile)

if __name__ == '__main__':
    main()


