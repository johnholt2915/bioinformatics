"""
This script is used to test a list of barcodes (of any length) for a variety of criteria.
The barcode list must in in a text file with a header, where each line is a new barcode.

Example:
Header
ACGTACGT
TGCATGCA
ATCGATCG
GTCAGTCA
.
.
.

The follow tests are available:
1) homopolymer stretch of longer than a (default = 2) specified length.
2) test all barcodes to see if they collide with other barcodes in the set in terms of base mismatch. The users can specify the distance in terms of base mismatch:
    Two barcodes with 1 base mismatch have at least 3 bases different.
3) substring match test in forward or reverse orientation.

"""
import os
import sys
from itertools import product
import optparse
from random import shuffle


def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-b', '--barcode-txt-file', dest='barcode_txt_file',default='')
    parser.add_option( '-p', '--homopolymer-test', dest='homopolymer_test',action='store_true',default=False)
    parser.add_option( '-s', '--substring-test', dest='substring_test',action='store_true',default=False)
    parser.add_option( '-u', '--universal-barcode-set', dest='universal_barcode_set',default='')
    parser.add_option( '-n', '--n-length-barcode', dest='n_length_barcode',default=8)
    parser.add_option( '-f', '--substring-test_filter',dest='substring_test_filter',default=5)
    parser.add_option( '-m', '--mismatch-test', dest='mismatch_test',action='store_true',default=False)
    parser.add_option( '-i', '--n-base-mismatch',dest='n_base',default=1)
    parser.add_option( '-t', '--homopolymer-thresh', dest='homopolymer_thresh',default=2,help='homopolymer segments longer than this value merits disgard of a barcode.')
    parser.add_option( '-o', '--output-file', dest='output_file',default='')
    parser.add_option( '-a', '--target-bcs-count', dest="target_bc_count",default=80)
    parser.add_option( '-z', '--seed-search-analysis', dest='seed_search_analysis',action='store_true',default=False)
    parser.add_option( '-c', '--remove-commonality-method', dest='remove_commonality_method',default=False,action='store_true')
    parser.add_option( '-r', '--rev-comp',dest='rev_comp',default=False,action='store_true')
    
    ( options, args ) = parser.parse_args()

    if not options.mismatch_test and not options.homopolymer_test and not options.substring_test:
        print "Exiting ... "
        sys.exit()
    return options

#def gen_barcodes(n_base):
#    return product('AGCT',repeat=n_base)

REV_COMP = {'A':'T','T':'A','G':'C','C':'G'}
def acceptable(bc1,bc_2,n_base,rc=False,rc_dict=REV_COMP):
    if rc:
        bc2 = "".join([rc_dict[nt] for nt in bc_2][::-1])
    else:
        bc2 = bc_2
    cmbs1 = [bc1[j:(j+i)] for i in range(n_base,n_base+1) for j in range(len(bc1)-i+1)]
    return not any([(bc2[j:(j+i)] in cmbs1) for i in range(n_base,len(bc2)) for j in range(len(bc2)-i+1)])

def len_homopolymer(bc,length):
    history = ''
    last = ''
    for nt in bc:
        if nt == last:
            if len(history) == length:
                return True
            else:
                history += nt
        else:
            last = nt
            history = nt
    if len(history) > 2:
        return True
    return False

def substring_test(bc,bcs,n_base,return_bool=False,return_all=False,rev_comp=False):
    if not return_all:
        failures = ''
        for b in bcs:
            if b != bc:
                if not acceptable(bc,b,n_base,rc=rev_comp):
                    failures += b+'|'
        if len(failures) > 0:
            if return_bool:
                return False
            else:
                return 'Failed,'+failures
        if return_bool:
            return True
        return 'Passed,NA'
    else:
        results = []
        for b in bcs:
            for base in range(8,0,-1):
                if not acceptable(bc,b,base,rc=rev_comp):
                    results.append(str(base))
                    break
        return ",".join(results)

class pseudo_generator:
    def __init__(self,txt_file,n=0,header=True):
        self.txt_file = txt_file
        self.n = n
        self.header = header
        if txt_file != '':
            self.fh = open(txt_file,'r')
            if header:
                self.header = self.fh.readline().split('\n')[0]
        else:
            self.gen_barcodes = product('AGCT',repeat=n)

    def next(self):
        if self.txt_file != '':
            return self.fh.readline().split('\n')[0]
        else:
            try:
                return "".join(self.gen_barcodes.next())
            except:
                pass

    def get_all_bcs(self):
        if self.txt_file != '':
            self.bc_universe = [bc.split('\n')[0] for bc in self.fh]
        else:
            self.bc_universe = ["".join(bc) for bc in self.gen_barcodes]
        return self.bc_universe

class barcodes_generator:
    def __init__(self,target,substring_filter):
        if target != 'inf':
            self.target = int(target)
        else:
            self.target = target
        self.ss_filter = int(substring_filter)
        #self.bc_length = int(bc_length)
        #self.gen_bcs = product('AGCT',repeat=self.bc_length)
        #self.bc_universe = [bc for bc in self.gen_bcs]
        self.ss_dict = {}
        self.bc_out_list = []

    def attempt_addition(self,bc_test,return_bool=False):
        tmp_dict = {}
        for i in range(self.ss_filter,len(bc_test)):
            for j in range(len(bc_test)-i+1):
                #print self.ss_dict
                #print bc_test[j:(j+i)]
                try:
                    self.ss_dict[bc_test[j:(j+i)]]
                    return False
                    #return False
                except:
                    tmp_dict[bc_test[j:(j+i)]] = 0
                    #self.ss_dict[bc_test[j:(j+i)]] = 0
                    continue
        self.ss_dict.update(tmp_dict)
        self.bc_out_list.append(bc_test)
        if self.target == len(self.bc_out_list):
            for bc in self.bc_out_list:
                print bc
            sys.exit()
        return True

def distance_check(bc,bc_out_list,measure='hammond'):
    line = []
    for b in bc_out_list:
        if measure == 'hammond':
            dist = sum([ch1 != ch2 for ch1,ch2 in zip(bc,b)])
        elif measure == 'edit':
            dist = distance(bc,b)
        #print "Hammond",bc,b,sum([ch1 != ch2 for ch1,ch2 in zip(bc,b)])
        #print "edit",bc,b,distance(bc,b)
        line.append(str(dist))
    #sys.exit()
    return ",".join(line)

def distance(a,b):
    """Pure python version to compute the levenshtein distance between a and b.
    The Levenshtein distance includes insertions, deletions, substitutions; 
    unlike the Hamming distance, which is substitutions only.
    ORIGINALLY FROM:  http://hetland.org/coding/python/levenshtein.py
    """
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
    return current[n]

def fails_mismatch(bc1,bc2,n_base):
    hamming_dist = 0
    for i in range(len(bc1)):
        if bc1[i] != bc2[i]:
            hamming_dist += 1
    if hamming_dist <= n_base:
        return [hamming_dist,True]
    else:
        return [hamming_dist,False]


def n_base_mismatch(bc,bcs,n_base):
    failures = ''
    for b in bcs:
        if b != bc:
            out = fails_mismatch(bc,b,n_base)
            if out[1]:
                failures += b+'('+str(out[0])+')|'
    if failures == '':
        return 'Passed,NA'
    return 'Failed,'+failures

def main(argv):
    options = argParser()
    if options.mismatch_test and options.homopolymer_test:
        to_do = '5'
    elif options.mismatch_test:
        to_do = '4'
    elif options.homopolymer_test and options.substring_test:
        to_do = '3'
    elif options.homopolymer_test:
        to_do = '1'
    elif options.substring_test:
        to_do = '2'
    if options.barcode_txt_file != '':
        bcs = open(options.barcode_txt_file,'r').read()
        bcs = bcs.split('\n')[1:-1]
        if len(bcs[0].split(',')) > 1:
            bcs = [l.split(',')[0] for l in bcs]
    else:
        bcs = []
    rev_comp = options.rev_comp
    length = int(options.homopolymer_thresh)
    n_bases = int(options.n_length_barcode)
    filter = int(options.substring_test_filter)
    universe_set = options.universal_barcode_set
    if options.output_file == '':
        outfilename = 'barcodes_'+str(length)+'_homopoly_thresh_'+str(n_bases)+'_base.txt'
    bc_hist = {bc:'' for bc in bcs}
    if to_do == '1':
        print "Barcode,Homopolymer_Test"
        for bc in bcs:
            if len_homopolymer(bc,length):
                print bc+',Failed'
            else:
                print bc+',Passed'
    elif to_do == '2':
        print "Barcode,Substring_Match_Test,Substring_Matched_To"
        print ","+",".join(bcs)
        for bc in bcs:
            print bc+','+substring_test(bc,bcs,n_bases,return_bool=False,return_all=True)
    elif to_do == '3':
        print "Barcode,Homopolymer_Test,Substring_Match_Test,Substring_Matched_To"
        for bc in bcs:
            if len_homopolymer(bc,length):
                print bc+',Failed,NA,NA'
            else:
                print bc+',Passed,'+substring_test(bc,bcs,filter,return_bool=False,rev_comp=rev_comp)
    elif to_do == '4':
        print "Barcode,Mismatch_Test,Mismatch_clashed_with"
        for bc in bcs:
            print bc+','+n_base_mismatch(bc,bcs,(int(options.n_base)*2))
    elif to_do == '5':
        print "Barcode,Homopolymer_Test,Mismatch_Test,Mismatch_clashed_with"
        for bc in bcs:
            if len_homopolymer(bc,length):
                print bc+',Failed,NA,NA'
            else:
                print bc+',Passed,'+n_base_mismatch(bc,bcs,(int(options.n_base)*2))
    else:
        print "unrecognized selection for your second argument... please pass 1,2 or 3"

if __name__ == '__main__':
    main(sys.argv)

