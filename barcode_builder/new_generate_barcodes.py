import os
import sys
from itertools import product
from random import shuffle
from random import random
import optparse
import math

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

def populate_ss_dict(bcs_universe):
    """This generates a dictionary of all substrings in a universe by providing a universe as input """
    d = {}
    for bc in bcs_universe:
        for j in range(2,len(bc)):
            for k in range(len(bc)-j+1):
                try:
                    d[bc[k:(j+k)]] += 1
                except:
                    d[bc[k:(j+k)]] = 1
    return d

def substring_frequency(bcs_universe,clash_max,score_dict,return_min=False):
    bc_score_dict = {}
    for bc in bcs_universe:
        bc_score_dict[bc] = sum([score_dict[bc[k:(clash_max+k)]] for k in range(len(bc)-clash_max+1)])
    if return_min:
        m_val = min(bc_score_dict.values())
        return [k for k,v in bc_score_dict.items() if v == m_val]
    return bc_score_dict

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

class pseudo_generator:
    def __init__(self,txt_file,n=0,header=False):
        self.txt_file = txt_file
        self.n = n
        self.header = header
        if txt_file != '':
            self.fh = open(txt_file,'w')
            if header:
                self.header = fh.readline().split('\n')[0]
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
        #print "pseudo_generator",len(self.bc_universe)
        return self.bc_universe

def all_four_nts(bc):
    return (len(set(list(bc))) == 4)

REV_COMP = {'A':'T','T':'A','G':'C','C':'G'}
def self_homologous(bc,len_sh,rev_comp=REV_COMP):
    splits = [bc[i:(i+(len_sh/2))] for i in range(len(bc)-(len_sh/2)+1)]
    #print splits
    for i in range(len_sh/2,len(splits)):
        for j in range((len_sh/2)-(i-(len_sh/2))):
            #print splits[j],splits[i+j],reverse_complement(splits[i+j])
            if splits[j] == reverse_complement(splits[i+j]):
                return True
    return False
    #for i in range(len_sh-1):
    #    print i
    #    for j in range(int(math.ceil(len(bc)/((len_sh/2)+i)))):
    #        print j
    #        print splitter[j]
    #        print splitter[j+(len_sh/2)+i]
    #        print reverse_complement(splitter[j+(len_sh/2)+i])
    #        print ""
    #        if splitter[j] == reverse_complement(splitter[j+(len_sh/2)+i]):
    #            print splitter[j],splitter[j+(len_sh/2)+i],reverse_complement(splitter[j+(len_sh/2)+i])
    #            return True
    #return False

def return_barcodes(n_base=4,len_barcodes=8,all_base=False,homopoly_len=2):
    pseudo_gen = pseudo_generator('',n=len_barcodes)
    if all_base:
        #tmp = [bc for bc in pseudo_gen.get_all_bcs() if len(set(list(bc[0:4]))) == n_base]
        #print "len(set(list(bc[0:4]))) ==",n_base,len(tmp)
        #tmp = [bc for bc in tmp if not len_homopolymer(bc,homopoly_len)]
        #print "not len_homopolymer(bc,",str(homopoly_len)+')',len(tmp)
        #tmp = [bc for bc in tmp if not self_homologous(bc,6)]
        #print "not self_homologous(bc,6)",len(tmp)
        #tmp = [bc for bc in tmp if all_four_nts(bc)]
        #print "all_four_nts(bc)",len(tmp)
        #return tmp
        return [bc for bc in pseudo_gen.get_all_bcs() if len(set(list(bc[0:4]))) == n_base and not len_homopolymer(bc,homopoly_len) and not self_homologous(bc,6) and all_four_nts(bc)]
    else:
        #tmp = [bc for bc in pseudo_gen.get_all_bcs() if len(set(list(bc[0:4]))) == n_base]
        #print "len(set(list(bc[0:4]))) ==",n_base,len(tmp)
        #tmp = [bc for bc in tmp if not len_homopolymer(bc,homopoly_len)]
        #print "not len_homopolymer(bc,",str(homopoly_len)+")",len(tmp)
        #tmp = [bc for bc in tmp if not self_homologous(bc,6)]
        #print "not self_homologous(bc,6)",len(tmp)
        #tmp = [bc for bc in tmp if all_four_nts(bc)]
        #print "all_four_nts(bc)",len(tmp)
        #return tmp
        return [bc for bc in pseudo_gen.get_all_bcs() if len(set(list(bc[0:4]))) == n_base and not len_homopolymer(bc,homopoly_len) and not self_homologous(bc,6)]

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

DIST_VALS = {}

def perform_distance_measure(bc_universe):
    seen = {}
    for bc1 in bc_universe:
        print bc1
        DIST_VALS[bc1] = {}
        for bc2 in bc_universe:
            try:
                seen[bc2+'_'+bc1]
            except:
                DIST_VALS[bc1][bc2] = distance(bc1,bc2)

def check_distance(bc,bc_list):
    dist_list = []
    for b in bc_list:
        dist_list.append(distance(bc,b))
    return float(sum(dist_list))/len(dist_list)

REV_COMP = {'A':'T','T':'A','G':'C','C':'G'}
def reverse_complement(bc,rc=REV_COMP):
    return "".join([rc[nt] for nt in bc][::-1])

def adjust_universe(bc,bcs_universe,ss_dict,unv_ss_dict,clash,rev=False):
    bcs_universe_out = []
    cmbs = []
    if rev:
        bc = reverse_complement(bc)
    ## remove all substrings of the barcode provided (bc) from unv_ss_dict
    ## There are no longer those substrings represented in the universe
    subtractions = 0
    for j in range(2,len(bc)):
        for k in range(len(bc)-j+1):
            ss_dict[bc[k:(j+k)]] = 0
            cmbs.append(bc[k:(j+k)])
    ## if any barcodes in the universe clash with at least one substring length 2 - 7, remove the barcodes
    ## from the universe and thus their substrings
    failures = 0
    for b in bcs_universe:
        failed = False
        #for ss in [c for c in cmbs if len(c) == clash]:
        for ss in cmbs:
            if ss in b:
                ## subtract all substrings from the universe substring dictionary
                for tmp_ss in [b[m:(l+m)] for l in range(2,len(b)) for m in range(len(b)-l+1)]:
                    subtractions += 1
                    unv_ss_dict[tmp_ss] -= 1
                failed = True
                failures += 1
                break
        if not failed:
            bcs_universe_out.append(b)
    #print "failures:",failures
    #print "subtractions:",subtractions
    return bcs_universe_out,ss_dict,unv_ss_dict

def fails_mismatch(bc1,bc2,n_base):
    hamming_dist = 0
    for r in range(len(bc1)):
        if bc1[r] != bc2[r]:
            hamming_dist += 1
    if hamming_dist < n_base:
        return True
    return False


def n_base_mismatch(bc,bcs,n_base):
    failures = ''
    for b in bcs:
        if b != bc:
            if fails_mismatch(bc,b,n_base):
                failures += b+'|'
    if failures == '':
        return False
    return True

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-b', '--base-mismatch', dest='base_mismatch',default=1)
    parser.add_option( '-s', '--substring-match-threshold', dest='substring_match_threshold',default=5)
    parser.add_option( '-a', '--all-bases-represented', dest="all_bases_represented",action='store_true',default=False)
    parser.add_option( '-f', '--first-four-only',dest='first_four_only',action='store_true',default=False)
    parser.add_option( '-e', '--epochs',dest='epochs',default=10)
    parser.add_option( '-o', '--output-directory', dest='outdir',default='.')
    parser.add_option( '-c', '--barcode-seed-set',dest='barcode_seed_set_file',default='')
    parser.add_option( '-p', '--max-homopolymer-len',dest='max_homopolymer_len',default=2)
    parser.add_option( '-n', '--n-base-barcode',dest='n_base_barcode',default=8)
    
    ( options, args ) = parser.parse_args()
    print "###### Barcode Generator Input Parameters #######"
    for o in options.__dict__:
        print o + "\t" + str(options.__dict__[o])
    return options

def main(argv):
    options = argParser()
    print "###### Barcode Generator Main Function #######"
    # adding plus 1 because my algorithm actually checks for 1 base mm less than the value passed as input.
    # my algorithm will accept the two 5-base barcodes below, even though it is possible for them to clash with 1 miscall each
    # ACGTG ==> AGGTG (transversion from C to G at 2nd index)
    # AGCTG ==> AGGTG (transversion from C to G at 3rd index)
    base_mismatch = 2*(int(options.base_mismatch))+1
    ss_match_fail = int(options.substring_match_threshold)
    all_bases_represented = options.all_bases_represented
    first_4_only = options.first_four_only
    hp_len = int(options.max_homopolymer_len)
    n_base_barcode = int(options.n_base_barcode)
    if first_4_only:# only attempt to add barcodes where first 4 bases have all 4 nts repsented at least 1 time.
        outermost_loop_list = [4]
    else:
        outermost_loop_list = range(4,1,-1)
    epochs = int(options.epochs)
    outdir = options.outdir
    ## lists and dictionaries for running the algorithm
    bcs_universe = []
    rev_bcs_universe = []
    ss_dict = {}
    rev_ss_dict = {}
    bcs_list = []
    rev_bcs_list = []
    rev_bcs_distance_check_list = []
    if options.barcode_seed_set_file != '': # this means that the user passed a list of barcodes as input in the form of a csv file with a header, each line is a new barcode
        bcs_list = [bc.split(',')[0] for bc in open(options.barcode_seed_set_file).read().split('\n')[1:-1]]
        rev_bcs_list = [bc for bc in bcs_list]
        for bc_test in bcs_list:
            for j in range(2,len(bc_test)):
                for k in range(len(bc_test)-j+1):
                    ss_dict[bc_test[k:(j+k)]] = 0
                    rev_ss_dict[bc_test[k:(j+k)]] = 0
            rev_bcs_distance_check_list.append(reverse_complement(bc_test))
    last_dir = os.path.join(outdir,'last_dir')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isdir(last_dir):
        os.mkdir(last_dir)
    ## Main Algorithm
    for n_base in outermost_loop_list:
        #print "bcs_universe with all",str(n_base),"bases in first four..."
        ## generate a universe of barcodes based on n_base
        bcs_universe = return_barcodes(n_base,len_barcodes=n_base_barcode,all_base=all_bases_represented,homopoly_len=hp_len)
        #print len(bcs_universe)
        #continue
        shuffle(bcs_universe)
        print "n_barcodes when n_base =",str(n_base)+":",len(bcs_universe)

        ## Setup
        if len(bcs_list) == 0:# if the barcode list is empty, we need to add one to the list and adjust universe
            print "STARTING BARCODE:",bcs_universe[0]
            bcs_list = [bcs_universe[0]]
            rev_bcs_list = [bcs_universe[0]]
            rev_bcs_distance_check_list = [reverse_complement(bcs_universe[0])]
            for j in range(2,len(bcs_list[0])):
                for k in range(len(bcs_list[0])-j+1):
                    ss_dict[bcs_list[0][k:(j+k)]] = 0
                    rev_ss_dict[bcs_list[0][k:(j+k)]] = 0
            bcs_universe.remove(bcs_list[0])
        rev_bcs_universe = [reverse_complement(bc) for bc in bcs_universe]
        for_rev= [False,True]
        last = []
        last_filt = 0

        ## Epochs
        for epoch in range(epochs):
            print "Epoch:",epoch+1
            print "Reverse barcode universe size:",len(rev_bcs_universe)
            addable_lists = [[],[]]
            addable_list_dists = [[],[]]
            ## Forward attempt 2,3,4,5 size filtering to find at most 1 barcode to add
            if for_rev[0]:
                for filt in range(2,ss_match_fail+1):
                    print "filter:",filt
                    throw_away_bcs= []
                    for bc_test in bcs_universe:
                        add_bc = True
                        #for ss_bc in [bc_test[k:(k+j)] for j in range(filt,len(bc_test)) for k in range(len(bc_test)-j+1)]:
                        #for ss_bc in [bc_test[j:(filt+j)] for j in range(len(bc_test)-filt+1)]:
                        for ss_bc in [bc_test[k:(k+j)] for j in range(filt,ss_match_fail+1) for k in range(len(bc_test)-j+1)]:
                            try:
                                ss_dict[ss_bc]
                                add_bc = False
                                break
                            except:
                                continue
                        if add_bc:
                            if not n_base_mismatch(bc_test,bcs_list,base_mismatch):
                                d = check_distance(bc_test,bcs_list)
                                addable_lists[0].append(bc_test)
                                addable_list_dists[0].append(d)
                            elif filt >= ss_match_fail:
                                print "removing",bc_test
                                throw_away_bcs.append(bc_test)
                        elif filt == ss_match_fail:
                            print "removing",bc_test
                            throw_away_bcs.append(bc_test)
                    print "Number of addable barcodes:",len(addable_lists[0])
                    if len(addable_lists[0]) > 0:
                        m = max(addable_list_dists[0])
                        addable_list = [addable_lists[0][i] for i,j in enumerate(addable_list_dists[0]) if j == m]
                        shuffle(addable_list)
                        bcs_list.append(addable_list[0])
                        print "Added",addable_list[0],"to Reverse List"
                        # after adding the barcode to the list of barcodes, then add the substrings to the dictionary
                        for j in range(2,len(bcs_list[-1])):
                            for k in range(len(bcs_list[-1])-j+1):
                                ss_dict[bcs_list[-1][k:(k+j)]] = 0
                        #print "New length of rev_ss_dict:",len(rev_ss_dict)
                        ## remove the added barcode.
                        bcs_universe.remove(bcs_list[-1])
                        if filt == ss_match_fail: # and all others which clash 5 length substring from the universe
                            addable_lists[0].remove(bcs_list[-1])
                            for code in addable_lists[0]:
                                for j in range(ss_match_fail-1):
                                    if code[j:(j+ss_match_fail)] in bcs_list[-1]:
                                        throw_away_bcs.append(code)
                                        break
                            for b in throw_away_bcs:
                                bcs_universe.remove(b)
                        break
                    else:
                        if filt == ss_match_fail:
                            print "NO MORE REVERSE TO BE FOUND!!!"
                            for_rev[0] = False
                        continue

            ## Reverse attempt
            if for_rev[1]:
                for filt in range(2,ss_match_fail+1):
                    print "reverse filter:",filt
                    throw_away_bcs= []
                    for bc_test in rev_bcs_universe:
                        add_bc = True
                        #for ss_bc in [bc_test[k:(k+j)] for j in range(filt,len(bc_test)) for k in range(len(bc_test)-j+1)]:
                        #for ss_bc in [bc_test[j:(j+filt)] for j in range(len(bc_test)-filt+1)]:
                        for ss_bc in [bc_test[k:(k+j)] for j in range(filt,ss_match_fail+1) for k in range(len(bc_test)-j+1)]:
                            try:
                                rev_ss_dict[ss_bc]
                                add_bc = False
                                break
                            except:
                                continue
                        if add_bc:
                            rc_bc_test = reverse_complement(bc_test)
                            if not n_base_mismatch(rc_bc_test,rev_bcs_list,base_mismatch) and rc_bc_test not in rev_bcs_list:
                                d = check_distance(bc_test,rev_bcs_distance_check_list)
                                addable_lists[1].append(rc_bc_test)
                                addable_list_dists[1].append(d)
                            elif filt == ss_match_fail:
                                #print "removing",bc_test
                                throw_away_bcs.append(bc_test)
                        elif filt == ss_match_fail:
                            #print "removing",bc_test
                            throw_away_bcs.append(bc_test)
                    print "Number of addable barcodes:",len(addable_lists[1])
                    if len(addable_lists[1]) > 0:
                        m = max(addable_list_dists[1])
                        addable_list = [addable_lists[1][i] for i,j in enumerate(addable_list_dists[1]) if j == m]
                        shuffle(addable_list)
                        print addable_list[0],"in rev_bcs_list?",addable_list[0] in rev_bcs_list
                        rev_bcs_list.append(addable_list[0])
                        print "Added",addable_list[0],"to Reverse List"
                        rev_bcs_distance_check_list.append(reverse_complement(addable_list[0]))
                        # after adding the barcode to the list of barcodes, then add the substrings to the dictionary
                        for j in range(2,len(rev_bcs_list[-1])):
                            for k in range(len(rev_bcs_list[-1])-j+1):
                                rev_ss_dict[rev_bcs_list[-1][k:(k+j)]] = 0
                        #print "New length of rev_ss_dict:",len(rev_ss_dict)
                        ## remove the added barcode. 
                        rev_bcs_universe.remove(reverse_complement(rev_bcs_list[-1]))
                        if filt == ss_match_fail:# and all others which clash 5 length substring from the universe
                            addable_lists[1].remove(rev_bcs_list[-1])
                            for code in addable_lists[1]:
                                rc_code = reverse_complement(code)
                                for j in range(len(bcs_list[0])-ss_match_fail+1):
                                    if rc_code[j:(j+ss_match_fail)] in rev_bcs_list[-1]:
                                        throw_away_bcs.append(rc_code)
                                        break
                            for b in throw_away_bcs:
                                rev_bcs_universe.remove(b)
                        break
                    else:
                        if filt == ss_match_fail:
                            print "NO MORE REVERSE TO BE FOUND!!!"
                            for_rev[1] = False
                        continue
            if not any(for_rev):
                break
            print ""
        print "couldn't find another barcode from universe... checking another universe..."
    ## make output file name
    #bc_[fr]_abr_ff_bmm#_ssmt#_ep#_seedset.csv
    abr = ""
    ff = ""
    seedset = ""
    if all_bases_represented:
        abr = "abr"
    if first_4_only:
        ff = "ff"
    if options.barcode_seed_set_file != '':
        seedset = "seedset"
    out_file_name = os.path.join(outdir,"_".join([rev_bcs_list[0],"reverse",abr,ff,"bmm"+str(options.base_mismatch),"ssmt"+str(ss_match_fail),"ep"+str(epochs),"hp"+str(options.max_homopolymer_len),seedset])+".csv")
    print "less",out_file_name
    out = open(out_file_name,'w')
    if seedset != "":
        out.write('barcodes,'+os.path.basename(out_file_name)+','+options.barcode_seed_set_file+"".join([',']*(len(rev_bcs_list[0])-2))+'\n')
    else:
        out.write('barcodes,'+os.path.basename(out_file_name)+"".join([',']*(len(rev_bcs_list[0])-1))+'\n')
    for i,bc in enumerate(rev_bcs_list):
        #print bc,bcs_list[i]
        out.write(bc+','+",".join(list(bc))+'\n')
    out.close()
    out_file_name = os.path.join(outdir,"_".join([bcs_list[0],"forward",abr,ff,"bmm"+str(options.base_mismatch),"ssmt"+str(ss_match_fail),"ep"+str(epochs),"hp"+str(options.max_homopolymer_len),seedset])+".csv")
    out = open(out_file_name,'w')
    if seedset != "":
        out.write('barcodes,'+os.path.basename(out_file_name)+','+options.barcode_seed_set_file+"".join([',']*(len(bcs_list[0])-2))+'\n')
    else:
        out.write('barcodes,'+os.path.basename(out_file_name)+"".join([',']*(len(bcs_list[0])-1))+'\n')
    for bc in bcs_list:
        out.write(bc+','+",".join(list(bc))+'\n')
    out.close()

if __name__ == '__main__':
    main(sys.argv)


