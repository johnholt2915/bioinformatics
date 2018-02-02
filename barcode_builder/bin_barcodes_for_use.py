from random import shuffle
import os
import sys

RB_DICT = {'A':'R','C':'R','G':'B','T':'B'}

def rb_conversion(bc,rb_dict=RB_DICT):
    return "".join([rb_dict[nt] for nt in bc])

def hamming_dist(bc,bcs,n_dist):
    addable = []
    for b in enumerate(bcs):
        if b != bc:
            if len([b[1] for i in range(len(bc)) if bc[i] != b[1][i]]) == n_dist:
                print bc,len([b[1] for i in range(len(bc)) if bc[i] != b[1][i]])
                addable.append(b)
    return addable

def select_seeds(bin,len_alternating_color,barcodes):
    alt_list = ['','']
    for i in range(len_alternating_color):
        if i % 2 == 0:
            alt_list[0] += 'R'
            alt_list[1] += 'B'
        else:
            alt_list[0] += 'B'
            alt_list[1] += 'R'
    output = []
    for b in bin:
        if bin[b][0:len_alternating_color] in alt_list:
            output.append(tuple([b,barcodes[b],len_alternating_color,bin[b]]))
    return output

def main(argv):
    bc_file = argv[1]
    min_alt_len = int(argv[2])
    bins = {4:{},3:{},2:{}}
    rb_bins = {4:{},3:{},2:{}}
    bcs_list = [line.split(',')[0].split('\n')[0] for line in open(bc_file,'r')][1:]
    for i,bc in enumerate(bcs_list):
        bin = len(set(list(bc[0:4])))
        bins[bin][i] = bc
        rb_bc = rb_conversion(bc)
        rb_bins[bin][i] = rb_bc
    out_list = [[],[],[]]
    loop_n = 0
    for i in range(4,1,-1):
        for j in range(len(bcs_list[0]),min_alt_len-1,-1):
            to_add = select_seeds(rb_bins[i],j,bcs_list)
            #print to_add
            #sys.exit()
            for tup in to_add:
                del rb_bins[i][tup[0]]
                out_list[loop_n].append(tup)
        loop_n += 1
    #print sum([len(o) for o in out_list])
    rb_dist = {'R':{i:0 for i in range(len(bcs_list[0]))},'B':{i:0 for i in range(len(bcs_list[0]))}}
    for o in out_list:
        for tup in o:
            for i,col in enumerate(tup[3]):
                rb_dist[col][i] += 1
    #print 'Red:',rb_dist['R']
    #print 'Blue:',rb_dist['B']
    checker_board = []
    cols = {}
    if rb_dist['R'][0] > rb_dist['B'][0]:
        next_col = start_col = 'R'
        cols = {0:'R',1:'B'}
    else:
        next_col = start_col = 'B'
        cols = {0:'B',1:'R'}
    #print "seed_col:",next_col
    found = False
    cb_height = 0
    for i in range(len(bcs_list[0]),min_alt_len-1,-1):
        start_height = i
        for j in range(len(out_list)):
            addable = []
            for l,tup in enumerate(out_list[j]):
                if tup[2] == i and tup[3][0] == next_col:
                    addable.append(tuple([l,tup]))
            #print addable
            if len(addable) != 0:
                shuffle(addable)
                #print addable
                #print "Adding:",addable[0]
                checker_board.append(addable[0][1])
                out_list[j].pop(addable[0][0])
                cb_height = len(checker_board)
                next_col = cols[cb_height % 2]
                found = True
                break
        if found:
            break
    #print ""
    #print "Starting Checker Board:",checker_board

    ## Build the largest Checker Board
    #for i in range(sum([len(o) for o in out_list])-1):
    for i in range(rb_dist[start_col][0]*2):
        #print "Finding Barcode:",i
        found = False
        for j in range(start_height,min_alt_len-1,-1):
            #print "Alternate Length:",j
            for k in range(len(out_list)):
                #print "Out list:",k,"Length:",len(out_list[k])
                addable = []
                for l,tup in enumerate(out_list[k]):
                    if tup[2] == j and tup[3][0] == next_col:
                        addable.append(tuple([l,tup]))
                if len(addable) != 0:
                    shuffle(addable)
                    #print "Adding:",addable[0]
                    checker_board.append(addable[0][1])
                    out_list[k].pop(addable[0][0])
                    cb_height += 1
                    next_col = cols[cb_height%2]
                    found = True
                    break
            if found:
                #print "Found a barcode!"
                #start_height = j
                break
            #else:
                #print "Looking for barcodes with shorter alternating color stretches..."
        #print ""
    #print ""
    #print "Checker Board Height:",len(checker_board)
    print 'Red:',rb_dist['R']
    print 'Blue:',rb_dist['B']
    #for t in checker_board:
    #    print t[0],t[1],t[2],t[3]
    
    ## Write output 
    #for i in range(len(out_list)):
    #    out_file_name = os.path.join(os.path.dirname(bc_file),os.path.basename(bc_file).split('.csv')[0]+'_rb_sorted_bin_'+str(i)+'.csv')
    #    outfile = open(out_file_name,'w')
    #    outfile.write('barcodes'+str("".join([',']*len(out_list[0][0])))+'\n')
    #    for bc in out_list[i]:
    #        outfile.write(str(bc[0])+":"+str(bc[1])+","+",".join(list(bc[1]))+','+str(bc[2])+'\n')
    #    outfile.close()
    #    print "open",out_file_name
    out_file_name = os.path.join(os.path.dirname(bc_file),os.path.basename(bc_file).split('.csv')[0]+'_checkerboard_minAlt'+str(min_alt_len)+'.csv')
    outfile = open(out_file_name,'w')
    outfile.write('barcodes'+str("".join([',']*len(checker_board)))+'\n')
    for bc in checker_board:
        outfile.write(str(bc[0])+":"+str(bc[1])+","+",".join(list(bc[1]))+','+str(bc[2])+'\n')
    outfile.close()
    print "open",out_file_name

if __name__ == '__main__':
    main(sys.argv)


