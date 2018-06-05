import os

dir = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/output_dir_revised'
os.chdir(dir)
dirs = os.listdir('.')
d_dict = {}
for d in dirs:
    if len(d.split('_kegg_')) == 2:
        k = d.split('_kegg_')[0]
    elif len(d.split('_reactome_')) == 2:
        k = d.split('_reactome_')[0]
    else:
        k = d     
    try:
        d_dict[k].append(d)
    except:
        d_dict[k] = [d]

# filter the folders with nothing in them ... 
hist = {}
for d in d_dict:
    try:
        hist[len(d_dict[d])] +=1
    except:
        hist[len(d_dict[d])] = 1
#>>> hist
#{1: 889, 2: 901, 3: 85}

for d in d_dict:
    d_dict[d] = [v for v in d_dict[d] if 'kegg' in v or 'react' in v]

hist2 = {}
for d in d_dict:
    try:
        hist2[len(d_dict[d])] +=1
    except:
        hist2[len(d_dict[d])] = 1
#>>> hist2
#{0: 889, 2: 986}


### find the genepairs that were most common pan-cancer
canc_dirs = [os.path.join('..',d) for d in os.listdir('../') if '.new' in d]

# make a dictionary where keys are parameters and values are the files with lists of genepairs
param_dict_gp = {}
for c in canc_dirs:
    for f in os.listdir(c):
        param = f.split('.new_')[1].split('.tsv')[0]
        try:
            param_dict_gp[param] += [os.path.join(c,f)]
        except:
            param_dict_gp[param] = [os.path.join(c,f)]

param_hist_gp = {i:0 for i in range(1,16)}
for p in param_dict_gp:
    param_hist_gp[len(param_dict_gp[p])] += 1

# to confirm that all tests were performed (all should be length 15...)
for i in sorted(param_hist_gp.keys()):
    print i,param_hist_gp[i]

def getGenePairs(f):
    gp_set = {}
    with open(f,'r') as fh:
        for line in fh:
            gp = "|".join(sorted(line.strip().split('\t')[2:]))
            gp_set[gp] = 0
    return gp_set.keys()

all_params_hist = {p:{} for p in param_dict_gp}
for p in param_dict_gp:
    print p
    for f in param_dict_gp[p]:
        c = os.path.basename(f).split('_pairs_')[1].split('_ndel_')[0]
        for gp in getGenePairs(f):
            try:
                all_params_hist[p][gp] += [c]
            except:
                all_params_hist[p][gp] = [c]

## SANITY CHECK!!!
missed = {}
for p in all_params_hist:
    print p
    for gp in all_params_hist[p]:
        if len(gp.split('|')) != 2:
            print gp
            try:
                missed[p].append(gp)
            except:
                missed[p] = [gp]


output_stats_dir = '../stats_out_revised'
if not os.path.isdir(output_stats_dir):
    os.mkdir(output_stats_dir)

def writeGenePairs(p,p_gp_dict):
    of = open(os.path.join(output_stats_dir,'gene_pairs_'+p+'.csv'),'w')
    for gp in p_gp_dict:
        line = "".join(gp.split('\n'))+'\t'+'\t'.join(p_gp_dict[gp])+'\n'
        of.write(line)
    of.close()
    print "wrote %s" % (os.path.join(output_stats_dir,'gene_pairs_'+p+'.csv'))

for p in all_params_hist:
    writeGenePairs(p,all_params_hist[p])

of = open(os.path.join('..','number_of_pan_cancer_genepairs_per_parameter_set.csv'),'w')
for i in all_params_hist:
    of.write(str(i)+'\t'+str(len(all_params_hist[i]))+'\n')

of.close()

# pickle the results thusfar
#import pickle
#with open('../checkpoint_1.pkl','wb') as fh:
#    pickle.dump([dir,d_dict,canc_dirs,param_dict_gp,all_params_hist],fh)
#
#with open('../checkpoint_1.pkl','rb') as fh:
#    dir,d_dict,canc_dirs,param_dict_gp,all_params_hist_gp = pickle.load(fh)


### find the pathwys pairs most common pan-cancer
param_dict_pw = {}
for d in d_dict:
    p = d.split('.new_')[1]
    try:
        param_dict_pw[p]['kegg'] += [t for t in d_dict[d] if 'kegg' in t]
        param_dict_pw[p]['reactome'] += [t for t in d_dict[d] if 'reactome' in t]
    except:
        param_dict_pw[p] = {'kegg':[],'reactome':[]}
        param_dict_pw[p]['kegg'] += [t for t in d_dict[d] if 'kegg' in t]
        param_dict_pw[p]['reactome'] += [t for t in d_dict[d] if 'reactome' in t]

def getPathways(tab3):
    pw_set = {}
    with open(tab3,'r') as fh:
        for line in fh:
            sep = line.split('\t')
            pw1 = sep[2]
            pw2 = sep[14]
            pw_key = "|".join(sorted([pw1,pw2]))
            pw_set[pw_key] = 0
    return pw_set.keys()


all_params_kegg_pws = {p:{} for p in param_dict_pw}
db = 'kegg'
for p in param_dict_pw:
    print p
    for d in param_dict_pw[p][db]:
        c = d.split('_pairs_')[1].split('_ndel_')[0]
        for pw in getPathways(d+'/table_3.txt'):
            try:
                all_params_kegg_pws[p][pw] += [c]
            except:
                all_params_kegg_pws[p][pw] = [c]

all_params_reactome_pws = {p:{} for p in param_dict_pw}
db = 'reactome'
for p in param_dict_pw:
    print p
    for d in param_dict_pw[p][db]:
        c = d.split('_pairs_')[1].split('_ndel_')[0]
        for pw in getPathways(d+'/table_3.txt'):
            try:
                all_params_reactome_pws[p][pw] += [c]
            except:
                all_params_reactome_pws[p][pw] = [c]

def writePathways(p,p_pw_dict,db):
    of = open(os.path.join(output_stats_dir,'pathways_'+p+'_'+db+'.csv'),'w')
    for pw in p_pw_dict:
        line = "".join(pw.split('\n'))+'\t'+'\t'.join(p_pw_dict[pw])+'\n'
        of.write(line)
    of.close()
    print "wrote %s" % (os.path.join(output_stats_dir,'pathways_'+p+'_'+db+'.csv'))

for p in all_params_kegg_pws:
    writePathways(p,all_params_kegg_pws[p],'kegg')
    writePathways(p,all_params_reactome_pws[p],'reactome')

of = open(os.path.join('..','number_of_pan_cancer_kegg_pathways_per_parameter_set.csv'),'w')
for i in all_params_kegg_pws:
    of.write(str(i)+'\t'+str(len(all_params_kegg_pws[i]))+'\n')

of.close()
of = open(os.path.join('..','number_of_pan_cancer_reactome_pathways_per_parameter_set.csv'),'w')
for i in all_params_reactome_pws:
    of.write(str(i)+'\t'+str(len(all_params_reactome_pws[i]))+'\n')

of.close()

        

### find commonly occurring gene pairs across all cancers for each parameter set
dir = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/stats_out_revised'
os.chdir(dir)

def split_files():
    other = []
    genes = []
    kegg = []
    reactome = []
    for i in os.listdir('.'):
        if 'gene_pair' in i:
            genes += [i]
        elif 'kegg' in i:
            kegg += [i]
        elif 'reactome' in i:
            reactome += [i]
        else:
            other += [i]
    return genes,kegg,reactome,other

gene_files,kegg_files,reactome_files,other_files = split_files()

files_dict = {}
for f in [g for g in gene_files if 'gene_pairs' in g]:
    key = f.split('.csv')[0].split('_pairs_')[1]
    files_dict[key] = {}
    print f
    with open(f,'r') as fh:
        for line in fh:
            sep = line.strip().split('\t')
            files_dict[key][sep[0]] = len(sep[1:])

## gather all gene pairs which were found in a number (n) cancers
param_hists_counts = {}
for p in files_dict:
    param_hists_counts[p] = {}
    print p
    for gp in files_dict[p]: 
        try:    
            param_hists_counts[p][files_dict[p][gp]] += 1
        except:
            param_hists_counts[p][files_dict[p][gp]] = 1

param_hists_gps = {}
for p in files_dict:
    param_hists_gps[p] = {}
    print p
    for gp in files_dict[p]: 
        try:                             
            param_hists_gps[p][files_dict[p][gp]] += [gp]
        except:
            param_hists_gps[p][files_dict[p][gp]] = [gp]

# create tables where columns are the number of cancers a gene pair was found in
for p in param_hists_gps:
    of = open('gene_pair_histograms_'+p+'.csv','w')
    try:
        longest = len(param_hists_gps[p][min(param_hists_gps[p])])
    except ValueError as e:
        print "skipping %s" % (p)
        of.close()
        continue
    header = "\t".join([str(i) for i in param_hists_gps[p]])+'\t'+"\n"
    of.write(header)
    counts_line = "\t".join([str(param_hists_counts[p][k]) for k in param_hists_gps[p]])+'\t'+'\n'
    of.write(counts_line)
    for i in range(longest):
        line = ""
        for j in sorted(param_hists_gps[p].keys()):
            try:
                line += str(param_hists_gps[p][j][i]) +"\t"
            except:
                line += "\t"
        of.write(line+"\n")
    of.close()
    print "wrote %s" % ('gene_pair_histograms_'+p+'.csv')

# what genepairs are found most frequently accross parameters
all_gps = {}
missing_pair = {}
for p in files_dict:
    print p
    for gp in files_dict[p]:
        if len(gp.split('|')) != 2:
            try:
                missing_pair[gp] += [p]
            except:
                missing_pair[gp] = [p]
        else:
            try:
                all_gps[gp] +=1
            except:
                all_gps[gp] = 1

hist_counts = {}
for gp in all_gps:
    try: 
        hist_counts[all_gps[gp]] += 1
    except:
        hist_counts[all_gps[gp]] = 1


hist = {}
for gp in all_gps:
    try:  
        hist[all_gps[gp]] += [gp]
    except:
        hist[all_gps[gp]] = [gp]

longest = max([len(v) for v in hist.values()])
of = open('../pan_parameter_histogram.csv','w')
header = "\t".join([str(i) for i in sorted(hist.keys())])+'\n'
of.write(header)
counts_line = "\t".join([str(hist_counts[k]) for k in hist])+'\t'+'\n'
of.write(counts_line)
for i in range(longest):
    line = ""
    for j in sorted(hist.keys()):
        try:
            line += str(hist[j][i]) + '\t'
        except:
            line += '\t'
    of.write(line+'\n')

of.close()

### Repeat this procedure for kegg and reactome pathway pairs
kegg_pathway_files_dict = {}
for f in [v for v in kegg_files if '_kegg.csv' in v]:
    key = f.split('_kegg.csv')[0].split('pathways_')[1]
    kegg_pathway_files_dict[key] = {}
    print f
    with open(f,'r') as fh:
        for line in fh:
            sep = line.strip().split('\t')
            kegg_pathway_files_dict[key][sep[0]] = len(sep[1:])

kegg_param_hists_counts = {}
for p in kegg_pathway_files_dict:
    kegg_param_hists_counts[p] = {}
    print p
    for pw in kegg_pathway_files_dict[p]: 
        try:   
            kegg_param_hists_counts[p][kegg_pathway_files_dict[p][pw]] += 1
        except:
            kegg_param_hists_counts[p][kegg_pathway_files_dict[p][pw]] = 1

kegg_param_hists_pws = {}
for p in kegg_pathway_files_dict:
    kegg_param_hists_pws[p] = {}
    print p
    for pw in kegg_pathway_files_dict[p]: 
        try:   
            kegg_param_hists_pws[p][kegg_pathway_files_dict[p][pw]] += [pw]
        except:
            kegg_param_hists_pws[p][kegg_pathway_files_dict[p][pw]] = [pw]

for p in kegg_param_hists_pws:
    of = open('kegg_pathway_histograms_'+p+'.csv','w')
    try:
        longest = len(kegg_param_hists_pws[p][min(kegg_param_hists_pws[p])])
    except ValueError as e:
        print "skipping %s" % (p)
        of.close()
        continue
    header = "\t".join([str(i) for i in kegg_param_hists_pws[p]])+'\t'+"\n"
    of.write(header)
    counts_line = "\t".join([str(kegg_param_hists_counts[p][k]) for k in kegg_param_hists_pws[p]])+'\t'+'\n'
    of.write(counts_line)
    for i in range(longest):
        line = ""
        for j in sorted(kegg_param_hists_pws[p].keys()):
            try:
                line += str(kegg_param_hists_pws[p][j][i]) +"\t"
            except:
                line += "\t"
        of.write(line+"\n")
    of.close()
    print "wrote %s" % ('kegg_pathway_histograms_'+p+'.csv')


all_kegg_pws = {}
missing_pair = {}
for p in kegg_pathway_files_dict:
    print p
    for pw in kegg_pathway_files_dict[p]:
        if len(pw.split('|')) != 2:
            try:
                missing_pair[pw] += [p]
            except:
                missing_pair[pw] = [p]
        else:
            try:
                all_kegg_pws[pw] +=1
            except:
                all_kegg_pws[pw] = 1

hist_counts = {}
for pw in all_kegg_pws:
    try:
        hist_counts[all_kegg_pws[pw]] += 1
    except:
        hist_counts[all_kegg_pws[pw]] = 1


hist = {}
for pw in all_kegg_pws:
    try:
        hist[all_kegg_pws[pw]] += [pw]
    except:
        hist[all_kegg_pws[pw]] = [pw]

longest = max([len(v) for v in hist.values()])
of = open('../pan_parameter_kegg_pathway_histogram.csv','w')
header = "\t".join([str(i) for i in sorted(hist.keys())])+'\n'
of.write(header)
counts_line = "\t".join([str(hist_counts[k]) for k in hist])+'\t'+'\n'
of.write(counts_line)
for i in range(longest):
    line = ""
    for j in sorted(hist.keys()):
        try:
            line += str(hist[j][i]) + '\t'
        except:
            line += '\t'
    of.write(line+'\n')

of.close()


##### REACTOME 
reactome_pathway_files_dict = {}
for f in [v for v in reactome_files if '_reactome.csv' in v]:
    key = f.split('_reactome.csv')[0].split('pathways_')[1]
    reactome_pathway_files_dict[key] = {}
    print f
    with open(f,'r') as fh:
        for line in fh:
            sep = line.strip().split('\t')
            reactome_pathway_files_dict[key][sep[0]] = len(sep[1:])

reactome_param_hists_counts = {}
for p in reactome_pathway_files_dict:
    reactome_param_hists_counts[p] = {}
    print p
    for pw in reactome_pathway_files_dict[p]:
        try:
            reactome_param_hists_counts[p][reactome_pathway_files_dict[p][pw]] += 1
        except:
            reactome_param_hists_counts[p][reactome_pathway_files_dict[p][pw]] = 1

reactome_param_hists_pws = {}
for p in reactome_pathway_files_dict:
    reactome_param_hists_pws[p] = {}
    print p
    for pw in reactome_pathway_files_dict[p]:
        try:
            reactome_param_hists_pws[p][reactome_pathway_files_dict[p][pw]] += [pw]
        except:
            reactome_param_hists_pws[p][reactome_pathway_files_dict[p][pw]] = [pw]

for p in reactome_param_hists_pws:
    of = open('reactome_pathway_histograms_'+p+'.csv','w')
    try:
        longest = len(reactome_param_hists_pws[p][min(reactome_param_hists_pws[p])])
    except ValueError as e:
        print "skipping %s" % (p)
        of.close()
        continue
    header = "\t".join([str(i) for i in reactome_param_hists_pws[p]])+'\t'+"\n"
    of.write(header)
    counts_line = "\t".join([str(reactome_param_hists_counts[p][k]) for k in reactome_param_hists_pws[p]])+'\t'+'\n'
    of.write(counts_line)
    for i in range(longest):
        line = ""
        for j in sorted(reactome_param_hists_pws[p].keys()):
            try:
                line += str(reactome_param_hists_pws[p][j][i]) +"\t"
            except:
                line += "\t"
        of.write(line+"\n")
    of.close()
    print "wrote %s" % ('reactome_pathway_histograms_'+p+'.csv')


all_reactome_pws = {}
missing_pair = {}
for p in reactome_pathway_files_dict:
    print p
    for pw in reactome_pathway_files_dict[p]:
        if len(pw.split('|')) != 2:
            try:
                missing_pair[p] += [pw]
            except:
                missing_pair[p] = [pw]
        else:
            try:
                all_reactome_pws[pw] +=1
            except:
                all_reactome_pws[pw] = 1

hist_counts = {}
for pw in all_reactome_pws:
    try:
        hist_counts[all_reactome_pws[pw]] += 1
    except:
        hist_counts[all_reactome_pws[pw]] = 1


hist = {}
for pw in all_reactome_pws:
    try:
        hist[all_reactome_pws[pw]] += [pw]
    except:
        hist[all_reactome_pws[pw]] = [pw]

longest = max([len(v) for v in hist.values()])
of = open('../pan_parameter_reactome_pathway_histogram.csv','w')
header = "\t".join([str(i) for i in sorted(hist.keys())])+'\n'
of.write(header)
counts_line = "\t".join([str(hist_counts[k]) for k in hist])+'\t'+'\n'
of.write(counts_line)
for i in range(longest):
    line = ""
    for j in sorted(hist.keys()):
        try:
            line += str(hist[j][i]) + '\t'
        except:
            line += '\t'
    of.write(line+'\n')

of.close()


##### select gene pairs and pathway pairs for survival analysis
import os
dir = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/stats_out_revised'
os.chdir(dir)

def read_histogram_file(f):
    with open(f,'r') as fh:
        print "reading",f
        head = fh.readline()
        if head[-2:] != '\t\n':
            head = head.strip()+'\t\n'
        head_dict = {}
        for i,v in enumerate(head.split('\t')):
            try:
                head_dict[i] = int(v)
            except:
                head_dict[i] = -1
        bins = {i:[] for i in range(min(head_dict.values()),max(head_dict.values())+1)}
        counts = fh.readline()
        counts_dict = {head_dict[i]:int(v) for i,v in enumerate(counts.strip().split('\t'))}
        for line in fh:
            for i,v in enumerate(line.split('\t')):
                bins[head_dict[i]] += [v]
    return bins


pan_files = [f for f in os.listdir('..') if 'pan' in f]
gps_parameter_histogram = read_histogram_file('../pan_parameter_histogram.csv')

## select the parameter combo with the most gene pairs
## (theoretically this will be super set of all other parameter combos)
gps_per_param = {}
with open('../number_of_pan_cancer_genepairs_per_parameter_set.csv','r') as fh:
    for line in fh:
        sep = line.strip().split('\t')
        gps_per_param[sep[0]] = int(sep[1])

max_ind = ""
max_val = 0
for i in gps_per_param:
    if gps_per_param[i] > max_val:
        max_ind = i
        max_val = gps_per_param[i]

max_file= [f for f in gene_files if max_ind in f and 'histogram' in f][0]

max_gps_parameter_histogram = read_histogram_file(max_file)
all_gps_param = {}
xy = {'x':gps_parameter_histogram,'y':max_gps_parameter_histogram}
for i in xy:
    print "determinating",i,"values"
    for j in xy[i]:
        for gp in xy[i][j]:
            try:
                all_gps_param[gp][i] = j
            except:
                all_gps_param[gp] = {i:j}

## write the coordinates for all genepairs with both x,y coordinates
of = open('../'+max_ind+'_coordinates.csv','w')
head = 'genepair\tx\ty\n'
of.write(head)
count = 0
for i in all_gps_param:
    if len(all_gps_param[i]) != 2:
        continue
    of.write(i+'\t'+str(all_gps_param[i]['x'])+'\t'+str(all_gps_param[i]['y'])+"\n")
    count +=1

of.close()

## Repeat the coordinate identification for pathways
kegg_parameter_histogram = read_histogram_file('../pan_parameter_kegg_pathway_histogram.csv')

pws_per_param = {}
with open('../number_of_pan_cancer_kegg_pathways_per_parameter_set.csv','r') as fh:
    for line in fh:
        sep = line.strip().split('\t')
        pws_per_param[sep[0]] = int(sep[1])


def split_files():
    other = []
    genes = []
    kegg = []
    reactome = []
    for i in os.listdir('.'):
        if 'gene_pair' in i:
            genes += [i]
        elif 'kegg' in i:
            kegg += [i]
        elif 'reactome' in i:
            reactome += [i]
        else:
            other += [i]
    return genes,kegg,reactome,other

gene_files,kegg_files,reactome_files,other_files = split_files()

max_ind = ""
max_val = 0
for i in pws_per_param:
    if pws_per_param[i] > max_val:
        max_ind = i
        max_val = pws_per_param[i]

max_file = [f for f in kegg_files if max_ind in f and 'histogram' in f][0]

max_pws_parameter_histogram = read_histogram_file(max_file)
all_pws_param = {}
xy = {'x':kegg_parameter_histogram,'y':max_pws_parameter_histogram}
for i in xy:
    print "determinating",i,"values"
    for j in xy[i]:
        for pw in xy[i][j]:
            try:
                all_pws_param[pw][i] = j
            except:
                all_pws_param[pw] = {i:j}


##### Combine the actual cancer info with each coordinate file
dir = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists'
os.chdir(dir)

with open('ndel_-2_percdel_0.02_mes_1.2_mot_0.06_coordinates.tsv','r') as fh:
    gp_coord = {}
    head = fh.readline()
    for line in fh:
        sep = line.strip().split('\t')
        gp_coord[sep[0]] = sep[1:3]

with open('stats_out_revised/gene_pairs_ndel_-2_percdel_0.02_mes_1.2_mot_0.06.tsv','r') as fh:
    gp_cancers = {}
    for line in fh:
        sep = line.strip().split('\t')
        gp_cancers[sep[0]] = sep[1:]

## Check and see that both have the same genes
diff = []
if len(gp_coord) > len(gp_cancers):
    for i,gp in enumerate(gp_coord):
        try:
            assert gp_cancers[gp]
        except:
            print gp,i,"not in gp_cancers"
            diff.append(gp)
else:
    for i,gp in enumerate(gp_cancers):
        try:
            gp_coord[gp]
        except:
            print gp,i,"not in gp_coord"

for i in diff:
    del gp_coord[i]

of = open('ndel_-2_percdel_0.02_mes_1.2_mot_0.06_coordinates_and_cancers.tsv','w')
for gp in gp_coord:
    gp_coord[gp].extend(gp_cancers[gp])
    of.write(gp+'\t'+"\t".join(gp_coord[gp])+'\n')

of.close()


#### combine cancer info into pathway coordinate file too
#1) gather all pathways (both kegg and reactome) and make a dictionary of dictionaries (2) of lists
# {pw : { cancers : [...] , p_values : [] } } for all pws
#2) grab all p_values across all table_3.txt (a dictionary of lists)
#3) read in the pan_parameter histogram values for kegg and reactome pw_pairs
#4) take the geometric mean of each pathway
#5) make table with following info:
"""
pw1 pw2 prevalence  geom_mean_1 geom_mean_2
"""
# The table just generated can be used to sort pw_pairs by their combined geometric mean and prevalence 
# pw_pairs with high prevalence and significance would be good candidates for survival analysis
dbs = {'kegg':kegg_files,'reactome':reactome_files}




### count frequency of prevalence for individual genes across all parameters and cancers
import os
dir = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists'
os.chdir(dir)

dirs = [d for d in os.listdir('.') if '.new' in d]

genes = {}
for d in dirs:
    for f in os.listdir(d):
        print f
        p = f.split('.new_')[1].split('.tsv')[0]
        with open(d+'/'+f,'r') as fh:
            for line in fh:
                sep = line.strip('\n').split('\t')
                for i in sep[2:]:
                    try:
                        genes[i]['params'] += [p]
                        if d not in genes[i]['cancs']:
                            genes[i]['cancs'] += [d]
                    except:
                        genes[i] = {'params':[p],'cancs':[d]}

of = open('individual_gene_histogram.tsv','w')
for g in genes:
    of.write("\t".join([g,str(len(genes[g]['params'])),str(len(genes[g]['cancs'])),str(",".join(genes[g]['cancs'])),str(len(genes[g]['params'])*len(genes[g]['cancs']))])+"\n")

of.close()

