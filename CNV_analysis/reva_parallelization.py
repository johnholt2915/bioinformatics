import os
import sys
import time
from itertools import combinations
import subprocess
from math import ceil

def room_for_data(dir):
    p = subprocess.Popen(["df",dir],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    out,err = p.communicate()
    for s in out.split(' '):
        try:
            sp = s.split('%')
            if len(sp) > 1:
                if int(sp[0]) > 95:
                    return False
                else:
                    return True
        except:
            continue
    return False
    
def gen_bsub(n_cores,chunk,canc_data_file,pp_gp_file,lsf_file,outdir,bsub_out_dir,workdir):
    parallelize_script = "/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data/scripts/test_specific_genes_3.R"
    out = open(lsf_file,'w')
    out.write('#BSUB -J '+canc_data_file+'_chunk_'+str(chunk)+'\n')
    out.write('#BSUB -W 3:00\n')
    out.write('#BSUB -n '+str(n_cores)+'\n')
    out.write('#BSUB -m manda\n')
    out.write('#BSUB -q premium\n')
    out.write('#BSUB -P DrugSensitivity\n')
    out.write('#BSUB -o '+os.path.join(bsub_out_dir,canc_data_file+'_chunk_'+str(chunk)+'_bsub.out')+'\n')
    out.write('#BSUB -e '+os.path.join(bsub_out_dir,canc_data_file+'_chunk_'+str(chunk)+'_bsub.err')+'\n')
    out.write('#BSUB -u john.holt@mssm.edu\n')
    out.write('#BSUB -L /bin/bash\n')
    out.write('#BSUB -R rusage[mem=4200]\n')
    out.write('#BSUB -R span[hosts=1]\n')
    out.write('module purge\n')
    out.write('module load R\n')
    out.write('cd '+workdir+'\n')
    out.write('Rscript '+parallelize_script+' '+pp_gp_file+' '+canc_data_file+' '+outdir+' '+str(n_cores)+'\n')
    out.write('exit\n')
    out.close()
    print lsf_file
    #cmd = 'bsub < '+lsf_file
    #print runCommand(cmd)
    
def runCommand(command):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out,err


def jobs_done():
    p = subprocess.Popen(["bjobs"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = p.communicate()
    if out == '':
        return True
    return False

def gzipper(chunk,outdir,canc_data_file):
    zippable_files = [os.path.join(os.path.join(outdir,canc_data_file),f) for f in os.listdir(os.path.join(outdir,canc_data_file)) if ".gz" not in f and "chunk_"+str(chunk)+'_' not in f]
    if len(zippable_files) > 0:
        for zippable in zippable_files:
            print "gzipping",zippable
            subprocess.call("gzip "+zippable,shell=True)
        print "done"

def main(argv):
    ######## ARGS ########
    ## a file containing all gene names considered for comparison (n-genes means n-lines in the file)
    gene_file = argv[1]
    genes = [l.split('\n')[0] for l in open(gene_file,'r')]
    n_genes = len(genes)
    n_genes = 24776
    n_gps = ((n_genes*n_genes)/2) - (n_genes/2)
    ## a directory where the final R processed data will be output (these files wil be gzipped to make the most room)
    outdir = argv[2]
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    processing_cancers = [d for d in os.listdir(outdir) if d != 'preprocessing_gp_dir']
    ## the number of cores to use for processing genes. This will dictate the size of the chunks generated and (hopefully) the processing time per cancer
    canc_data_files = argv[3].split(',')#os.path.basename(argv[3])
    n_cores = int(argv[4])

    ## Per hour processing rates for each cancer
    proc_rates = {'ucec.cnv_clinical_ucec.new': 370523, 'kirp.cnv_clinical_kirp.new': 327854, 'lgg.cnv_clinical_lgg.new': 459662, 'kirc.cnv_clinical_kirc.new': 296039, 'stad.cnv_clinical_stad.new': 291969, 'blca.cnv_clinical_blca.new': 284861, 'paad.cnv_clinical_paad.new': 349989, 'prad.cnv_clinical_prad.new': 284887, 'esca.cnv_clinical_esca.new': 297802, 'lihc.cnv_clinical_lihc.new': 295085, 'coadread.cnv_clinical_coadread.new': 263153, 'skcm.cnv_clinical_skcm.new': 353356, 'ov.cnv_clinical_ov.new': 262881, 'brca.cnv_clinical_brca.new': 216454, 'sarc.cnv_clinical_sarc.new': 390861}
    chunk_size = (max([proc_rates[canc] for canc in canc_data_files]))*(n_cores-1)
    print "Chunk size:",chunk_size
    #sys.exit()
    
    ######## Directory Setup ########
    workdir = '/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data'
    os.chdir(workdir)
    if not os.path.isdir(os.path.join(workdir,"parallelize")):
        os.mkdir(os.path.join(workdir,"parallelize"))
    ## this folder will contain the gene pairs output from this script before any R processing occurs
    pp_gp_dir = os.path.join(workdir,"parallelize/preprocessing_gp_dir")
    if not os.path.isdir(pp_gp_dir):
        os.mkdir(pp_gp_dir)
    backup_pp_gp_dir = os.path.join(outdir,"preprocessing_gp_dir")
    if not os.path.isdir(backup_pp_gp_dir):
        os.mkdir(backup_pp_gp_dir)
    lsf_dir = os.path.join(workdir,"parallelize/parallelized_lsfs")
    if not os.path.isdir(lsf_dir):
        os.mkdir(lsf_dir)
    bsub_out_dir = os.path.join(outdir,'bsub_logs')
    if not os.path.isdir(bsub_out_dir):
        os.mkdir(bsub_out_dir)
    
    ######## Main Algo ########
    ## generate chunks of all pairwise combinations of genepairs one chunk at a time
    ## process each chunk with an R script utilizing 20 cores
    gp_gen = combinations(genes,2)
    n_chunks = int(ceil((float(n_gps)/(n_genes-1))/n_cores))
    for chk in range(n_chunks):
        #pp_gp_files = []
        lsf_files = []
        for canc_data_file in canc_data_files:
            if room_for_data(pp_gp_dir):
                #pp_gp_files.append(os.path.join(pp_gp_dir,canc_data_file+"_gps_chunk_"+str(chk)+".csv"))
                pp_gp_file = os.path.join(pp_gp_dir,canc_data_file+"_gps_chunk_"+str(chk)+".csv")
                lsf_files.append(os.path.join(lsf_dir,canc_data_file+"_gps_chunk_"+str(chk)+".lsf"))
            elif room_for_data(backup_pp_gp_dir):
                #pp_gp_file.append(os.path.join(backup_pp_gp_dir,canc_data_file+"_gps_chunk_"+str(chk)+".csv"))
                pp_gp_file = os.path.join(backup_pp_gp_dir,canc_data_file+"_gps_chunk_"+str(chk)+".csv")
                lsf_file.append(os.path.join(lsf_dir,canc_data_file+"_gps_chunk_"+str(chk)+".lsf"))
            else:
                print "No Room for data ... exiting at chunk",chk,"having processed",(chk-1)*chunk_size,"gene pairs for",canc_data_file
                print "df",pp_gp_dir
                print "df",backup_pp_gp_dir
                sys.exit()
        pair_count = 1
        print "writing gene pairs to",pp_gp_file
        gps_out = open(pp_gp_file,'w')
        while pair_count%chunk_size != 0:
            gps_out.write(",".join(gp_gen.next())+'\n')
            pair_count += 1
        gps_out.close()
        for i,canc_data_file in enumerate(canc_data_files):
            print "n_cores:",n_cores
            print "chunk:",chk
            print "canc_data_file:",canc_data_file
            print "gene_pair_file:",pp_gp_file
            print "lsf_file:",lsf_files[0]
            print "outdir:",outdir
            print "bsub_out_dir:",bsub_out_dir
            gen_bsub(n_cores,chk,canc_data_file,pp_gp_file,lsf_files[i],outdir,bsub_out_dir,workdir)
        sys.exit()
        n_waits = 0
        while not jobs_done():
            # gzip data output from last chunk to make room for more ...
            #for canc_data_file in canc_data_files:
            #    gzipper(chk-1,outdir,canc_data_file)
            # 5 minutes of sleep before checking if jobs are done...
            print "Waited",n_waits*5,"minutes"
            time.sleep(5*60)
            n_waits += 1
        break
    

if __name__ == '__main__':
    main(sys.argv)

