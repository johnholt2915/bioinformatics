import os
import sys
import glob

###############################
decryption_key = "vPL7ickF5ETF5cSjYlGcWBvELW4wCCnv"
n_lsfs = 2

###############################

def lsf_header(fh,out_dir,step_name,job_name,threads,wall_time,rusage,acc='acc_cncrseq')
    fh.write('#!/bin/bash\n')
    fh.write('#BSUB -J %s\n' % (job_name))
    fh.write('#BSUB -W %s\n' % (str(wall_time)))
    fh.write('#BSUB -n %d\n' % (threads))
    fh.write('#BSUB -m manda\n')
    fh.write('#BSUB -q premium\n')
    fh.write('#BSUB -P %s\n' % (acc))
    fh.write('#BSUB -o %s\n' % (os.path.join(out_dir,str(step)+'_bsub.out')))
    fh.write('#BSUB -e %s\n' % (os.path.join(out_dir,str(step)+'_bsub.err')))
    fh.write('#BSUB -u john.holt@mssm.edu\n')
    fh.write('#BSUB -R rusage[mem=%d]\n' % (rusage))
    fh.write('#BSUB -R span[host=1]\n')
    fh.write('#BSUB -L /bin/bash\n')

def decrpytion_lsf(sample,out_dir,key=decryption_key):
    of = open(os.path.join(out_dir,'decrypt.lsf'),'w')
    lsf_header(of,out_dir,'decrypt',out_dir.split('Sample_')[1]+"_decrypt",4,'1:00',4200)
    output_sample_name = os.path.basename(sample).split('.enc')[0]
    of.write('openssl enc -d -aes-256-cbc -in %s -out %s -k %s' % (sample,output_sample_name,key))
    of.close()
    return output_sample_name

def bam2fastq_lsf(out_dir,bam):
    of = open(os.path.join(out_dir,'bam2fastq.lsf'),'w')
    sample = out_dir.split('Sample_')[1]
    lsf_header(of,out_dir,'bam2fastq',sample+'_bam2fq',4,'12:00',16800)
    of.write('cd %s\n' $ (out_dir))
    of.write('module purge all\n')
    of.write('module load picard\n')
    r1_fastq = bam.split('.bam')[0]+"_R1.fastq"
    r2_fastq = bam.split('.bam')[0]+"_R2.fastq"
    of.write('echo "java -jar $PICARD SamToFastq I=%s FASTQ=%s SECOND_END_FASTQ=%s VALIDATION_STRINGENCY=LENIENT"\n' % (bam,r1_fastq,r2_fastq))
    of.write('java -jar $PICARD SamToFastq I=%s FASTQ=%s SECOND_END_FASTQ=%s VALIDATION_STRINGENCY=LENIENT\n' % (bam,r1_fastq,r2_fastq))
    of.close()
    return r1_fastq,r2_fastq

def wgs_pipeline_lsfs(out_dir,r1s,r2s,acc='acc_cncrseq'):
    sample = out_dir.split('Sample_')[1]
    xml = open(os.path.join(out_dir,sample+'.xml'),'w')
    xml.write('<?xml version="1.0" ?>\n')
    xml.write('<cohort library="%s" platform="illumina" platform_unit="">\n' % (sample))
    for i in range(len(r1s)):
        xml.write('  <fastq id="%s" sample="%s">\n' % (sample,sample))
        xml.write('    <fastq1 file="%s"/>\n' % (r1s[i]))
        xml.write('    <fastq2 file="%s"/>\n' % (r2s[i]))
        xml.write('  </fastq>\n')
    xml.write('</cohort>\n')
    xml.close()
    print "wrote",os.path.join(out_dir,sample+".xml")
    makefile = open(os.path.join(out_dir,"Makefile"),'w')
    makefile.write('PROJECT=%s\n' % (sample))
    makefile.write('INPUT=%s.xml\n' % (sample))
    makefile.write('PIPELINE=GENOME\n')
    makefile.write('JOB_ACCT=%s\n' % (acc))
    makefile.write('JOB_QUEUE=premium\n')
    makefile.write('SCATTERGATHER=32\n')
    makefile.write('THREADS=20\n')
    makefile.write("RUN_ARGS=-jobNative '-m manda'\n")
    makefile.write('REFERENCE=hg19\n')
    makefile.write('WALLTIME=$$(( 96 * 60 * 60 ))\n')
    makefile.write('include $(NGSMAKE_ROOT)/workflow.mk\n')
    makefile.close()
    print "wrote",os.path.join(out_dir,"Makefile")
    lsf = open(os.path.join(out_dir,'runngs.lsf'),'w')
    lsf.write('#!/bin/bash\n')
    lsf.write('#BSUB -J %s\n' % (sample))
    lsf.write('#BSUB -W 96:00\n')
    lsf.write('#BSUB -n 4\n')
    lsf.write('#BSUB -m manda\n')
    lsf.write('#BSUB -q premium\n')
    lsf.write('#BSUB -P %s\n' % (acc))
    lsf.write('#BSUB -o %s\n' % (os.path.join(out_dir,'ngs_bsub.out')))
    lsf.write('#BSUB -e %s\n' % (os.path.join(out_dir,'ngs_bsub.err')))
    lsf.write('#BSUB -u john.holt@mssm.edu\n')
    lsf.write('#BSUB -L /bin/bash\n')
    lsf.write('cd %s\n' % (out_dir))
    lsf.write('module purge all\n')
    lsf.write('module load ngs/3.2.0\n')
    lsf.write('make alignmentQCSummary.tsv variantMetrics\n')
    lsf.close()
    print "wrote",os.path.join(out_dir,'runngs.lsf')

def gzip_fastqs_lsf(out_dir,fq_list,read):
    gzip_submitter = open(os.path.join(out_dir,'gzip_'+read.lower()+'s_submitter.sh'),'w')
    gzip_submitter.write('#!/bin/bash\n')
    gzip_out = []
    for fq in fq_list:
        fq_id = fq.split('_'+read.lower()+'_')[1].split('_'+read+'.fastq')[0]
        of = open(os.path.join(out_dir,'gzip_'+read.lower()+'_'+fq_id+'.lsf'),'w')
        gzip_submitter.write('bsub < %s\n' % (os.path.join(out_dir,'gzip_'+read.lower()+'_'+fq_id+'.lsf')))
        sample = out_dir.split('Sample_')[1]
        lsf_header(of,out_dir,'gzip_'+read.lower()+'_'+fq_id,sample+'_gz_'+read.lower()+'_'+fq_id,4,'6:00',8400)
        of.write('cd %s\n' % (out_dir))
        of.write('gzip %s\n' % (fq))
        of.close()
        gzip_out = fq +".gz"
    gzip_submitter.close()
    return gzip_out

def rounder(val):
    if (val - int(val)) < 0.5:
        return int(val)
    else:
        return int(val)+1

def split_corr_fastqs(sample,n_paired_end_fastqs):
    r1_out = sample+"_r1_"
    r2_out = sample+"_r2_"
    ab = 'abcdefghijklmnopqrstuvwxyz'
    r1s = []
    r1s_corr = []
    r2s = []
    r2s_corr = []
    tmp_len = n_paired_end_fastqs
    for i in range((n_paired_end_fastqs/len(ab))+1):
        if tmp_len >= len(ab):
            inner_loop = 26
        else:
            inner_loop = tmp_len % len(ab)
        for j in range(inner_loop):
            r1s.append(r1_out+ab[i]+ab[j])
            r1s_corr.append(r1s[-1]+"_R1.fastq")
            r2s.append(r2_out+ab[i]+ab[j])
            r2s_corr.append(r2s[-1]+"_R2.fastq")
        tmp_len -= inner_loop
    return r1s, r2s, r1s_corr, r2s_corr

def splitfastq_lsf(out_dir,r1_fastq,r2_fastq,bam_enc,n_lines=400000000):
    bam_size = os.path.getsize(bam_enc)
    n_paired_end_fastqs = (8 * bam_size) / 42724671487.904564
    n_paired_end_fastqs = rounder(n_paired_end_fastqs)
    sample = r1_fastq.split('_R1.fastq')[0]
    r1_split_fqs, r2_split_fqs, corr_r1_fqs, corr_r2_fqs = split_corr_fastqs(sample,n_paired_end_fastqs)
    sample = os.path.basename(r1_fastq).split('_R1.fastq')[0]
    r1_out = sample+"_r1_"
    r2_out = sample+"_r2_"
    of = open(os.path.join(out_dir,'split_fastq_r1.lsf'),'w')
    lsf_header(of,out_dir,'split_fastq_r1',sample+'_split_r1',4,'4:00',4200)
    of.write('cd %\n' % (out_dir))
    of.write('split -l %d %s %s\n' % (n_lines,r1_fastq,r1_out))
    for i in range(len(r1s)):
        of.write('mv %s %s\n' % (r1s[i],r1s_corr[i]))
        of.write('mv %s %s\n' % (r2s[i],r2s_corr[i]))
    return r1s_corr, r2s_corr

def get_batched_files(raw_dir):
    batches = [os.path.join(b) for b in os.listdir(raw_dir)]
    samples_processed = []
    for b in batches:
        samples = glob.glob(os.path.join(b,"batch*/Raw/*/Sample_*"))
        for s in samples:
            lsfs = [f for f in os.listdir(s) if f[-4:] == '.lsf']
            if len(lsfs) == n_lsfs:
                samples_process.append(s)
    return samples_processed

def get_batch(download_dir,raw_dir,processed_files,batch_size=10):
    files = [os.path.join(download_dir,f) for f in os.listdir(download_dir)]
    files.sort(key=lambda x: os.path.getmtime(x))
    batch = files[:batch_size]
    return batch

def setup_batch_folder(batch_to_process,raw_dir,processed_dir):
    processed_batches = [int(d.split('batch')[0]) for d in os.listdir(raw_dir) if "batch" in d]
    next_batch = max(processed_batches) + 1
    batch_dir = os.path.join(raw_dir,'batch'+str(next_batch),'Raw','DNA.pipeline.WGS')
    processed_batch_dir = os.path.join(processed_dir,'batch'+str(next_batch))
    if not os.path.isdir(batch_dir):
        os.makedirs(batch_dir)
    if not os.path.isdir(processed_batch_dir):
        os.mkdir(processed_batch_dir)
    for s in batch_to_process:
        d = os.path.join(batch_dir,"Sample_"+os.path.basename(s).split('.enc')[0])
        os.mkdir(d)
    print "setup",batch_dir,"complete"
    return batch_dir,processed_batch_dir

def vcf_copy_sh(out_dir,processed_dir):
    of = open(os.path.join(out_dir,'vcf_copy.sh'),'w')
    of.write('#!/bin/bash\n')
    of.write('cd %s\n' % (out_dir))
    of.write('for i in $(find . -name *combined.vcf); do\n')
    of.write('  echo "copying $i"\n')
    of.write('  cp $i %s/.\n' % (processed_dir))
    of.write('done\n')

def create_lsf_scripts(batch,batch_dir,processed_dir):
    for s in batch:
        out_dir = os.path.join(batch_dir,"Sample_"+os.path.basename(s))
        if not os.path.exists(out_dir):
            print out_dir,"does not exist for",s
            continue
        bam = decryption_lsf(s,out_dir)
        r1_fastq,r2_fastq = bam2fastq_lsf(out_dir,bam)
        r1s_out,r2s_out = splitfastq_lsf(out_dir,r1_fastq,r2_fastq,s)
        gzip_r1s = gzip_fastqs_lsfs(out_dir,r1s_out,'R1')
        gzip_r2s = gzip_fastqs_lsfs(out_dir,r2s_out,'R2')
        wgs_pipeline_lsfs(out_dir,gzip_r1s,gzip_r2s)
        vcf_copy_sh(out_dir,processed_dir)

def main(argv):
    download_dir = argv[1]
    raw_dir = argv[2]
    processed_dir = argv[3]
    batch_size = 10

    processed_files = get_batched_files(processed_dir,raw_dir)
    batch_to_process = get_batch(download_dir,raw_dir,processed_files,batch_size)
    raw_batch_folder,processed_batch_folder = setup_folders(batch_to_process,raw_dir,processed_dir)
    create_lsf_scripts(batch_to_process,raw_batch_folder,processed_batch_folder)
    print "Done"

if __name__ == '__main__':
    main(sys.argv)

