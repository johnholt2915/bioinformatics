from Bio.Blast import NCBIWWW, NCBIXML
import math
import random
import re
import datetime
import optparse
import sys
import os
import operator
import xml.etree.ElementTree

def argParser():
    parser = optparse.OptionParser()
    parser.add_option( '-f', '--fastq', dest='fastq')
    parser.add_option( '-s', '--sample-size', dest='samp_size')
    ( options, args ) = parser.parse_args()

    try:
        if options.samp_size == None:
            raise Exception
        if options.fastq == None:
            raise Exception
    except Exception:
        print "Usage: blastReads.py -f <fastq-file> -s <sample-size"
        sys.exit()

    return options

def getRandomLines(filename,n_lines):
    filesize = os.path.getsize(filename)
    fastq = open(filename,'r')
    fasta_lines = {}
    count = 0
    output = []
    while count < n_lines:
        offset = random.randrange(filesize)
        #print "offset:",offset
        fastq.seek(offset)
        fastq.readline()
        found = False
        while not found:
            line = fastq.readline()
            split = line.split(":")
            if len(split) == 10 and split[0][0] == '@':
                found = True
                seq_id = re.sub("@",">",line)
                try:
                    assert fasta_lines[seq_id]
                except:
                    fasta_lines[seq_id] = fastq.readline()
                    output.extend([seq_id,fasta_lines[seq_id]])
                    count += 1
    fasta_file_out = outputFileName(filename,"fa")
    fh = open(fasta_file_out,'w')
    fh.write("".join(output))
    fh.close()
    return "".join(output)

# older implementation ... takes linear time with respect to the file size, plus additional time to generate the indeces.
def grabLines(fastq,indeces):
    fasta = []
    with open(fastq,"r") as fq:
        line_num = 0
        found = 0
        print "found: "+str(found)
        line = fq.readline()
        seq_id = ""
        while len(indeces[found:]) > 0:
            if line_num%4 == 0:
                seq_id = line
            if line_num%4 == 1:
                convert_index = (line_num+3)/4 # because fastq files are 4 standardized line per read, the indeces are 4 times smaller than actual number of lines in the file ...
                if indeces[found] == convert_index:
                    seq_id = re.sub("@",">",seq_id)
                    fasta.extend([seq_id,line])
                    found += 1
                    print "found: "+str(found)
            line_num += 1
            line = fq.readline()
    return "".join(fasta)

def runBlast(fasta):
    xml_out = outputFileName(fasta,"xml")
    #cmd = "blastn -query %s -db nr -out %s -evalue 0.001 -outfmt 5"
    #subprocess.call(cmd,shell=True)
    result_handle = NCBIWWW.qblast("blastn","nt",fasta,expect=0.01)
    return result_handle

def saveXML(blast_xml,output_file):
    xml_file = open(output_file,"w") 
    xml_file.write(blast_xml.read())
    xml_file.close()

def calculateAlignedTo(align_str):
    """This function takes an alignment.title string from the qblast xml output and determines if the species is human, carp, or other."""
    human = ["human","homo","sapiens"]
    carp = ["carp","cyprinus","carpio"]
    words = align_str.split('|')[4]
    words = words.split()
    for word in words:
        word = word.replace(',','')
        word = filter(str.isalpha,word)
        if word.lower() in human:
            return "Human"
        elif word.lower() in carp:
            return "Carp"
    return "Other"

def writeAlignments(fastq):
    outfile = outputFileName(fastq,"alignments")
    xml_out = outputFileName(fastq,"xml")
    blast_xml = open(xml_out,'r')
    records = NCBIXML.parse(blast_xml)
    species_dict = {"Human":0,"Carp":0,"Other":0}
    
    return species_dict

def parseBlast(blast_xml,fastq,save=True):
    """Display only the top 10 alignments for each query sequence."""
    outfile = outputFileName(fastq,"alignments")
    xml_out = ""
    if save:
        # saving the blast results to xml file.
        xml_out = outputFileName(fastq,"xml")
        print "saving BLAST results to",xml_out
        saveXML(blast_xml,xml_out)
        print "SAVED RESULTS TO",xml_out
    else:
        xml_out = blast_xml
    # loading the blast results from xml file
    print "loading alignments from",blast_xml
    blast_xml = open(xml_out,"r")
    print "LOADED ALIGNMENTS"
    records = NCBIXML.parse(blast_xml)
    species_dict = {"Human":0,"Carp":0,"Other":0}
    print "writing alignments to",outfile
    with open(outfile,"w") as op:
        op.write("Alignment records for "+fastq+" with e-value threshold of 0.01\n")
        op.write("\tAlignment date: "+str(datetime.datetime.now())+"\n")
        op.write("\n\n")
        for record in records:
            op.write("********************************************************************************\n")
            op.write(str(record.query)+"\n")
            tmp_species_dict = {"Human":0,"Carp":0,"Other":0}
            align_count = 1
            # loop through all alignments and determine what species they came from using 'calculateAlignedTo()'
            # save these values in tmp_species_dict
            for alignment in record.alignments:
                if align_count <= 3:
                    op.write("** ALIGNMENT "+str(align_count)+" **\n")
                    #if alignment.hsps[0].score < 150.0:
                    #    continue
                    aligned_to = calculateAlignedTo(str(alignment.title))
                    tmp_species_dict[aligned_to] += 1
                    count = 1
                    for hsp in alignment.hsps:
                        if count <= 10:
                            op.write("********HSP "+str(count)+"********\n")
                            op.write("Aligned To: "+str(alignment.title)+"\n")
                            op.write("Evalue: "+str(hsp.expect)+"\n")
                            op.write("Score: "+str(hsp.score)+"\n")
                            op.write(str(hsp.query)+"\n")
                            op.write(str(hsp.match)+"\n")
                            op.write(str(hsp.sbjct)+"\n")
                            op.write("\n")
                            count += 1
                        else:
                            break
                    align_count += 1
                else:
                    break
            # determine the key within tmp_species_dict which has the most alignments
            # increment species_dict for the max_key
            max_key = max(tmp_species_dict.iteritems(), key=operator.itemgetter(1))[0]
            species_dict[max_key] += 1
    return species_dict

def outputFileName(fastq,appendix):
    fq_split = fastq.split(".")
    fq_split[-1] = appendix
    return ".".join(fq_split)

def makeWordHist(fastq,word_hist,samp_size,count=False):
    outfile = outputFileName(fastq,'hist.txt')
    print word_hist
    if count:
        alignment_file = outputFileName(fastq,"alignments")
        word_hist['total'] = 0
        with open(alignment_file,'r') as fh:
            for line in fh:
                if "Aligned To:" in line:
                    word_hist['total'] += 1
                    sep = line.split('|')[4]
                    sep2 = sep.split()
                    for word in sep2:
                        word = word.replace(',','')
                        if word == '':
                            continue
                        try:
                            assert int(word)
                            continue
                        except:
                            pass
                        try:
                            word_hist[word] += 1
                        except:
                            word_hist[word] = 1
    of = open(outfile,'w')
    sum = 0
    for key in word_hist:
        sum += int(word_hist[key])
    for key in word_hist:
        of.write(key+","+str(float(word_hist[key])/sum)+"\n")
    of.close()

def main():
    options = argParser()
    fastq = options.fastq
    samp_size = int(options.samp_size)
    xml_file = outputFileName(fastq,"xml")
    alignment_hist = {}
    if os.path.exists(xml_file):
        print xml_file,"exists..."
        alignment_hist = parseBlast(xml_file,fastq,save=False)
    else:
        print "Sampling",samp_size,"reads from",fastq
        fasta = getRandomLines(fastq,samp_size)
        print fasta
        #sys.exit()
        # pass fasta sequences to run blast
        print "running BLAST... This may take sometime..."
        blast_xml = runBlast(fasta)
        print "ran blast... now writing results"
        # write and return results to present all alignments and the species to which they are aligned
        alignment_hist = parseBlast(blast_xml,fastq,save=True)

    # word histogram of sequences aligned to
    makeWordHist(fastq,alignment_hist,samp_size)

if __name__ == '__main__':
    main()

