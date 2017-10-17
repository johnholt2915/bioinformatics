from igv import IGV
import sys
import os
'''
Author: John Holt
Date: April 10th, 2017
Purpose: This script is designed to loop through a file of regions (chromosome\tposition\tbam_url format) line by line and to take screen shots of these regions in igv.
'''

def visualizeRegions(igv,regions,outdir):
    d = {}
    with open(regions,'r') as fh:
        fh.readline()
        for line in fh:
            sep = line.split()
            chr = sep[0]
            pos = sep[1]
            bam_url = sep[2]
            try:
                d[chr][pos].append(bam_url)
            except:
                try:
                    d[chr][pos] = [bam_url]
                except:
                    d[chr] = dict()
                    d[chr][pos] = [bam_url]
    print d

    for chr in d:
        for pos in d[chr]:
            load = raw_input("Type 'y' to continue or 'n' to quit:")
            igv.clear()
            if load == 'y':
                window_min = str(int(pos)-20)
                window_max = str(int(pos)+90)
                rev_q = ''
                if '_revcomp' in d[chr][pos][0]:
                    rev_q = '_revcomp'
                    window_min = str(int(pos)-70)
                    window_max = str(int(pos)+40)
                igv.go(chr+":"+window_min+"-"+window_max)
                # load all urls for each position of interest
                for url in d[chr][pos]:
                    igv.load(url)
                    print "loaded",url
                igv.save(outdir+"/"+chr+"_"+pos+rev_q+".png")
                print "saved screenshot to",outdir+"/"+chr+"_"+pos+rev_q+".png"
            else:
                print "quitting"
                sys.exit()

def main():
    # load the file of regions in (chromosome\tposition\tbam_url...) format as well as the output folder for screen shots.
    regions = sys.argv[1]
    outdir = sys.argv[2]
    if not os.path.isdir(outdir):
        print "mkdir %s" % outdir
        os.mkdir(outdir)

    # connect to local (running) igv.
    igv = IGV()
    try:
        igv.clear()
        igv.genome('hg19')
    except:
        print "start IGV first, then run this script!"
        sys.exit()

    # two urls are passed as input to compare mulitple bam files.
    #visualizeRegions(igv,regions,url1a,url1b,url2a,url2b,outdir)
    visualizeRegions(igv,regions,outdir)

if __name__ == '__main__':
    main()


