##################################
#                                #
# Last modified 2017/04/17       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam
import math

def run():

    if len(sys.argv) < 4:
        print('usage: python %s bedfilename chrField BAMfilename outputfilename [-nomulti] [-readLength min max] [-mismatchesMD M] [-mismatches M] [-uniqueBAM] [-singleFieldRegion] [-end2only] [-end1only] [-stranded +|-] [-namePrefix string]' % sys.argv[0])
        print('\tthe script will divide multireads by their multiplicity')
        print('\tuse the uniqueBAM option if the BAM file contains only unique alignments and you have no NH tags in your file; this option should never be used for files with multimappers')
        print('\totherwise the presence of NH tags is assumed')
        print('\tthe script only works on sorted and index BAM files')
        sys.exit(1)
    
    bed = sys.argv[1]
    fieldID = int(sys.argv[2])
    BAM = sys.argv[3]
    outfilename = sys.argv[4]

    doEnd1Only = False
    doEnd2Only = False
    if '-end1only' in sys.argv and '-end2only' in sys.argv:
        print('both -end1only and -end2only option specified, a logical impossiblity, exiting')
        sys,exit(1)

    if '-end1only' in sys.argv:
        doEnd1Only = True
        print('will only consider the first end of read pairs')

    if '-end2only' in sys.argv:
        doEnd2Only = True 
        print('will only consider the second end of read pairs')
    
    doStranded=False
    if '-stranded' in sys.argv:
        doStranded=True
        thestrand = sys.argv[sys.argv.index('-stranded')+1]
        print('will only consider', thestrand, 'strand reads')

    doNamePrefix = False
    if '-namePrefix' in sys.argv:
        doNamePrefix = True
        namePrefix = sys.argv[sys.argv.index('-namePrefix') + 1]
 
    doSFR = False
    if '-singleFieldRegion' in sys.argv:
        doSFR = True

    noMulti=False
    if '-nomulti' in sys.argv:
        noMulti=True
        print('will discard multi-read alignments')

    doReadLength=False
    if '-readLength' in sys.argv:
        doReadLength=True
        minRL = int(sys.argv[sys.argv.index('-readLength')+1])
        maxRL = int(sys.argv[sys.argv.index('-readLength')+2])
        print('will only consider reads between', minRL, 'and', maxRL, 'bp length')
        ORLL = 0

    doMaxMMMD=False
    if '-mismatchesMD' in sys.argv:
        doMaxMMMD=True
        maxMM = int(sys.argv[sys.argv.index('-mismatchesMD')+1])
        print('Will only consider alignments with', maxMM, 'or less mismatches')

    doMaxMM=False
    if '-mismatches' in sys.argv:
        doMaxMM=True
        maxMM = int(sys.argv[sys.argv.index('-mismatches')+1])
        print('Will only consider alignments with', maxMM, 'or less mismatches')

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print('will treat all alignments as unique')
        doUniqueBAM = True
        TotalReads = 0
        pass
        
    regionDict={}

    outfile = open(outfilename, 'w')

    samfile = pysam.Samfile(BAM, "rb" )
    lineslist = open(bed)
    i=0
    for line in lineslist:
        i+=1
        if i % 10000 == 0:
            print(i, 'regions processed')
        if line[0]=='#':
            continue
        fields = line.strip().split('\t')
        if len(fields) < fieldID+2:
            continue
        chr = fields[fieldID]
        left = int(fields[fieldID+1])
        right = int(fields[fieldID+2])
        if left >= right:
            print('problem with region, skipping:', chr, left, right)
            continue
        reads=0
        try:
            for alignedread in samfile.fetch(chr, left, right):
                works = True
                break
        except:
            reads = 0
            print('problem with region:', chr, left, right, 'assigning 0 value')
            if doSFR:
                outline = chr + ':' + str(left) + '-' + str(right) +' \t' + str(reads)
            else:
                outline = line.strip() +' \t' + str(reads)
            if doNamePrefix:
                outline = namePrefix + outline
            outfile.write(outline + '\n')
            continue
        for alignedread in samfile.fetch(chr, left, right):
            if doEnd2Only and alignedread.is_read1:
                continue
            if doEnd1Only and alignedread.is_read2:
                continue
            if doStranded:
                if alignedread.is_reverse:
                    s = '-'
                else:
                    s = '+'
                if s != thestrand:
                    continue
            if doUniqueBAM:
                multiplicity = 1
            else:
                multiplicity = alignedread.opt('NH')
            if doMaxMM:
                mismatches = 0
                for (m,bp) in alignedread.cigar:
                    if m == 8:
                        mismatches+=1
                if mismatches > maxMM:
                    continue
            if doMaxMMMD:
                MM = alignedread.opt('MD')
                mismatches = 0
                if MM.isdigit():
                    pass
                else:
                    for s in range(len(MM)):
                        if MM[s].isalpha():
                            mismatches+=1
                if mismatches > maxMM:
                    continue
            if doReadLength:
                if len(alignedread.seq) > maxRL or len(alignedread.seq) < minRL:
                    continue
            if noMulti and multiplicity > 1:
                continue
            reads += 1./multiplicity
#            print multiplicity, 1./multiplicity, alignedread.opt('NH'), str(alignedread)

        if doSFR:
            if doStranded:
                outline = chr + ':' + str(left) + '-' + str(right) + '|' + thestrand + '\t' + str(reads)
            else:
                outline = chr + ':' + str(left) + '-' + str(right) +' \t' + str(reads)
        else:
            outline = line.strip() +' \t' + str(reads)
        if doNamePrefix:
            outline = namePrefix + outline
        outfile.write(outline + '\n')
         
    outfile.close()
   
run()

