#!/usr/bin/env python

import sys
import os
import csv

# Take in genome.fa, RNA_reads.sra, RPF_reads.sra, and /Footer output
# directory as arguments in command line.
genome = sys.argv[1]
genomeShortcut = os.path.splitext(os.path.basename(genome))[0]

RNA = sys.argv[2]
RNAShortcut = os.path.splitext(os.path.basename(RNA))[0]
                               
RPF = sys.argv[3]
RPFShortcut = os.path.splitext(os.path.basename(RPF))[0]

outputDir = sys.argv[4]
outputShortcut = os.path.splitext(os.path.basename(outputDir))[0]

# Run bowtie to index genome
# Only needs to be done once per genome
# to improve efficiency we can make this happen
# Written by Shalyn
def run_bowtie(data):
    print ("Indexing... ")
    os.system("bowtie-build " + data + " " + "INDEX")
    print ("Done!\n")

# Changes into Tophat Out directory
# Converts accepted_hits.bam file into accepted_hits.sam file, moving them
# To the appropriate output directory
# Written by Shalyn
def convert(name):
    print("Converting .bam into .sam...")
    os.system("cd " + name + "_thout")
    os.system("samtools view -h -o accepted_hits.sam accepted_hits.bam")
    os.system("mv accepted_hits.bam ..") # +outputDir)
    os.system("mv accepted_hits.sam ..") #+ #outputDir)
    os.system("cd ..")
    print("Done!\n")

# Run Tophat on given RNA or RPF data, output to _thout directory
# Written by Shalyn
def run_tophat(data,name):
    print ("Running tophat on " + name+ "...")
    os.system("tophat -o "+name+"_thout INDEX "+data)
    convert(name)
    print("Done!\n")

def findRegions():
    regionFile = open(genomeShortcut + "-regionResults.txt", "w")
    with open(RPF) as regionFileName:
        for line in regionFileName:
            gatherLine = line.split()
            if(gatherLine[0] == '@SQ'):
                secRegionIndex = -1
                thirdRegionIndex = -1
                for blocks in xrange(len(gatherLine)):
                    if(gatherLine[blocks][0:2] == 'SN'):
                        secRegionIndex = blocks
                    elif(gatherLine[blocks][0:2] == 'LN'):
                        thirdRegionIndex = blocks
                if(secRegionIndex == -1 or thirdRegionIndex == -1):
                    continue
                secondRegion = gatherLine[secRegionIndex]
                regionFile.write(gatherLine[secRegionIndex][3:] + "\t")
                thirdRegion = gatherLine[thirdRegionIndex]
                regionFile.write(gatherLine[thirdRegionIndex][3:] + "\t")
                regionFile.write("\n")
    regionFile.close()

# Create file to write start/stop coordinates to of
# regions with continuous coverage
# Written by Naina
def create_blocks():
    print ("Creating coordinates of areas with continuous reads...")
    os.system("bedtools bamtobed -i accepted_hits.bam > reads.bed")
    os.system("sort -k1,1 -k2,2n reads.bed > reads.sort.bed")
    os.system("bedtools merge -i reads.sort.bed -d 10 > "+
              RPFShortcut + "-tempBlocks.txt")
    with open(RPFShortcut + "-tempBlocks.txt", "r") as f:
        with open(RPFShortcut + "-blocks.txt", 'a') as blocks:
            for line in f:
                if line.strip():
                    blocks.write("\t".join(line.split()[1:3]) + "\n")
    print("Done!\n")
    # Now have tempBlocks.txt file in directory
    
# Find read counts for each index in genome given
# Written by Alex
def find_counts():
    print("Finding read counts...")
    os.system("samtools sort accepted_hits.bam acceptedhits.sorted")
    os.system("bedtools genomecov -d -ibam acceptedhits.sorted.bam "+
            "-g "+ genomeShortcut +"-regionResults.txt > " + RPFShortcut + "-readCounts.txt")
    print("Done!\n")
    # Now have readCounts.txt file in directory
    
# Find read counts for each index in genome given
# Written by Alex/Shalyn
def new_counts():
    print("Finding read counts...")
    os.system("samtools sort " + RPFShortcut + '-28reads.bam' + RPFShortcut + "-28reads.sorted")
    os.system("bedtools genomecov -d -ibam "+RPFShortcut+"-28reads.sorted.bam "+ 
            "-g " + genomeShortcut +"-regionResults.txt > "+ RPFShortcut + "-readCounts.txt")
    print("Done!")
    # Now have readCounts.txt file in directory    
    
# Find peaks in reads given certain threshold
# Written by Jacob
def findPeaks(decimalThreshold):
    firstPassPeaks = open(RPFShortcut + "-possibleORFPeaks.txt", "w")
    with open(RPFShortcut + "-readCounts.txt") as findAvg:
        total = 0
        count = 0
        for line in findAvg:
            getReads = line.split()
            #Get reads for that base
	    if(int(getReads[2]) > 0):
            	total += int(getReads[2])
            	count += 1
	print total
	print count
        peakAverage = total/count
	print "AVERAGE"
	print peakAverage
	print "AVERAGE"
        findAvg.close()
    with open(RPFShortcut + "-readCounts.txt") as findAvg:
	lineAddCount = 0
        for line in findAvg:
            findPeaks = line.split()
            #Peaks defined as 50% higher than average
            threshold = decimalThreshold
            if(int(findPeaks[2]) >= threshold * peakAverage):
                #Found a peak
		lineAddCount += 1
                for index in xrange(len(findPeaks)):
                    firstPassPeaks.write(findPeaks[index] + "\t")
                firstPassPeaks.write("\n")
	print "LINES IN FILE"
	print lineAddCount
	print "LINES IN FILE"
        firstPassPeaks.close()
        findAvg.close()

def findORFs(charSearchRange):
    with open(genome) as theGenome:
        orgGenome = theGenome.read()
        genomeLength = len(orgGenome)
        with open(RPFShortcut + "-possibleORFPeaks.txt") as firstPassPeaks:
            #Go through each potential peak to look for the base at that site
	    lineNumber = 0
	    foundOrf = open(RPFShortcut + "-foundORFs.txt", "w")
    	    errorFile = open(RPFShortcut + "-orfFindErrors.txt", "w")
            for line in firstPassPeaks:
		lineNumber += 1
		print lineNumber
                findLoc = line.split()
                #Range to look - in case peak not completely accurate
                charRange = charSearchRange
                #Corrects for non-genome chars (headers) in TopHat file
                posCorrection=countNonGenomeChars(orgGenome,int(findLoc[1]))
                startIndex = int(findLoc[1]) + posCorrection
                startSearch = startIndex - charRange
                #In bounds
                if(startSearch < 0):
                    startSearch = 0
                endSearch = startIndex + charRange
                if(endSearch >= genomeLength):
                    endSearch = genomeLength - 1
                orfSearch = ""
                #Create substring to search within
                for nuc in xrange(startSearch,endSearch+1):
                    orfSearch += orgGenome[nuc]
                orfInRange = orfSearch.find("ATG")
                #Found an ATG in range!
                if(orfInRange != -1):
                    startSiteOffset = orfInRange - charRange
                    orfPos = startIndex + startSiteOffset
                    #Location is where in genome it is
                    foundOrf.write("ATG" + "\t" + str(orfPos) +
                                   "\t" + findLoc[2])
                    foundOrf.write("\n")
                else:
                    #Create an error log
                    for index in xrange(len(findLoc)):
                        errorFile.write(findLoc[index] + "\t")
                    errorFile.write("\n")
	    foundOrf.close()
	    errorFile.close()
	    print "DONE FINDING ORFS"
            firstPassPeaks.close()
        theGenome.close()
                    
def searchForORFInRange(orfSearch, startIndex, charRange, findLoc):
    orfInRange = orfSearch.find("ATG")
    #Found an ATG in range!
    foundOrf = open(RPFShortcut + "-foundORFs.txt", "w")
    errorFile = open(RPFShortcut + "-orfFindErrors.txt", "w")
    if(orfInRange != -1):
        startSiteOffset = orfInRange - charRange
        orfPos = startIndex + startSiteOffset
        #Location is where in genome it is
        foundOrf.write("ATG" + "\t" + str(orfPos) + "\t" + findLoc[2])
	foundOrf.write("\n")
        #Search didn't find an ORF
    else:
        #Create an error log
        for index in xrange(len(findLoc)):
            errorFile.write(findLoc[index] + "\t")
            errorFile.write("\n")
    foundOrf.close()
    errorFile.close()


def countNonGenomeChars(genome, endPos):
    charCount = 0
    pos = 0
    isCounting = False
    while(pos < endPos or isCounting == True):
        if(genome[pos] == '>'):
            isCounting = True
        #Count must occur before testing whether end character reached
        #So that the end character and beginning character are both counted
        if(isCounting == True):
            charCount += 1
        if(isCounting == True and genome[pos] == ']'):
            isCounting = False
        pos += 1    
    return charCount

def sortByPercent(pathToTextFile, indexOfSortableQuantity, decimalSortPercent):
    #text file is thing you want sorted
    #index is the thing you want sorted (eg read values) in that file
    with open(pathToTextFile) as sortingFile:
        sortingFullList = []
        for line in sortingFile:
            breakUpLine = line.split()
            lineList = []
            lineList.append(int(breakUpLinep[indexOfSortableQuantity]))
            for index in xrange(len(breakUpLine)):
                if(index != indexOfSortableQuantity):
                    lineList.append(breakUpLine[index])
            sortingFullList.append(lineList)
        percentToSort = decimalSortPercent
        numberOfLists = ceiling(1,percentToSort)
        sortedFullList = sorted(sortingFullList)
        percentLists = open(pathToTextFile + "-sortedbypercent.txt", "w")
        itemsToAdd = ceiling(percentToSort,len(sortedFullList))
        currIndex = len(sortedFullList)
        for x in xrange(numberOfLists):
            newList = []
            currLength = 0
            for value in xrange(currIndex-1,-1,-1):
                if(currLength <= itemsToadd):
                    newList.append(sortedFullList[value])
                    currIndex -= 1
                    currLength ++ 1
                else:
                    break
            percentLists.write(newList)
            percentLists.write("\n")


def ceiling(num1, num2):
    if(int(num1/num2) == num1*1.0/num2):
        return int(num1/num2)
    else:
        return int(num1/num2 + 1)
    
def create_28():     
    with open(RPFShortcut + "-28reads.txt", "w") as fo:
        with open("accepted_hits.sam") as f:
            for line in f:
                noSpace = line.split()
		if(noSpace[0][0] != "@"):
		    if(noSpace[5] == "28M"):
                   	newLine = "\t".join(noSpace)
                    	fo.write(newLine)
			fo.write("\n")
                    
def sam2bam():
    os.system("samtools view -bT "+genome+" "+ RPFShortcut+"-28reads.sam"+" > "+ RPFShortcut+"-28reads.bam")

def every_three():
    start = [ ]
    end = [ ]
    line_count = 0
    
    # Turn Naina's coordinates into two arrays with start and end coordinates
    with open(RPFShortcut + "-blocks.txt") as h:
        for line in h:
            block_coords = line.split()
            start[line_count] = block_coords[0]
            end[line_count] = block_coords[1]
            line_count+=1

    # Get genome from main read in
    genomeArray = open(genome).read()
    isStopCodon=False
    stopCodonList = ["TAG","TAA","TGA"]
    # Read in Jacob's file
    with open(RPFShortcut + "-foundORFs.txt") as g:
        for line in g:
            foundBlock = False
            stopCodonFound = False
            getORFCoord = line.split()
            orfStart = int(getORFCoord[1])
            #Length is number blocks found
            for start_coord in xrange(len(start)):   
                corrStartCoord = start[start_coord] + countNonGenomeChars(genomeArray,start[start_coord])
                corrEndCoord = end[start_coord] + countNonGenomeChars(genomeArray, end[start_coord])       
                if (orfStart >= start[corrStartCoord] and orfStart <= end[corrEndCoord]):
                    orig_start = orfStart
                    foundBlock = True
                    goodPeaks = 0
                    totalTranslocs = 0
                    while (orfStart + 3 < end[start_coord] and stopCodonFound == False):
                        #Find temp_start # in Alex's index
                        offset1 = orfStart+1
                        offset2 = orfStart+2
                        offset3 = orfStart+3
                        alexReads = open(RPFShortcut + "-readCounts.txt")
                        totalTranslocs += 1
                        if (alexReads[offset3] > alexReads[offset2] 
                                and alexReads[offset3] > alexReads[offset1]):
                            goodPeaks += 1
                        if (genomeArray[offset1:offset1+3] in stopCodonList):
                            stopCodonFound = True
                            goodStopCodons = open(RPFShortcut + "-stopcodons.txt", "w")
                            goodStopCodons.write(orig_start + "\t" + offset3)
                            goodStopCodons.write("\n")
                            goodStopCodons.close()
                            alexReads.close()
                        orfStart += 3
                    if(not stopCodonFound):
                        badPeakFile = open(RPFShortcut + "ORF_SCErr.txt", 'w')
                        badPeakFile.write(orig_start + "\t")
                        badPeakFile.write("Stop codon not found within bounds\n")
                        badPeakFile.close()
                if(not foundBlock):
                    badPeakFile = open(RPFShortcut + "ORF_SCErr.txt",'w')
                    badPeakFile.write(orig_start + "\t")
                    badPeakFile.write("Stop codon did not match cont read\n")
                    badPeakFile.close()
    genomeArray.close()

def main():
    run_bowtie(genome)
    run_tophat(RNA, RNAShortcut)
    run_tophat(RPF, RPFShortcut)
    findRegions()
    create_blocks()
    find_counts()
    findPeaks(2.5)
    findORFs(15)
    create_28()
    sam2bam()
    new_counts()
    every_three()

main()
