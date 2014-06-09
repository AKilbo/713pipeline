#!/usr/bin/env python

import sys
import os
import csv
import time
from Bio.Seq import Seq
from operator import itemgetter

# Take in genome.fa, RNA_reads.fastq, RPF_reads.fastq, and output
# directory as arguments in command line.
genome = sys.argv[1]
genomeShortcut = os.path.splitext(os.path.basename(genome))[0]

RNA = sys.argv[2]
RNAShortcut = os.path.splitext(os.path.basename(RNA))[0]
                               
RPF = sys.argv[3]
RPFShortcut = os.path.splitext(os.path.basename(RPF))[0]

outputDir = sys.argv[4]
outputShortcut = os.path.splitext(os.path.basename(outputDir))[0]

#Use from Tophat directories or output directory
def changeToHomeDirec():
    os.chdir(os.path.dirname(os.getcwd()))

#Use from home
def changeToOutputDirec():
    os.chdir(os.path.join(os.getcwd(),outputDir))

#Use from home
def changeToTHRNA():
    os.chdir(os.path.join(os.getcwd(),RNAShortcut + "_thout"))

#Use from home
def changeToTHRPF():
    os.chdir(os.path.join(os.getcwd(),RPFShortcut + "_thout"))

#From http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary
def createOutputDir():
    try:
        changeToOutputDirec()
        changeToHomeDirec()
    except:
        os.mkdir(outputDir)



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
    os.chdir(name + "_thout")
    os.system("samtools view -h -o accepted_hits.sam accepted_hits.bam")
    changeToHomeDirec()
    print("Done!\n")

# Run Tophat on given RNA or RPF data, output to _thout directory
# Written by Shalyn
def run_tophat(data,name):
    print ("Running tophat on " + name+ "...")
    os.system("tophat -o "+name+"_thout INDEX "+data)
    convert(name)
    print("Done!\n")

#moves the tophat output into a separate directory
#Written by Shalyn
def moveHitsFiles():
    os.chdir(RPFShortcut + "_thout")
    os.system("mv " + "accepted_hits.sam " + os.path.dirname(os.getcwd()))
    os.system("mv " + "accepted_hits.bam " + os.path.dirname(os.getcwd()))
    changeToHomeDirec()
    os.system("mv " + "accepted_hits.sam "+os.path.join(os.getcwd(),outputDir))
    os.system("mv " + "accepted_hits.bam "+os.path.join(os.getcwd(),outputDir))
    print "MOVED IT"

#finding regions from the sam file from tophat
#Written by Jacob
def findRegions():
    print os.getcwd()
    changeToOutputDirec()
    regionFile = open(genomeShortcut + "-regionResults.txt", "w")
    with open("accepted_hits.sam") as regionFileName:
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
    print os.getcwd()

# Create file to write start/stop coordinates to of
# regions with continuous coverage
# Written by Naina
def create_blocks():
    print ("Creating coordinates of areas with continuous reads...")
    os.system("bedtools bamtobed -i accepted_hits.bam > reads.bed")
    os.system("sort -k1,1 -k2,2n reads.bed > reads.sort.bed")
    os.system("bedtools merge -i reads.sort.bed -s > "+
              RPFShortcut + "Blocks.txt")
    print("Done!\n")

#Finding ORFs from the blocks 
#Written by Jacob   
def findALLtheorfs():
    timeStart = time.time()
    changeToHomeDirec()
    os.system("cp " + genome + " " + os.path.join(os.getcwd(),outputDir))
    changeToOutputDirec()
    #Blocks
    filterOutSmallReads(RPFShortcut + "Blocks.txt")
    #Blocks with length >60 nt
    blocks = open(RPFShortcut+"Blocks-filtered.txt")
    #Reads at each base
    print "made it this far"
    orgGenome = open(genome).read()
    foundORFs = open(RPFShortcut + "-orfs.txt", "w")
    #Offset dictionary for a given abse
    offsetDict = getOffsets(orgGenome, genomeShortcut + "-regionResults.txt")
    print offsetDict
    print "great"
    stopCodonList = ["TAA","TAG","TGA"]
    lineCount = 0
    found = 0
    for line in blocks:
        startDist = 0
        startsFound = 0
        blocksFiltered = line.split()
        if(blocksFiltered[3] == "+"):
            lineCount += 1
	    #Corrects base for chromosome and headers
            start=orgGenome.find(blocksFiltered[0])+int(blocksFiltered[1])+1+findRightOffset(offsetDict,int(blocksFiltered[1]))
            end=orgGenome.find(blocksFiltered[0])+int(blocksFiltered[2])+findRightOffset(offsetDict,int(blocksFiltered[2]))
            genomeToSearch = str(orgGenome[start:end])
            origSearch = str(orgGenome[start:end])
            genomeLen = len(origSearch)
            if(blocksFiltered[3] == "-"):
                genomeToSearch = str(Seq(genomeToSearch).reverse_complement())
            ORFlist = []
            while(genomeToSearch.find("ATG") != -1):
                subLen = len(genomeToSearch)
                startsFound += 1
                startDist += genomeToSearch.find("ATG")
                startPos = genomeToSearch.find("ATG")
                currPos = startPos + 3
                orfCoords = []
                while(currPos + 2 <= end):
                    if(str(genomeToSearch[currPos:currPos+3]) in stopCodonList):
                        endPos = currPos + 2
                        orfSeq = str(genomeToSearch[startPos:endPos+1])
                        orfLen = len(orfSeq)
                        if(blocksFiltered[3] == "+"):
                            actStart = startDist+(startsFound-1)+start
                            actEnd = actStart + orfLen - 1
                        else:
                            actStart = start + genomeLen -1 - startDist+(startsFound-1)
                            actEnd = actStart - orfLen
                        orfCoords.append(orfLen)
                        orfCoords.append(str(actStart))
                        orfCoords.append(str(actEnd))
                        orfCoords.append(str(start))
                        orfCoords.append(str(end))
                        orfCoords.append(blocksFiltered[3])
                        orfCoords.append(blocksFiltered[0])
                        removeHeaderStart = int(blocksFiltered[1]) + startDist+(startsFound-1)
                        orfCoords.append(removeHeaderStart)
                        orfCoords.append(removeHeaderStart + orfLen-1)
                        orfCoords.append(int(blocksFiltered[1]))
                        orfCoords.append(int(blocksFiltered[2]))
                        found += 1
			break
                    else:
                        currPos += 3
                if(len(orfCoords) > 0):
                    ORFlist.append(orfCoords)
                genomeToSearch = genomeToSearch[startPos+1:end]
            if(len(ORFlist) > 0):
                bestOrf = max(ORFlist)
                for index in xrange(len(bestOrf)):
                    foundORFs.write(str(bestOrf[index]) + "\t")
                foundORFs.write("\n")
	    print found
	    print lineCount
    foundORFs.close()
    naina.close()
    print time.time() - timeStart

#filters small reads
#Written by Jacob
def filterOutSmallReads(filePath):
    fileToFilter = open(filePath)
    filteringFile = open(os.path.splitext(os.path.basename(filePath))[0] + "-filtered.txt", "w")
    for line in fileToFilter:
        breakLineUp = line.split()
        if(int(breakLineUp[2])-int(breakLineUp[1]) > 60):
            filteringFile.write("\t".join(breakLineUp))
            filteringFile.write("\n")
    fileToFilter.close()
    filteringFile.close()

#finds the right offset
#Written by Jacob
def findRightOffset(offsetDict, query):
    min = -1
    for key in offsetDict:
        if(query >= key and key > min):
            min = key
    return offsetDict[min]
    
# rounds it up
#Written by Jacob
def ceiling(num1, num2):
    if(int(num1/num2) == num1*1.0/num2):
        return int(num1/num2)
    else:
        return int(num1/num2 + 1)

# calculates and returns the offsets
#Written by Jacob
def getOffsets(genome, regionFilePath):
    regionDict = dict()
    regionLists = open(regionFilePath)
    for line in regionLists:
        regionSplit = line.split()
        regionName = regionSplit[0]
        regionLoc = genome.find(regionName)
        if(regionLoc != -1):
            while(genome[regionLoc] != "]"):
                regionLoc += 1
            firstNonHeaderPos = regionLoc + 2
            offset = countNonGenomeChars(genome,firstNonHeaderPos)
        regionDict[firstNonHeaderPos] = offset
    regionLists.close()
    return regionDict

# counts the non genome characters
#Written by Jacob
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

#Filters the ORFs
#Written by Jacob
def filterOrfs():
    filteredOrfs = open(RPFShortcut + "-orfs-filtered.txt", "w")
    ORFs = open(RPFShortcut + "-orfs.txt")
    genome = open(genome).read()
    stopCodonList = ["TAA","TAG","TGA"]
    for line in ORFs:
        splits = line.split()
        start = int(splits[1])
        end = int(splits[2])
        blockLen = abs(int(splits[4]) - int(splits[3]) + 1)
        orfLen = splits[0]
        difference = blockLen - int(orfLen)
        gen = genome[start:end+1]
        if(splits[5] == "-"):
            gen = Seq(genome[end:start]).reverse_complement()
        if(str(gen[0:3]) == "ATG" and
           str(gen[len(gen)-3:len(gen)]) in stopCodonList):
           if(difference > 25 and difference < 40):
               filteredOrfs.write("\t".join(splits))
               filteredOrfs.write("\n")
    filteredOrfs.close()
    ORFs.close()

#function that takes in jacobs orf file and 
#does some simple math to report where the untranslated regions of the blocks are
def findutrs():
    print os.getcwd()
    utrdoc = open('untranslatedregions.txt','w')
    with open(RPFShortcut+"-orfs-filtered.txt") as getHeaders:
        for line in getHeaders:
            (length, startorf, endorf, startblock,endblock,strand,geneome) = line.split()
            utrstart = startblock
            utrend = startorf
            secondutr = endorf
            secondutrend = endblock
            utrdoc.write(str(utrstart) +"\t" + str(utrend) + "\t" + strand + "\t" + geneome +"\n" )
            utrdoc.write(str(secondutr) +"\t" + str(secondutrend) + "\t" + strand + "\t" + geneome +"\n" )
    utrdoc.close()


#Written by Jacob
def makeExonSeqFile():
    jacobFile = open(RPFShortcut + "-orfs-filtered.txt")
    genomeString = open(genome).read()
    exonOutput = open(RPFShortcut + "-orfs-exons.fasta", "w")
    exonList = []
    seqOutput = open(RPFShortcut + "-orfs-seqs.fasta", "w")
    for line in jacobFile:
        lineList = []
        splitLine = line.split()
        lineList.append(genomeString.find(splitLine[6]))
        lineList.append(splitLine[6])
        lineList.append(int(splitLine[7]))
        lineList.append(int(splitLine[8]))
        lineList.append(int(splitLine[9]))
        lineList.append(int(splitLine[10]))
        exonList.append(lineList)
    sortedExons = sorted(exonList, key=itemgetter(0,2))
    origChrom = sortedExons[0][1]
    for index in xrange(len(sortedExons)):
        lineToWrite = sortedExons[index]
        currChrom = lineToWrite[1]
        #if(currChrom != origChrom):
       # exonOutput.write("\n")
       #     origChrom = currChrom
      #  exonOutput.write(">")
        exonOutput.write(lineToWrite[1] +"_"+ str(index) + " ")
	exonOutput.write(str(1) + " ")
        exonOutput.write(str(lineToWrite[3]-lineToWrite[2]))
        exonOutput.write("\n")
	exonOutput.write("\n")
    for index in xrange(len(sortedExons)):
        currLine = sortedExons[index]
        seqOutput.write(">")
        seqOutput.write(currLine[1]+ "_" + str(index) + "\n")
        seqOutput.write(genomeString[currLine[2]:currLine[3]])
        seqOutput.write("\n")
    seqOutput.close()
    exonOutput.close()
    jacobFile.close()


#integrates GlimmerHMM with the genome file and the exons file 
#to produce the .gff file in the output dir
#Written by Naina
def glimmerRun(exons,seqs):
    print os.getcwd()
    os.system("~/GlimmerHMM/train/make")
    os.system("~/GlimmerHMM/train/trainGlimmerHMM " + seqs + " "+ exons);

def main():    
    createOutputDir()
    print "DIRECTORY was CHANGED"
    run_bowtie(genome)
    print "BOWTIE is DONE"
    run_tophat(RNA, RNAShortcut)
    run_tophat(RPF, RPFShortcut)
    moveHitsFiles()
    print "TOPHAT DONE"
    findRegions()
    create_blocks()
    findALLtheorfs()
    filterOrfs()
    makeExonSeqFile()
    exons = "~/"+outputDir+"/"+RPFShortcut+"-orfs-exons.fasta"
    seqs = "~/"+outputDir+"/"+RPFShortcut+"-orfs-seqs.fasta"
    glimmerRun(exons,seqs)

main()
