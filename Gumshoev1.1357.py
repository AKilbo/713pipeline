#!/usr/bin/env python

import sys
import os
import csv
import time
from Bio.Seq import Seq
from operator import itemgetter

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

#From http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary
def createOutputDir():
    try:
	os.stat(outputDir)
    except:
        os.mkdir(outputDir)

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

def moveHitsFiles():
    os.chdir(RPFShortcut + "_thout")
    os.system("mv " + "accepted_hits.sam " + os.path.dirname(os.getcwd()))
    os.system("mv " + "accepted_hits.bam " + os.path.dirname(os.getcwd()))
    changeToHomeDirec()
    os.system("mv " + "accepted_hits.sam "+os.path.join(os.getcwd(),outputDir))
    os.system("mv " + "accepted_hits.bam "+os.path.join(os.getcwd(),outputDir))
    print "MOVED IT"

def sortIntoStrands():
    print("Sorting into strands..")
    print os.getcwd()
    with open(RPFShortcut + "-posStrandBlocks.sam", "w") as posStrand:
        with open(RPFShortcut + "-negStrandBlocks.sam", "w") as negStrand:
            with open("accepted_hits.sam") as profReadsFile:
                for line in profReadsFile:
                    readSplit = line.split()
                    if(readSplit[0][0] != "@"):
			flagInt = int(readSplit[1])
                        isNegStrand = (flagInt>>4)&1
                        if(isNegStrand == 1):
			    negStrand.write("\t".join(readSplit))
			    negStrand.write("\n")
			elif(isNegStrand == 0):
			    posStrand.write("\t".join(readSplit))
			    posStrand.write("\n")
			else:
			    with open(RPFShortcut + "-strandErrors.txt", "w") as strandErrors:
  				strandErrors.write("\t".join(readSplit))
				print "I got 99 problems and a strand error is one"
			    strandErrors.close()
		posStrand.close()
		negStrand.close()
    print os.getcwd()

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
    #with open(RPFShortcut + "Blocks.txt", "r") as f:
        #with open(RPFShortcut + "-blocks.txt", 'a') as blocks:
            #for line in f:
                #if line.strip():
                    #blocks.write("\t".join(line.split()[0:5]) + "\n")
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

def createPosBlocks():
    print "Creating positive strand blocks.."
    os.system("mv " + RPFShortcut+ "-posStrandBlocks.sam " + os.path.dirname(os.getcwd()))
    changeToHomeDirec()
    os.system("samtools view -bT Spar_ultrascaf.fa " + RPFShortcut + "-posStrandBlocks.sam > " + RPFShortcut + "-posStrandBlocks.bam")
    os.system("mv " + RPFShortcut +  "-posStrandBlocks.sam "+os.path.join(os.getcwd(),outputDir))
    os.system("mv " + RPFShortcut + "-posStrandBlocks.bam "+os.path.join(os.getcwd(),outputDir))
    print "Done!"
    
# Find read counts for each index in genome given
# Written by Alex/Shalyn
def new_counts():
    print("Finding read counts...")
    os.system("samtools sort " + RPFShortcut + '-28reads.bam' + RPFShortcut + "-28reads.sorted")
    os.system("bedtools genomecov -d -ibam "+RPFShortcut+"-28reads.sorted.bam "+ 
            "-g " + genomeShortcut +"-regionResults.txt > "+ RPFShortcut + "-readCounts.txt")
    print("Done!")
    # Now have readCounts.txt file in directory

def findALLtheorfs():
    timeStart = time.time()
    changeToHomeDirec()
    os.system("cp " + genome + " " + os.path.join(os.getcwd(),outputDir))
    changeToOutputDirec()
    #Blocks
    filterOutSmallReads(RPFShortcut + "Blocks.txt")
    #Blocks with length >60 nt
    naina = open(RPFShortcut+"Blocks-filtered.txt")
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
    for line in naina:
        startDist = 0
        startsFound = 0
        nainaSplit = line.split()
        if(nainaSplit[3] == "+"):
            lineCount += 1
	    #Corrects base for chromosome and headers
            start=orgGenome.find(nainaSplit[0])+int(nainaSplit[1])+1+findRightOffset(offsetDict,int(nainaSplit[1]))
            end=orgGenome.find(nainaSplit[0])+int(nainaSplit[2])+findRightOffset(offsetDict,int(nainaSplit[2]))
            genomeToSearch = str(orgGenome[start:end])
            origSearch = str(orgGenome[start:end])
            genomeLen = len(origSearch)
            if(nainaSplit[3] == "-"):
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
                        if(nainaSplit[3] == "+"):
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
                        orfCoords.append(nainaSplit[3])
                        orfCoords.append(nainaSplit[0])
                        removeHeaderStart = int(nainaSplit[1]) + startDist+(startsFound-1)
                        orfCoords.append(removeHeaderStart)
                        orfCoords.append(removeHeaderStart + orfLen-1)
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

def findRightOffset(offsetDict, query):
    min = -1
    for key in offsetDict:
        if(query >= key and key > min):
            min = key
    return offsetDict[min]
    

def ceiling(num1, num2):
    if(int(num1/num2) == num1*1.0/num2):
        return int(num1/num2)
    else:
        return int(num1/num2 + 1)

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

def filterOrfs():
    print os.getcwd()
    moreGoddamnFiltering = open(RPFShortcut + "-orfs-filtered.txt", "w")
    IHateAllOrfsWithAFieryBurningPassion = open(RPFShortcut + "-orfs.txt")
    sillyGenome = open(genome).read()
    stopCodonList = ["TAA","TAG","TGA"]
    lineCount = 0
    goodEnds = 0
    goodDiff = 0
    goodOrfs = 0
    goodNeg = 0
    goodPos = 0
    negLine = 0
    posLine = 0
    for freakinLine in IHateAllOrfsWithAFieryBurningPassion:
        lineCount += 1
        splittingIsTheDevil = freakinLine.split()
        if(splittingIsTheDevil[5] == "-"):
            negLine += 1
        else:
            posLine += 1
        start = int(splittingIsTheDevil[1])
        end = int(splittingIsTheDevil[2])
        blockLen = abs(int(splittingIsTheDevil[4]) - int(splittingIsTheDevil[3]) + 1)
        orfLen = splittingIsTheDevil[0]
        difference = blockLen - int(orfLen)
        freakinGenome = sillyGenome[start:end+1]
        if(splittingIsTheDevil[5] == "-"):
            print sillyGenome[end:start]
        if(splittingIsTheDevil[5] == "-"):
            freakinGenome = Seq(sillyGenome[end:start]).reverse_complement()
        if(str(freakinGenome[0:3]) == "ATG" and
           str(freakinGenome[len(freakinGenome)-3:len(freakinGenome)]) in stopCodonList):
           goodEnds += 1
           if(splittingIsTheDevil[5] == "-"):
               goodNeg += 1
           else:
               goodPos += 1
           if(difference > 25 and difference < 40):
               goodDiff += 1
               moreGoddamnFiltering.write("\t".join(splittingIsTheDevil))
               moreGoddamnFiltering.write("\n")
               goodOrfs += 1
    moreGoddamnFiltering.close()
    IHateAllOrfsWithAFieryBurningPassion.close()

#function that takes in jacobs orf file and does some simple math to report where the untranslated regions of the blocks are
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
        lineList.append(int(splitLine[1]))
        lineList.append(int(splitLine[2]))
        exonList.append(lineList)
    sortedExons = sorted(exonList, key=itemgetter(0,2))
    origChrom = sortedExons[0][1]
    for index in xrange(len(sortedExons)):
        lineToWrite = sortedExons[index]
        currChrom = lineToWrite[1]
        if(currChrom != origChrom):
            exonOutput.write("\n")
            origChrom = currChrom
      #  exonOutput.write(">")
        exonOutput.write(lineToWrite[1] + " ")
	exonOutput.write(str(lineToWrite[2]) + " ")
        exonOutput.write(str(lineToWrite[3]))
        exonOutput.write("\n")
    for index in xrange(len(sortedExons)):
        currLine = sortedExons[index]
        seqOutput.write(">")
        seqOutput.write(currLine[1] + "\n")
        seqOutput.write(genomeString[currLine[4]:currLine[5]+1])
        seqOutput.write("\n")
    seqOutput.close()
    exonOutput.close()
    jacobFile.close()

def glimmerRun(exons,seqs):
    #integrates GlimmerHMM with the genome file and the exons file 
    #to produce the .gff file in the output dir
    print os.getcwd()
    os.system("~/GlimmerHMM/train/make")
    os.system("~/GlimmerHMM/train/trainGlimmerHMM " + seqs + " "+ exons);

#Produce file exons.txt that give exon coordinates for glimmer
def getChromLengths():
    names = []
    length = []
    lastPos = 0

    # Create lists of chromosome names and their length
    with open(genome) as getHeaders:
        for line in getHeaders:
            if line.startswith(">"):
                sepLine = line.split()
                temp = list(sepLine[0])
                names.append("".join(temp[6:]))
                tempLen = list(sepLine[1])
                length.append("".join(tempLen[0:len(tempLen)-2]))
    #print names
    #print length

    count = 0
    currPos = int(length[count])
    #seen= []
    #totalSum = 0
    #for i in xrange(len(length)):
    #    totalSum+=int(length[i])
    
    #Create Exon file for Glimmer
    with open(RPFShortcut+"-orfSortedIncr.txt") as ORFS:
        for line in ORFS:
            sepLine = line.split()
            start = int(sepLine[1])
            end = int(sepLine[2])
            #strand = sepLine[5]
            #seen.append(start)
            chrom=sepLine[6]

            # while current start coordinate is greater than our
            # current location in the genome, add the next chromosome to currPos
            exons = open("exons.txt","w")
            while count<len(length) and start>currPos:
                lastPos = currPos
                currPos = currPos + int(length[count+1])
                count+=1
                exons.write("\n")

            newStart = start - lastPos
            newEnd = end - lastPos

            if count < len(length):
                if currPos == totalSum and seen[0]==start:
                    exons.close()
                    ORFS.close()
                if newStart>0:
                    exons.write(names[count]+" "+str(newStart)+" "+str(newEnd)+"\n")

#Sort the coordinates file in increasing order
def sortCoordFileIncr():
    with open(RPFShortcut+"-orfs-filtered.txt") as coords:
        data = []
        data=coords.readlines()
        #for line in coords:
        #    line=line.split()
        #    data.append(line)

    data.sort(key=lambda s:(int(s[1])))
    data.sort(key=lambda s:(int(s[6])))

    with open (RPFShortcut+"-orfSortedIncr.txt","a") as output:
        for line in data:
            output.write(line)
            #output.write("\t".join(line))
            #output.write("\n")

def main():
    createOutputDir()
    print "DIRECTORY was CHANGED"
    #run_bowtie(genome)
    print "BOWTIE is DONE"
    run_tophat(RNA, RNAShortcut)
    run_tophat(RPF, RPFShortcut)
    moveHitsFiles()
    print "TOPHAT DONE"
    findRegions()
    create_blocks()
    #find_counts()
    findALLtheorfs()
    filterOrfs()
    #findutrs()
    makeExonSeqFile()
    #sortCoordFileIncr()
    #getChromLengths()
    exons = "~/"+outputDir+"/"+RPFShortcut+"-orfs-exons.fasta"
    seqs = "~/"+outputDir+"/"+RPFShortcut+"-orfs-seqs.fasta"
    glimmerRun(exons,seqs)
main()
