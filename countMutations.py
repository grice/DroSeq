# Modified from the SHAPE-MaP pipeline to call mutations
# in miRNA cutting assay data. Functions were kepts and
# main loop was wrapped in an if __name__ == "__main__"
# statement
#
# Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
# Counts up mutations from sequencing data.
# Copyright Steven Busan 2014

#---------------------------------------------------------------------------------------------------
# GPL statement:
#
# This file is part of Shapemapper.
#
# ShapeMapper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ShapeMapper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ShapeMapper.  If not, see <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------------------

import sys, os, traceback, copy, re, time

#import conf
import numpy
import matplotlib.pyplot as plot

# def safeDivide(num,den):
# # convert input numbers to floats and divide
# # If denominator is 0, output -999 to indicate error.
#     result = -999.0
#     try:
#         result = float(num)/float(den)
#     except ZeroDivisionError:
#         pass
#     return result
# 
# def parseFasta(fastaFile):
# # given an open fasta file, return a dict of sequence names
# # and associated sequence strings
#     seq = {}
#     lines = fastaFile.readlines()
#     currentSeqName = ""
#     for line in lines:
#         if line[0] == ">":
#             # sequence name found
#             currentSeqName = line.strip()[1:]
#             seq[currentSeqName] = ""
#         else:
#             # remove whitespace and append to current sequence
#             seq[currentSeqName] += "".join(line.split())
#     return seq
# 
# # 
def parseCigarString(cigarString, read, quals, refSeq, startIndex):
    # define regular expression patterns for parsing CIGAR string
    item = re.compile(r"[0-9]+[M|I|D|S]")
    number = re.compile(r"([0-9]+)")
    splitCigarString = item.findall(cigarString.strip())
    cigarList = [number.split(x)[1:] for x in splitCigarString]
    #print "CIGAR: %s\nsplit CIGAR: %s"%(cigarString,str(cigarList))
    readCopy = read        
    alignedQuals = {}
    lastQual= "#"
    refIndex = startIndex # begin at leftmost matching position
    events = {} # keys will be nuc indices with respect to reference sequence
    #inserts = {}
    for i in range(len(cigarList)):
        region = cigarList[i]
        regionLength = int(region[0])
        regionType = region[1]
        # alignment match (could be match or mismatch, use raw read to resolve)
        if regionType == "M":
            matchRegion = read[:regionLength]
            qualsRegion = quals[:regionLength]
            read = read[regionLength:]
            quals = quals[regionLength:]
            for regionIndex in range(len(matchRegion)):
                nuc = matchRegion[regionIndex]
                qual = qualsRegion[regionIndex]
                lastQual = qual
                alignedQuals[refIndex] = qual

                # catch insertion at end of reference sequence
                if refIndex > len(refSeq) -1:
                    continue
                elif nuc != refSeq[refIndex]:
                    # nuc mismatch found
                    #print "Mismatch at ref index %i"%refIndex
                    events[refIndex] = nuc
                else:
                    events[refIndex] = "|"
                refIndex += 1
        # insertion
        elif regionType == "I":
            insertedSeq = read[:regionLength]
            #inserts[refIndex] = insertedSeq
            read = read[regionLength:]
            quals = quals[regionLength:]
            # ignore for now (insertions are usually unhelpful)
        # deletion
        elif regionType == "D":
            #print "Deletion %i nucs long at ref index %i"%(regionLength, refIndex)
            for deletionIndex in range(regionLength):
                events[refIndex] = "-"
                alignedQuals[refIndex] = lastQual # don't have a phred score for a deletion (no basecall associated with position), so just use last nearby score
                refIndex += 1
        # padding
        elif regionType == "S":
            read = read[regionLength:]
            quals = quals[regionLength:] # missing from v12b and all pipeline versions before - results in misaligned phred scores and errors combining read pairs
            if i == len(cigarList)-1: # rightmost end of read
                for offsetIndex in range(regionLength):
                    padIndex = refIndex+offsetIndex
                    if padIndex >= 0 and padIndex < len(refSeq):
                    #if True:
                        events[padIndex] = "s"
                        alignedQuals[padIndex] = "#"
            elif i == 0: # leftmost end of read 
                for offsetIndex in range(regionLength+1):
                    padIndex = refIndex-offsetIndex
                    if padIndex > 0 and padIndex < len(refSeq):
                    #if True:
                        events[padIndex] = "s"
                        alignedQuals[padIndex] = "#"
    sortedKeys  = sorted(events.keys())
    printEvents = ""
    printQuals  = ""
    evalName    = ""
    lowQualFlag = False
    mutCount    = 0
    for i in xrange(min(sortedKeys),max(sortedKeys)+1):
        try:
            printEvents += events[i]
            printQuals += alignedQuals[i]

            if events[i] != "|":
                mutCount += 1
                qscore = ord(alignedQuals[i]) -33 
                if qscore < 20:
                    lowQualFlag = True
                evalName += "{0}{1}".format(i, events[i])
        except KeyError:
            printEvents += " "
            printQuals += " "
                
    if mutCount == 0:
        evalName = None
    #print printEvents
    #print printQuals
    #print evalName, lowQualFlag
    return printEvents, printQuals, evalName, mutCount, lowQualFlag

if __name__ == '__main__':
    cigar = "4M1D5M"
    read   = "AGCAAAGGC"
    refSeq = "AGTAGTAGGC"
    quals  = "GGGGGGGGGG"
    start  = 0
    parseCigarString(cigar, read, quals, refSeq, start)
