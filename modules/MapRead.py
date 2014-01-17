#!/usr/bin/env python

import os 
import sys
import re
import difflib
from MapRead import *
from Sequence import *
from Utilities import *
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs


#from ..gtTools.seq import split *


##########################################################################################################################################
#####################################################  MAPPED READ  CLASS START ############################################################
##########################################################################################################################################
##########################################################################################################################################
#####################################################  MAPLINE-FILE CLASS START ##########################################################
########print##################################################################################################################################

class MapReads:

    def __init__(self,mapFile,strandSpecific):
 
        if strandSpecific == "OPPOSITE":    self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): False,("-","-"): False}
        elif strandSpecific == "MATCH":         self.strandKey = {("+","-"): False,("-","+"): False,("+","+"): True,("-","-"): True}
        else:                               self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): True,("-","-"): True}

    
        self.samStrand = {'0':'+','16':'-','+':'0','-':'16'}


        self.fileName = mapFile; self.gapped = False

        self.minOVERLAP=5


        try:                self.handle = open(self.fileName)
        except IOError:     errorQuit('ERROR: File '+self.fileName+' does not exist')
        
        self.open = True

        self.fileOpen = True
        


        #myLine = self.handle.readline().split()
        self.rawLine = self.handle.readline().split()
        myLine=self.rawLine
        self.samHeader=[]
        if len(self.rawLine) == 0:
            self.open = False
            self.format = None
            return
        if self.rawLine[0][0] == "@":
            while self.rawLine[0][0] == "@":
                self.samHeader.append("\t".join(self.rawLine))
                #print "\t".join(self.rawLine)
                self.rawLine = self.handle.readline().split()



        if self.fileName.split(".")[-1] == "mapping":
            self.seqLen = len(self.rawLine[1]) - 1
            self.format = "MAPPING"
            if len(self.rawLine[2].split(":")[0].split("|"))==6:
                self.getNextRead = self.getMappingFeatureRead
                self.printData   = self.printFeatureData
        
        elif self.fileName.split(".")[-1] == "sam":
            self.format = "SAM"
            self.seqLen = len(self.rawLine[9])
            if len(self.rawLine[2].split(":")[0].split("|"))==6:
                self.getNextRead = self.getSamFeatureRead
                self.printData   = self.printFeatureData
            elif self.rawLine[2][0:3]=="chr":
                if len(self.rawLine[5].split("M"))==2 and len(self.samHeader)>1:
                    for h in self.samHeader: print h
                self.getNextRead = self.getSamGenomeRead
                self.printData   = self.printGenomeData
        else:
            errorQuit(".mapping or .sam extension required")


            ### FIX THIS AND MAKE IT MORE CONCISE ###!!!!
            #self.readID,self.readSeq,self.subs,self.qual  = self.rawLine[0],self.rawLine[1],self.rawLine[6],self.rawLine[8]
            
#        if self.format == "SAM":
#            print "YO"
#            cigarSpot = myLine[5].split("N"); self.refData = myLine[2].split("|")
#            if len(cigarSpot) > 1 and len(self.refData)>1:
#                self.style = 'gapFeature'
#                self.seqlen = len(myLine[9]); self.rawLine = myLine; self.next = self.nextGappedLine; self.next()
                #print myLine 

#            elif len(self.refData) == 1:
#                self.rName,self.samStrand,self.fName,self.fPos,self.cigar,self.read,self.qual,self.subs = myLine[0],myLine[1],myLine[2],int(myLine[3]),myLine[5],myLine[9],myLine[10],int(myLine[11].split(":")[-1])
#                self.geneID,self.refStrand,self.chr,self.seqLen = None,"+",self.fName,len(self.read)-1
#                self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
#                self.next = self.nextGenomeLine

     
 #           else:
 #               self.rName,self.samStrand,self.fName,self.fPos,self.blank,self.cigar,self.blank,self.blank,self.blank,self.read,self.qual,self.subs,self.blank = myLine
 #               self.fPos = int(self.fPos)-1; self.subs= int(self.subs.split(":")[-1]); self.refData = self.fName.split("|")
 #               self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]
 #               self.seqLen = len(self.read) - 1 
 #               self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
 #               self.next = self.nextSamLine

        
##############################################################################################################
############################################   READ RELOCATION  ##############################################
##############################################################################################################

    

    def relocate(self,locs,fStart):
        read_remain = self.seqLen
        offset = 0; outLoc = ()
        pDist = 0 
        for pair in locs:
            pDist = pair[1]-pair[0]+1
            if pDist <= fStart:
                fStart -= pDist
            elif pDist > fStart + read_remain:
                return outLoc + (pair[0]+fStart,pair[0]+fStart+read_remain)
            else:
                outLoc += (pair[0]+fStart,pair[1])
                read_remain -= (pDist-fStart) 
                fStart = 0

    def cigarRelocate(self,locs,fStart,cigar):
        ## TEMPORARY CODE ###
        if len(locs)==1:
            cLocs=[[int(c.split("N")[0]),int(c.split("N")[1])] for c in ("0N"+cigar[0:-1]).split("M")]
            return tuple([item for pair in [[(fStart+locs[0][0])+cLocs[i][0]+(sum([sum(cLocs[j]) for j in range(i)]))-1,(fStart+locs[0][0])+(sum([sum(cLocs[j]) for j in range(i+1)]))-2] for i in range(len(cLocs))] for item in pair])
    
    def genomeRelocate(self,fStart,cigar):
        cLocs=[[int(c.split("N")[0]),int(c.split("N")[1])] for c in ("0N"+cigar[0:-1]).split("M")]
        return tuple([item for pair in [[fStart+cLocs[i][0]+(sum([sum(cLocs[j]) for j in range(i)]))-1,fStart+(sum([sum(cLocs[j]) for j in range(i+1)]))-2] for i in range(len(cLocs))] for item in pair])
         

    

    def getSamGenomeRead(self):
        #print "\t".join(self.rawLine)
        if len(self.rawLine)==12: self.readID,self.readSeq,self.subs,self.qual,self.genomeData=self.rawLine[0],self.rawLine[9],int(self.rawLine[11].split(":")[-1]),self.rawLine[10],[]
        elif len(self.rawLine)==13: self.readID,self.readSeq,self.subs,self.qual,self.genomeData=self.rawLine[0],self.rawLine[9],0,self.rawLine[10],[]
        else:
            print "WTF"; sys.exit()
        samData=self.rawLine[0:len(self.rawLine)-1]
        while self.rawLine[0] == self.readID:
            self.genomeData.append((self.rawLine[1],self.rawLine[2],self.rawLine[3],self.rawLine[4],self.rawLine[5],self.genomeRelocate(int(self.rawLine[3]),self.rawLine[5])))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break
    def printGenomeData(self):

        GENE_UNIQUE="0.0,0.0"

        if len(self.genomeData)==1:  GENOME_UNIQUE="1.0,0.0"
        else:                        GENOME_UNIQUE="0.0,"+str(1.0/len(self.genomeData))

        
        for g in self.genomeData:
            if len(g[5])==2:
                CLASS="CL:i:INTERGENIC:REG="+str(g[5][0])+","+str(g[5][1])
            else:
                CLASS="CL:i:INTERGENIC:JXN="+str(g[5][1])+","+str(g[5][2])
            samStartData= [self.readID,g[0],g[1],g[2],g[3],g[4],"*",0,0,self.readSeq,self.qual,'MISMATCHES:i:'+str(self.subs)]
            samStartData.extend(['GT:i:'+GENOME_UNIQUE,'TT:i:'+GENE_UNIQUE,CLASS,'GN:i:None','AN:i:None','FM:i:INTERGENIC','SN:i:True'])
            print "\t".join([str(s) for s in samStartData])


    def getSamFeatureRead(self):
        self.readID,self.readSeq,self.subs,self.qual,self.mapData,self.mapNum  = self.rawLine[0],self.rawLine[9],0,self.rawLine[10],[],1
        while self.rawLine[0] == self.readID:
            feature=self.rawLine[2].split(":")
            data=self.cigarRelocate([[int(x) for x in r.split("-")] for r in  feature[1].split("|")],int(self.rawLine[3]),self.rawLine[5])
            
            self.mapData.append((feature[0].split("|"), self.cigarRelocate([[int(x) for x in r.split("-")] for r in  feature[1].split("|")],int(self.rawLine[3]),self.rawLine[5]),self.samStrand[self.rawLine[1]],True))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break

        if len(self.mapData)>1:
            print "HI"
            sys.exit()
        return 


    def getMappingFeatureRead(self):
        self.readID,self.readSeq,self.subs,self.qual,self.mapData,self.mapNum  = self.rawLine[0],self.rawLine[1],self.rawLine[6],self.rawLine[8],[],1

        while self.rawLine[0] == self.readID:
            #print self.rawLine
            feature=self.rawLine[2].split(":")
            #mapData.append((self.rawLine[2].split(":")[0].split("|"), self.relocate([ [int(x) for x in r.split("-")] for r in self.rawLine[2].split(":")[1].split("|")], int(self.rawLine[3])),self.rawLine[5]))
            tmpLocs=[[int(x) for x in r.split("-")] for r in feature[1].split("|")]
            self.mapData.append((feature[0].split("|"), self.relocate( tmpLocs, int(self.rawLine[3])),self.rawLine[5],self.strandKey[self.rawLine[5],feature[0].split("|")[4]],tmpLocs))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break
        #self.mapType = "UNIQUE"
        if len(self.mapData)>1: self.disambiguateMaps()
       


    def printFeatureData(self):
        
        ### RATE GENE/GENOME UNIQUENESS ###

        if len(set([tuple(m[0][0:3]) for m in self.mapData]))==1: GENE_UNIQUE=str(1.0/len(self.mapData))+",0.0"
        else:                                                     GENE_UNIQUE="0.0"+str(1.0/len(self.mapData))
        
        if len(set([m[1] for m in self.mapData]))==1:             GENOME_UNIQUE=str(1.0)+",0.0"
        else:                                                     GENOME_UNIQUE="0.0,"+str(1.0/len(self.mapData))

        for m in self.mapData:

            ### FIX SMALL SPLICING OFFSETS ###
            if len(m[1])>2:
                if m[1][1]-m[1][0]<self.minOVERLAP:
                    if m[2]=="+":
                        self.readSeq,self.qual=self.readSeq[(m[1][1]-m[1][0])+1::],self.qual[(m[1][1]-m[1][0])+1::]
                    else:
                        self.readSeq,self.qual=self.readSeq[0:len(self.readSeq)-(m[1][1]-m[1][0]+1)],self.qual[0:len(self.readSeq)-(m[1][1]-m[1][0]+1)]
                    m = (m[0],m[1][2::],m[2],m[3],m[4])
                if m[1][-1]-m[1][-2]<self.minOVERLAP:
                    if m[2]=="+":
                        self.readSeq,self.qual=self.readSeq[0:len(self.readSeq)-(m[1][-1]-m[1][-2]+1)],self.qual[0:len(self.readSeq)-(m[1][-1]-m[1][-2]+1)]
                    else:
                        self.readSeq,self.qual=self.readSeq[(m[1][-1]-m[1][-2])+1::],self.qual[(m[1][-1]-m[1][-2])+1::]
                    m=(m[0],m[1][0:len(m[1])-2],m[2],m[3],m[4])
            

            SPLICESTR='None'

            
            if m[0][5]=="EXON":
                m[0][5]="EXON:COORDS="+"-".join([str(s) for s in m[4][0]])
            elif m[0][5]=="ITRN":
                if m[1][0] >= m[4][1][0] and m[1][1] <= m[4][1][1]: m[0][5]="ITRN:COORDS="+"-".join([str(s) for s in m[4][1]])
                elif m[1][0] >= m[4][1][0]: m[0][5]="IJXN:EX-INT="+"-".join([str(s) for s in m[1]])
                elif m[1][1] <= m[4][0][1]: m[0][5]="IJXN:INT-EX="+"-".join([str(s) for s in m[1]])
                elif m[1][0] <= m[4][1][0] and m[1][1] >= m[4][1][0]:   m[0][5]="IJXN:EX-INT="+"-".join([str(m[4][0][1]),str(m[4][1][0])])
                elif m[1][0] <= m[4][1][1] and m[1][1] >= m[4][2][0]:   m[0][5]="IJXN:INT-EX="+"-".join([str(m[4][1][1]),str(m[4][2][0])])
                else:
                    print "WTF"
                    print m
                    print ""
                    sys.exit()

            elif len(m[1])>2:   m[0][5]=m[0][5]+":JXN="+"-".join([str(m[1][i-1])+','+str(m[1][i]) for i in range(2,len(m[1])-1,2)])
                

            elif len(m[1])==2 and (m[0][5] == "KJXN" or m[0][5] =="NJXN"):  
                if m[1][1]==m[4][0][1]:     m[0][5]="EXON:END="+str(m[1][1])
                elif m[1][0]==m[4][1][0]:   m[0][5]="EXON:START="+str(m[1][0])
                else:
                    print "WTF EXON HUH"
                    sys.exit()
            elif m[0][5] != "FILTER":
                print "UNKNOWN THING"
                sys.exit()


            
            cigar="".join([str(m[1][i]-m[1][i-1]+1)+"M" if i%2==1 else str(m[1][i]-m[1][i-1]-1)+"N" for i in xrange(1,len(m[1]))])
            
            samStartData= [self.readID,self.samStrand[m[2]],m[0][3],m[1][0],'255',cigar,"*",0,0,self.readSeq,self.qual,'MISMATCHES:i:'+str(self.subs)]

 

            #samStartData.extend(['GT:i:'+GENOME_UNIQUE,'SP:i:'+SPLICESTR,'TT:i:'+GENE_UNIQUE,'CL:i:'+m[0][5],'GN:i:'+m[0][0],'AN:i:'+m[0][1],'FM:i:'+m[0][2],'SN:i:'+str(m[3])])
            samStartData.extend(['GT:i:'+GENOME_UNIQUE,'TT:i:'+GENE_UNIQUE,'CL:i:'+m[0][5],'GN:i:'+m[0][0],'AN:i:'+m[0][1],'FM:i:'+m[0][2],'SN:i:'+str(m[3])])
            
            print "\t".join([str(s) for s in samStartData])



                



    def disambiguateMaps(self):

        
        ### 0) REMOVE DUPLICATES
        self.mapData.sort();    self.mapData = self.mapData[0:1]+[self.mapData[i] for i in range(1,len(self.mapData)) if self.mapData[i]!=self.mapData[i-1]];
        if len(self.mapData)==1: return
        ### 1) CHECK FOR MULTI SENSE ### 
        if len(set([m[3] for m in self.mapData]))>1:   self.mapData=[m for m in self.mapData if m[3]]
        ### 2) Prioritize Map Codes ###
        mapCodes = [m[0][-1] for m in self.mapData]
        if "FILTER" in mapCodes:
            self.mapData=[m for m in self.mapData if m[0][-1]=="FILTER"]
            if len(self.mapData)>1:
                tmpData=[m for m in self.mapData if m[0][1] != "CHR"]
                if len(tmpData)>1:  tmpData=[t for t in tmpData if t[0][1][0:6] != "MT-RNR"]
                if len(tmpData)>0: self.mapData=tmpData

        elif "EXON" in mapCodes: self.mapData=[m for m in self.mapData if m[0][-1]=="EXON"]
        elif "KJXN" in mapCodes: self.mapData=[m for m in self.mapData if m[0][-1]=="KJXN"]
        if len(self.mapData)==1: return 
        ### 3) Merge Multi-Genes ###
        if len(set([m[1] for m in self.mapData]))==1 and len(set([m[0][3] for m in self.mapData]))==1:
            self.mapData=[([",".join(list(set([self.mapData[i][0][j] for i in range(len(self.mapData))]))) for j in range(len(self.mapData[0][0]))], self.mapData[0][1],self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
            if len(self.mapData)==1: return

        ### 4) Merge Mutlti-Spots --- ###

        if len(set([tuple(m[0]) for m in self.mapData])) == 1:
            if len(set([(m[1][1:len(m[1])-1]) for m in self.mapData])) == 1:
                tmp_starts=(max([m[1][0] for m in self.mapData]),min([m[1][0] for m in self.mapData])); tmp_ends=(max([m[1][-1] for m in self.mapData]),min([m[1][-1] for m in self.mapData]))
                if tmp_starts[0]-tmp_starts[1] == tmp_ends[0]-tmp_ends[1]:
                    tmp_offset = tmp_starts[0]-tmp_starts[1]
                    if tmp_offset*2 < len(self.readSeq)/2.0:
                        self.mapData = [(self.mapData[0][0],(tmp_starts[0],)+self.mapData[0][1][1:len(self.mapData[0][1])-1]+(tmp_ends[1],),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        self.qual = self.readSeq[tmp_offset:len(self.readSeq)-tmp_offset],self.qual[tmp_offset:len(self.readSeq)-tmp_offset]
            else:
                if len(set([m[1][-1] for m in self.mapData])) == 1:
                    tmpLocs=[]; k=-1;
                    while len(set([m[1][k] for m in self.mapData]))==1:
                        tmpLocs.append(self.mapData[0][1][k])
                        k-=1
                    if len(tmpLocs)%2==1:   tmpLocs.append(max([m[1][k] for m in self.mapData]))
                    
                    tmpLocs.reverse(); trimDist=len(self.readSeq)-sum([(tmpLocs[i]-tmpLocs[i-1])+1 for i in range(1,len(tmpLocs),2)])
                    if trimDist < len(self.readSeq)/2.0:
                        self.mapData=[(self.mapData[0][0],tuple(tmpLocs),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        if self.mapData[0][2] == "-": self.readSeq=self.readSeq[0:len(self.readSeq)-trimDist]; self.qual=self.qual[0:len(self.qual)-trimDist]
                        else:                         self.readSeq=self.readSeq[trimDist::]; self.qual=self.qual[trimDist::]
                        return 
                if len(set([m[1][0] for m in self.mapData])) == 1:
                    tmpLocs = []; k=0; 
                    while len(set([m[1][k] for m in self.mapData]))==1:
                        tmpLocs.append(self.mapData[0][1][k])
                        k+=1
                    if len(tmpLocs)%2==1:   tmpLocs.append(max([m[1][k] for m in self.mapData]))
                    
                    trimDist=len(self.readSeq)-sum([(tmpLocs[i]-tmpLocs[i-1])+1 for i in range(1,len(tmpLocs),2)])

                    if trimDist < len(self.readSeq)/2.0:
                        self.mapData=[(self.mapData[0][0],tuple(tmpLocs),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        if self.mapData[0][2] == "+": self.readSeq=self.readSeq[0:len(self.readSeq)-trimDist]; self.qual=self.qual[0:len(self.qual)-trimDist]
                        else:                         self.readSeq=self.readSeq[trimDist::]; self.qual=self.qual[trimDist::]
                        return
            
                        
                



'''

    def nextMappingLine(self):
        
        try: 
            self.rName,self.seq,refFeature,fPos,ref,strand,self.subs,self.locs,self.qual = self.rawLine
           
            ## TO BE REMOVED ##
        #    if refFeature.split("|")[0] == "FILTER":
        #        tmpData=refFeature.split("|")
        #        if tmpData[1] == "RIBO":
        #            refFeature = "NR_RIBO|"+tmpData[2]+"|ribosome|chrR|+|FILTER:0-"+str(int(fPos)+200)
        #        elif tmpData[1] == "MITO":
        #            refFeature = "NR_MITO|"+tmpData[2]+"|mitochondria|chrM|+|FILTER:0-"+str(int(fPos)+200)
            
            refFeature = refFeature.split(":")
            geneID,geneAlt,geneGroup,refChr,refStrand,refType = refFeature[0].split("|")
            #self.loc =self.relocate([ [int(r.split('-')[0]),int(r.split('-')[1])] for r in refFeature[1].split("|")],int(fPos))
            
            
            self.geneLoc = (geneAlt,geneID,geneGroup,refChr,refStrand,refType)
            self.refLoc  = (refChr,strand,self.relocate([ [int(r.split('-')[0]),int(r.split('-')[1])] for r in refFeature[1].split("|")],int(fPos)))
            self.data = (self.refLoc,self.geneLoc)
            
            
            
            #print self.data
            self.rawLine = self.handle.readline().split()

        except ValueError:
            self.open = False; self.rName = None
       

    def nextGappedLine(self):
        try:
            tmpLine = self.handle.readline().split()
            self.rName,self.samStrand,refFeature,self.fPos,self.blank,self.cigar,self.blank,self.blank,self.blank,self.read,self.qual,self.subs,self.blank = tmpLine
            refFeature = refFeature.split(":")
            self.gapLocs = cigarToLoc( int(self.fPos)-1, self.cigar)
            if len(refFeature) > 1:
                geneID,geneAlt,geneGroup,refChr,refStrand,refType = refFeature[0].split("|")
                genomeLocs = (refChr,[int(refFeature[1].split("-")[0])+genePos for genePos in self.gapLocs])
                geneLocs = (refFeature[0],self.gapLocs)
                self.data = (geneLocs,genomeLocs,tmpLine)
            else:
                print "WTF"
                sys.exit()
            self.rawLine= self.handle.readline().split()
        except ValueError:
            self.open = False; self.rName = None



    def nextSamLine(self):
            
        try:
            self.rName,self.samStrand,self.fName,self.fPos,self.blank,self.cigar,self.blank,self.blank,self.blank,self.read,self.qual,self.subs,self.blank = self.handle.readline().split()
            self.fPos = int(self.fPos) - 1; self.subs= int(self.subs.split(":")[-1]); self.refData = self.fName.split("|")
            self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]
            self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
        except ValueError:
            self.open = False; self.rName = None

    def nextGenomeLine(self):
        myLine = self.handle.readline().split()
        try:
            self.rName,self.samStrand,self.fName,self.fPos,self.cigar,self.read,self.qual,self.subs = myLine[0],myLine[1],myLine[2],int(myLine[3]),myLine[5],myLine[9],myLine[10],int(myLine[11].split(":")[-1])
            self.geneID,self.refStrand,self.chr,self.seqLen = None,"+",self.fName,len(self.read)-1
            self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
        except IndexError:
            self.open = False; self.rName = None



class InfoTemplate:
    def __init__(self,TYPE,data):

        if TYPE == "REFLOC":
            self.chr,self.strand,self.locs = data

        if TYPE == "GENEMAP":
            self.type,self.id,self.hugo,self.strand  = data
















 






























#!/usr/bin/env python

import os
import sys
import difflib
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs

from Sequence import *

#from ..gtTools.seq import split *


##########################################################################################################################################
#####################################################  MAPPED READ  CLASS START ############################################################
##########################################################################################################################################



class MapRead:
    def __init__(self,mapLines,strandRule):
 
        
        self.strandRule = strandRule
        if strandRule == "OPPOSITE":
            self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): False,("-","-"): False}
        elif strandRule == "MATCH":
            self.strandKey = {("+","-"): False,("-","+"): False,("+","+"): True,("-","-"): True}
        else:
            self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): True,("-","-"): True}

       
        self.validRef = dd(lambda: False); self.validRef["chr"] = True
        self.line = mapLines
        self.minOverlap = 5
        

        self.switchStrand = {'0':'+','16':'-','+':'0','-':'16'}
        
        if mapLines.format == "MAPPING":
            self.loadRead = self.pullMappingLine

        
        elif mapLines.format == "SAM":
            self.loadRead = self.pullSamLine
            #if mapLines.refType == "INTERGENIC":
            #    self.loadRead = self.pullGenomeLine
            #else:
            #    self.loadRead = self.pullSamLine
        

                
    



    def pullMappingLine(self):

        self.name = self.line.rName

        self.readLines = [self.line]

         

        self.name = self.line.rName
        self.subs = self.line.subs
        self.seq  = self.line.seq
        self.qual = self.line.qual
        
        #self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        #self.sense = False
        
        self.readLines = [self.line.data]
        self.refLocs =   [self.line.refLoc]
        self.geneLocs =  [self.line.geneLoc]
        self.line.next()
        while self.name == self.line.rName:
            if self.line.data != self.readLines[-1]:
                self.readLines.append(self.line.data)
                if self.line.refLoc != self.refLocs[-1] and self.line.geneLoc != self.geneLocs[-1]:
                    self.refLocs.append(self.line.refLoc); self.geneLocs.append(self.line.geneLoc)
            self.line.next()
        
        self.processLines2()
        #print self.name
        #print self.refLocs
        #print self.geneLocs
        #self.processLines()

      
       
    def pullSamLine(self):

        self.name = self.line.rName; # self.subs = self.line.subs; self.seq = self.line.seq; self.qual = self.line.qual
        #if self.line.gapped:
            
        #self.name, self.read, self.subs, self.qual, self.firstStrand = self.line.rName, self.line.read, self.line.subs, self.line.qual, self.line.samStrand
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        #self.ref = self.read
        self.readLines = [self.line.data]
        while self.name == self.line.rName:
            if self.line.data != self.readLines[-1]:
                self.readLines.append(self.line.data)
            self.line.next()
        if self.line.style == 'gapFeature':
            self.processGapLines()
#            readMaps.append([(self.line.chr,self.line.refStrand,self.hgPos),(self.line.geneID,self.genePos),(self.switchStrand[self.line.samStrand],self.line.refStrand)])
#            self.line.next()
#        if len(readMaps) == 1:
#            self.hgLoc, self.geneLoc, self.sense = readMaps[0][0], readMaps[0][1], self.strandKey[readMaps[0][2]] 
#            self.hgUniq, self.geneUniq = True, True 
#        else:
#            if len(readMaps)<12:
#                self.multiMaps = readMaps
#                self.disambiguate()
#            else:
#                self.invalid = True

















    def pullSamLine2(self):

        self.name, self.read, self.subs, self.qual, self.firstStrand = self.line.rName, self.line.read, self.line.subs, self.line.qual, self.line.samStrand
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        self.ref = self.read
        readMaps = []
        while self.name == self.line.rName:

            self.relocateTuple(self.key[self.line.fName], self.line.fLocs) 
            readMaps.append([(self.line.chr,self.line.refStrand,self.hgPos),(self.line.geneID,self.genePos),(self.switchStrand[self.line.samStrand],self.line.refStrand)])
            self.line.next()
        if len(readMaps) == 1:
            self.hgLoc, self.geneLoc, self.sense = readMaps[0][0], readMaps[0][1], self.strandKey[readMaps[0][2]] 
            self.hgUniq, self.geneUniq = True, True 
        else:
            if len(readMaps)<12:
                self.multiMaps = readMaps
                self.disambiguate()
            else:
                self.invalid = True
                return 
        
    
    def pullGenomeLine(self):

        self.name, self.read, self.subs, self.qual, self.firstStrand = self.line.rName, self.line.read, self.line.subs, self.line.qual, self.line.samStrand
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        self.ref = self.read
        
        readMaps = []

        while self.name == self.line.rName:
            readMaps.append([(self.line.chr,self.line.refStrand,self.line.fLocs),(self.line.geneID,self.line.fLocs),(self.switchStrand[self.line.samStrand],self.line.refStrand)])
            self.line.next()
        if len(readMaps) == 1:
            self.hgLoc, self.geneLoc, self.sense = readMaps[0][0], readMaps[0][1], self.strandKey[readMaps[0][2]] 
            self.hgUniq, self.geneUniq = True, True 
        else:
            if len(readMaps)<12:
                self.multiMaps = readMaps
                self.disambiguate()
            else:
                self.invalid = True
                return 
        
    






########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################




    def strandCheck(self):
        if self.strandRule:
            strandScrs = []
            strandScrs = [ self.strandKey[r[0][1],r[1][4]] for r in self.readLines]
            if len(set(strandScrs)) > 1:
                self.readLines = [self.readLines[i] for i in range(len(strandScrs)) if strandScrs[i]]
                return True
            else:
                return strandScrs[0]
        return True

    def filterCheck(self):
        FILTER_MAPS = [(r[1][0],r[1][1],r[1][2]) for r in self.readLines if r[1][-1]=="FILTER"]
        self.refTYPE = "FILTER"
        self.geneSTATUS = "FILTER"
        if len(FILTER_MAPS)>0:
             self.filtered = True
             self.filterGenes = list(set([x[1]+"|"+x[0] for x in FILTER_MAPS]))
             self.filterNames = list(set([x[2] for x in FILTER_MAPS]))
             #print self.filterGenes
             #print self.filterNames
             
             return True 

        
    def removeFalse(self,TYPE_FILTERS = ["EXON","KJXN","SPLICE=AG-AT","SPLICE=CT-AC"],GROUP_FILTERS=["exCANON_CANON_exCANON","exCANON_NOVEL_exCANON"]):



        ## 0) REMOVE DUPLICATES ##


        self.readLines.sort()
        dupList = [0] + [i for i in range(1,len(self.readLines)) if self.readLines[i] != self.readLines[i-1]]
        if len(dupList) < len(self.readLines):
            self.readLines = [self.readLines[i] for i in dupList]
        
        ## 1) FILTER FOR INVALID CHRS ##
     
        chrVal  = [i for i in range(len(self.readLines)) if self.validRef[self.readLines[i][0][0][0:3]]]
        if len(chrVal) < len(self.readLines) and len(chrVal) > 0:
            self.readLines = [self.readLines[i] for i in chrVal]
            if len(self.readLines) == 1:    return 

        ## 2) FILTER FOR NOVEL/KNOWN JXNS --- VALTYPES HEIRARCHY  ##
        
        types   = [r[1][5] for r in self.readLines]
        if len(set(types))>1:
            for T in TYPE_FILTERS:
                if T in types:
                    valTypes = [i for i in range(len(types)) if types[i] == T]
                    self.readLines = [ self.readLines[i] for i in valTypes]
                    break
                    
        groups = [r[1][2] for r in self.readLines]
        
        if len(set(groups)) >1:
            for G in GROUP_FILTERS:
                if G in groups:
                    valGroups = [i for i in range(len(groups)) if groups[i] == G]
                    self.readLines = [ self.readLines[i] for i in valGroups]
                    break
            
            
        if len(self.readLines) == 1:    return 
       
        ## 2b) FILTER FOR SAME PATH ## 
        if len(set([(r[0][2][0],r[0][2][-1]) for r in self.readLines]))==1 and (self.readLines[0][0][2][0],self.readLines[0][0][2][-1]) in [r[0][2] for r in self.readLines]:
            self.readLines = [[(r[0][0],r[0][1],(r[0][2][0],r[0][2][-1])),r[1],r[2]]  for r in self.readLines]


        ## 3) FILTER FOR POOR NAMING ###

        if len(set([r[0] for r in self.readLines])) == 1 and (len(set([r[1][0] for r in self.readLines])) ==1 or len(set([r[1][1] for r in self.readLines]))==1):
                self.readLines = self.readLines[0:1]
                return

        ## 4) FILTER FOR OVERLAPS 
        if len(set([r[1][1] for r in self.readLines])) == 1:
            leftCnt = len(set([r[0][2][0:len(r[0][2])-2] for r in self.readLines])); leftLen= max([r[0][2][-1]-r[0][2][-2] for r in self.readLines])
            rightCnt = len(set([r[0][2][2::] for r in self.readLines])); rightLen= max([r[0][2][1]-r[0][2][0] for r in self.readLines])
            if (leftCnt ==1 and leftLen < self.minOverlap) or (rightCnt == 1 and rightLen < self.minOverlap):
                self.readLines = self.readLines[0:1]
                return



    def trimOffsets(self):
        if not self.hgUniq: return
        if self.locs[1]-self.locs[0]  < self.minOverlap:
            OFFSET=self.locs[1]-self.locs[0]+1
            self.locs = self.locs[2::]
            self.seq=self.seq[OFFSET::]; self.ref=self.ref[OFFSET::]; self.qual = self.qual[OFFSET::]
        if self.locs[-1]-self.locs[-2] < self.minOverlap:
            OFFSET = self.locs[-1]-self.locs[-2]+1
            self.locs = self.locs[0:len(self.locs)-2]
            self.seq = self.seq[0:len(self.seq)-OFFSET]; self.ref = self.ref[0:len(self.ref)-OFFSET]; self.qual = self.qual[0:len(self.qual)-OFFSET]



    @staticmethod 
    def determine_splice_string(geneObj):
        if geneObj.refTYPE == "KJXN":
            SPLICE_TYPE = "exCANON_CANON_exCANON"
        elif geneObj.refTYPE == "NJXN":
            SPLICE_TYPE = "exCANON_NOVEL_exCANON"
        elif geneObj.refTYPE == "ITRN":
            SPLICE_TYPE = "INTRON_EXON_BOUNDARY"
        else:
            SPLICE_TYPE = geneObj.geneGroup
       
        if "exNOVEL" in SPLICE_TYPE.split("_") or "intNOVEL" in SPLICE_TYPE.split("_"):
            SPLICE_DATA = geneObj.refType.split("=")[1]+"|"+SPLICE_TYPE
        else:
  ii  kfkfk          SPLICE_DATA = "CANON"+"|"+SPLICE_TYPE

        return SPLICE_DATA

    def determineSpliceSites(self):

        if not self.hgUniq or len(self.locs)<3:
            self.spliced = False 
        else:
            self.spliced = True
            SPLICE_SITE=",".join([str(self.locs[i-1])+"-"+str(self.locs[i]) for i in range(2,len(self.locs),2)])

            if self.geneUniq:
                self.spliceSite = "|".join([self.chr,self.geneID,self.geneAltID,self.determine_splice_string(self),SPLICE_SITE])
            else:
                myIDS=",".join([m.geneID for m in self.multiGenes])
                myALTS=",".join([m.geneAltID for m in self.multiGenes])
                self.spliceSite = "|".join([self.chr,myIDS,myALTS,self.determine_splice_string(self),SPLICE_SITE])
                    
                


    def processLines2(self):

        self.strandCheck()

    def processLines(self):
        
        self.hgUniq, self.geneUniq, self.filtered, self.ambiguous, self.repetitive,self.spliced = False,False,False,False,False,False
        self.multiLocs, self.multiGenes = None, None
        
        self.mapTYPE     = "UNIQUE"
        self.featureTYPE = None 
        self.spliceTYPE  = "UNSPLICED"

        
        if self.filterCheck(): return 
        self.sense = self.strandCheck()

        if len(self.readLines) > 1: self.removeFalse()
        self.trimOffsets()


       

        HG_MAPS   = len(set([r[0] for r in self.readLines]))
        GENE_MAPS =  len(set([r[1] for r in self.readLines]))

        if GENE_MAPS == 1:
            UNIQUE_GENE_WEIGHT = 1.0 / len(self.readLines)
            AMBIG_GENE_WEIGHT  = 0.0
        else:
            UNIQUE_WEIGHT     = 0.0
            AMBIG_GENE_WEIGHT = 1.0

    
        for r in self.readLines:
            print r



        

    


        self.chr,self.strand,self.locs = self.readLines[0][0]
        self.geneID,self.geneAltID,self.geneGroup,self.refChr,self.refStrand,self.refTYPE = self.readLines[0][1] 
        self.ref = self.readLines[0][2]

        if len(set([r[0] for r in self.readLines])) == 1: self.hgUniq   = True
        if len(set([r[1] for r in self.readLines])) == 1: self.geneUniq = True
        

        self.determineSpliceSites()




        if not self.hgUniq:
            self.multiLocs =  [InfoTemplate("HGLOC",(r[0],r[1],r[2])) for r in list(set([R[0] for R in self.readLines]))]
            if len(set([r[2] for r in self.readLines])) == 1:  self.repetitive = True
            else:                                              self.ambiguous = True
            
        if not self.geneUniq:
            self.multiGenes = [ InfoTemplate("GENE",(r[0],r[1],r[2],r[3],r[4],r[5])) for r in list(set([R[1] for R in self.readLines]))]


    def processGapLines(self):
        MIN=0; self.gapData = None
        if len(self.readLines[0][0][0].split("|"))>1:
            ### CASE 1 --- GAPPED TO GENES ###
            if len(self.readLines) == 1 and self.readLines[0][0][1][3]-self.readLines[0][0][1][2] > MIN and self.readLines[0][0][1][1]-self.readLines[0][0][1][0] > MIN:
                pD = self.readLines[0]
                self.gapData = " ".join([str(s) for s in [pD[1][0],pD[1][1][1],pD[1][1][2],pD[0][0],pD[0][1][1],pD[0][1][2]]])
                

        
            


       

        
        


      


#####################################################################################################################################################
    






    def samString(self):
        myChr,myStrand,myPos = self.hgLoc
        cigar = "".join([str(myPos[i]-myPos[i-1]+1)+"M" if i%2==1 else str(myPos[i]-myPos[i-1]-1)+"N" for i in xrange(1,len(myPos))])
        return "\t".join([self.name,self.switchStrand[myStrand],myChr,str(myPos[0]),'255',cigar,'*','0','0',self.read,self.qual,'NM:i:'+str(self.subs)])



















class InfoTemplate:
    def __init__(self,TYPE,data):
        if TYPE == "GENE":
            self.geneID,self.geneAltID,self.geneGroup,self.refChr,self.refStrand,self.refType = data
        
        elif TYPE=="HGLOC":
            self.chr,self.strand,self.locs = data 





























































'''










