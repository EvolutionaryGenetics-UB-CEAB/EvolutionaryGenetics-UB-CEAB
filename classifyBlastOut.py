# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 11:36:32 2020

@author: cpegueroles

"""

import argparse,sys
import gffutils

parser = argparse.ArgumentParser(description="script to classify the blast matches obtained after a blast search as exonic, intronic and intergenic. To do so, the script requires a blast seach output (outfmt 6) and a gff file with the annotations (gzipped or not) ")

parser.add_argument("-b", "--blast", dest="blast", required=True, help="blast output in outfmt 6")
parser.add_argument("-g", "--gff", dest="gff", required=True, help="genome annotations in gff format")
parser.add_argument("-o", "--out", dest="out", required=True, help="out file")

args = parser.parse_args()
blast = str(args.blast)
gff = str(args.gff)
out = open(str(args.out),'w')

########################### create a db from gtf
db = gffutils.create_db(gff, ':memory:', merge_strategy="create_unique", disable_infer_genes=True)
db = db.update(db.create_introns(merge_attributes=False)) 

########################### FUNCTIONS
def read_blast(blast): #read blast out and store it in a dictionary
    blast_dict={}
    i=1
    for row in open(blast, 'r').readlines():
        row = row.rstrip().split("\t")
        if len(row) !=12:
            sys.exit("\tWARNING! Please check your blast out file! It should be outfmt 6")
        ID=['HIT_%s' % i]
        row=row+ID
        key= row[0]
        val= []
        if not key in blast_dict:
            blast_dict[key]=val
        blast_dict[key].append(row)
        i+=1
    return blast_dict

def locate_hit(blast_dict,myfeature,locate_dict):
    #blast_dict=read_blast(blast)
    found=0; intron=0
    for feature in db.features_of_type(myfeature): #exon, mRNA, CDS, gene, region
       #print (feature)
       print(feature[0],feature[3], feature[4], feature.id, myfeature) #scafold, start, end
       for k, v in blast_dict.items():
           for x in v:
               #x[1]: scafold; x[8]: start blast hit; x[9]: end blast hit
               #feature[0]: scafold; feature[3]: start feature; feature[4]: end feature
               #check whether the cordinates of the blast results are ordered (smallest coordinate in first column)
               if x[8]<x[9]:
                   blast_start=x[8]; blast_end=x[9]
               else:
                   blast_start=x[9]; blast_end=x[8]
               if x[1]==feature[0]: #check the scafold
                   if int(blast_start)<=feature[3] and int(blast_end)<=feature[4] and int(blast_end)>=feature[3] and int(blast_start)<=feature[4]: #caseII
                       #print(x,feature.id)
                       found=1
                   if int(blast_start)<=feature[3] and int(blast_end)>=feature[4]: #caseI
                       #print(x,feature.id)
                       found=1
                   if int(blast_start)>=feature[3] and int(blast_end)<=feature[4]: #caseIII
                       #print(x,feature.id)
                       found=1; intron=1
                   if int(blast_start)>=feature[3] and int(blast_end)>=feature[4] and int(blast_end)>=feature[3] and int(blast_start)<=feature[4]: #caseIV
                       #print(x,feature.id)
                       found=1
               if found==1 and myfeature !='intron':
                   #print(x,feature)
                   key=k
                   x.append(feature.id)
                   if not key in locate_dict:
                       locate_dict[key]=[]
                   locate_dict[key].append(x)
               if intron==1 and myfeature =='intron':
                   key=k
                   x.append(feature.id)
                   if not key in locate_dict:
                       locate_dict[key]=[]
                   locate_dict[key].append(x)
               found=0; intron=0
                     
    return locate_dict

def intergenic(blast_dict,locate_dict):
    myblastHitdict={}      
    for k, v in blast_dict.items(): 
        for x in v:
            key=x[12]
            val= []
            if not key in myblastHitdict:
                myblastHitdict[key]=val
            myblastHitdict[key].append(x)
            
    mylocateHitdict={}      
    for k, v in locate_dict.items(): 
        for x in v:
            key=x[12]
            val= []
            if not key in mylocateHitdict:
                mylocateHitdict[key]=val
            mylocateHitdict[key].append(x)
            
    
    for k, v in myblastHitdict.items():
        if not k in mylocateHitdict.keys():
            v[0].append('intergenic') 
            mylocateHitdict[key].append(v[0]) 
                
    return mylocateHitdict

########################## MAIN
locate_dict={}
blast_dict=read_blast(blast)
mydict=locate_hit(blast_dict,'gene', locate_dict)
mydict=locate_hit(blast_dict,'pseudogene', locate_dict)
mydict=locate_hit(blast_dict,'exon', locate_dict)
mydict=locate_hit(blast_dict,'intron', locate_dict)

mydict=intergenic(blast_dict,locate_dict)


for k, v in mydict.items():
    v=[list(x) for x in set(tuple(x) for x in v)]   #to eliminate repeated elements 
    for x in v:
        #print('\t%s' % x)
        out.write('%s\n'% '\t'.join(x))
