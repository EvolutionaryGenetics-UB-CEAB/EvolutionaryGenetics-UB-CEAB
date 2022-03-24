#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 10:25:45 2021

@author: cpegueroles
"""


import argparse, re, os

parser = argparse.ArgumentParser(description="script to count the located blast hits (as exonic, intronic and intergenic")

parser.add_argument("-i", "--input", dest="input", required=True, help="located blast hits; the results from running classifyBlastOut.py")
parser.add_argument("-o", "--out", dest="out", required=True, help="out file")

args = parser.parse_args()
in1 = str(args.input)
out = open(str(args.out),'w')

#1. discard non-uniq matchs
#2. count:
    #EXONS: exon/id-
    #INTRONS: intron/no label for pseudogens
    #GENES: gene
    #INTERGENIC: intergenic

#example:
#CLocus_11210	NW_018114504.1	89.831	59	6	0	1	59	901252	901194	2.03e-13	76.8	HIT_2410	gene-LOC109986035
#CLocus_21184	NW_018114430.1	94.915	59	3	0	1	59	2251066	2251124	2.02e-18	93.5	HIT_3178	gene-LOC110002752	id-LOC110002752

fileName = os.path.basename(in1)
sp=fileName.split(".")[0]

########################### FUNCTIONS
def uniq_hits_dict(in1): #read blast out and store it in a dictionary
    in_dict={}
    multiple=[]
    for row in open(in1, 'r').readlines():
        split_row = row.rstrip().split("\t")
        key= split_row[0]
        val= ''
        if not key in in_dict:
            in_dict[key]=val
            in_dict[key]=row
        else:       
            #in_dict.pop(key, None)
            multiple.append(key)
        #if key in in_dict:
            #in_dict[key]=row
            
    for k in multiple:
        in_dict.pop(k, None)     #to remove loci with multiple hits
        
    i=len(set(multiple))
        
    return (in_dict,i)

def count_elements(unique_dictionary):
    count_dict={'Gene':'', 'Intergenic':'', 'Exon':'', 'Intron':''}
    gene=0; intergenic=0; exon=0; intron=0; exonintron=0; pseudoIntrons=0
    for k, v in unique_dictionary.items():
        if re.search('gene', v):
            gene +=1 #correct
            if not re.search('exon', v) and not re.search('intron', v) and not re.search('id-', v):
                pseudoIntrons +=1 #correct
    for k, v in unique_dictionary.items():
        if re.search('intergenic', v):
            intergenic +=1  #correct          
    for k, v in unique_dictionary.items():
        if re.search('exon', v):
            exon +=1 #correct; contains also introns and id-; but there are empty lines!!!           
    for k, v in unique_dictionary.items():
        if re.search('id-', v) and not re.search('exon', v):
            exon +=1 #correct; NO introns NEITHER id-; but there are empty lines!!!
    for k, v in unique_dictionary.items():
        if re.search('intron', v):
            intron +=1
    for k, v in unique_dictionary.items():
        if re.search('exon', v) or re.search('id-', v):
            if re.search('intron', v) :
                exonintron +=1 #there are also labels exon-id

    #print(exon);print(intron);print(pseudoIntrons);print(exonintron);
    count_dict['Gene']=gene
    count_dict['Intergenic']=intergenic
    count_dict['Exon']=exon
    count_dict['Intron']=intron + pseudoIntrons - exonintron
    
    #print(exon+(intron + pseudoIntrons - exonintron))
    
    return count_dict

#######main

parseFile= uniq_hits_dict(in1)
print('%s\t' % (sp))
print("number of multiple blast hits: %s" % (parseFile[1]))
print("number of unique blast hits: %s" % (len(parseFile[0].keys())))
#out.write("number of unique blast hits: %s\n" % (len(uniq_hits_dict(in1).keys())))
#uniq_hits_dict is CORRECT

out.write('%s\t' % (sp))
for k, v in count_elements(parseFile[0]).items():
    print('%s\t%s' % (k,v))
    #out.write('%s\t%s\t' % (k,v))
    out.write('%s\t' % (v))
out.write('\n')

