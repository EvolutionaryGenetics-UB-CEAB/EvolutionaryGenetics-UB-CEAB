This repository contains the pipeline to classify blast matched into four genic categories (genic, intergenic, intronic, exonic) and count them. The pipeline contains two python scripts:

- classifyBlastOut.py: this script aims to classify the blast matches obtained after a blast search as exonic, intronic and intergenic. To do so, the script requires a blast seach output (outfmt 6) and a gff file with the annotations (gzipped or not)

Dependencies and requirements: 
- Python3
- gffutils: https://pythonhosted.org/gffutils/

Usage: 
python classifyBlastOut.py -b BLAST.out -g GFF.gff.gz -o OUT

- countingClasifiedHits_general.py: script to parse the output of classifyBlastOut.py and count the number of unique blast hits and the number of hits in each genic categories (genic, intergenic, intronic, exonic).

IMPORTANT: those blast hits located in both exonic and intronic regions are classified as exonic, independently if they belong to different transcripts of the same gene or to different genes.

Usage: 
python countingClasifiedHits_general.py -i IN -o OUT



