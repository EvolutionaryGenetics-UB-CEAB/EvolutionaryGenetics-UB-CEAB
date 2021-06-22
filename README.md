This repository contains the script classifyBlastOut.py

This script aims to classify the blast matches obtained after a blast search as exonic, intronic and intergenic. To do so, the script requires a blast seach output (outfmt 6) and a gff file with the annotations (gzipped or not)

Dependencies and requirements: 

- Python3
- gffutils: https://pythonhosted.org/gffutils/

Usage:

python classifyBlastOut.py -b BLAST.out -g GFF.gff.gz -o OUT
