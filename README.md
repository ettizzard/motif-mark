# Motif Mark

## Goals
The motivation for this project is to generate a simple image file that visualizes the presence of gene sequence motifs along an entire set of genes via an object-oriented python script. Visualization of gene motifs is valuable for genomic analysis, as it can lead to a better understanding of gene function and aid in identifying patterns within a genome. 



## Requirements
* A ".fasta" file with gene names at the beginning of its header lines and with case sensitivity corresponding to intronic and exonic regions (i.e., exonic bases are capitalized, and introns are in lowercase)
* A ".txt" file with each motif of interest on its own line - motifs with degenerate bases are supported!
* A conda environment with [pycairo](https://anaconda.org/conda-forge/pycairo) installed and activated.

## Usage
Running motif-mark-oop.py involves a simple one-line command. 
```
python motif-mark-oop.py -f {fasta file} -m {motif text file}
```
After a successful run, the script will output a .png image file with your visualized motifs!
