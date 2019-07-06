# Variant-Visualizer
### A visualization tool for visualizing the variations in any patient’s DNA sequence produced by Next Generation Sequencing, the tool can be further developed to compare various patient’s DNA and identify patterns.

### Purpose:
* This tool takes input a FASTA file and a variation (VCF) file and plots the DNA sequence With the variations.
* The tool allows the user to zoom into the plot using a slider or by left click and zoom out using the slider or by right click.
* The tool provides a miniature view to keep track of the part of the chromosome the user is in or the variations under the scope.
* The user can check the variations in any state by hovering over the variations.

#### Libraries used: Tkinter, Matplotlib, Biopython, Allel, Numpy, Pandas.

### User Guidelines:
1. The Script need a FASTA file and a VCF file in the same folder for it to work.
2. The Human genome is taken as the sample and can be downloaded from the NCBI website and given as default input, 
   user need to change the name of the FASTA file before starting.
3. The FASTA file can be broken down into chromosomes of separate files.
4. The large genome file need not to be processed every time.
5. The opening interface plots all the chromosomes with the variation plotted to choose from.
6. The user needs to choose a chromosome and click on the chromosome, the prompt shows the chromosome clicked on.
7. On closing the interface, the plotting of the selected chromosome is done.
8. The user can hover over the variations to check the actual position and variations.
9. The user needs to click on the variation to expand/zoom into the plot.
10. On choosing the variation, and extreme zooming the actual positions and the ATGC sequence can be seen.
11. The user can still hover and use slider, LEFT, RIGHT to zoom in/out and traverse the DNA strip in both directions.
---

The files needed can be downloaded from NCBI website.

Downloaded the sample Human Genome : 

* [ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/]

Reference Index of Human Genome : 
* [ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/]

---
#### NB: The script won't work without the required files in the same folder.
* _Required FASTA file: genome.fa_
* _VCF file: Any Variant file which reference the input fasta file._
---

###### Opening Interface with all the Chromosomes
![Chromosome Interface](https://github.com/nirmalya-mukherjee/Variant-Visualizer/blob/master/pics/chromosome%20interface.png)

###### Interface to zoom in, zoom out and hover over the variants with a miniature view
![Interface](https://github.com/nirmalya-mukherjee/Variant-Visualizer/blob/master/pics/Hypo%201.png)

###### Actual position with hovering and a miniature view
![Interface Actual positions](https://github.com/nirmalya-mukherjee/Variant-Visualizer/blob/master/pics/Actual%20position%201.png)

###### This is the first version of the tool developed by Nirmalya Mukherjee.

