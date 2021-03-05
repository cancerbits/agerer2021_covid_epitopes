# agerer2021_covid_epitopes
Code repository accompanying the paper by Agerer et al. (2021), Sci. Immunol.

# Data preparation and loading

1. Download the GEO directory into the folder where all the GIT scripts are located. cd to that directory.
2. Rename the GEO directory as GEOupload if it is not already called as such. 
3. In the command line, type 'bash prepdirs.sh'. this will create the following directories: 

sars042- with files ready for analysis 

sars060- with files ready for analysis

plots- to dump all figures

(Note for non linux/unix environment users: the prepdirs.sh script just separates the files from each patient into individual directories, and renames each file so that they are called exactly 'features.tsv.gz',, 'barcodes.tsv.gz' and 'matrix.mtx.gz'. This is a circumstancial requirement for the R package Seurat, and can easily do this manually without the script.)


4. The R code is not automated. Start R, make sure you have all the required packages described, and run the code chunk by chunk. 

