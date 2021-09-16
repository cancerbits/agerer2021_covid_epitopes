# Supplementary code repository for: SARS-CoV-2 escapes CD8+ T cell surveillance via mutations in MHC-I restricted epitopes
Benedikt Agerer, Maximilian Koblischke, Venugopal Gudipati, Luis Fernando Montaño-Gutierrez, Mark Smyth, Alexandra Popa, Jakob-Wendelin Genger, Lukas Endler, David M. Florian, Vanessa Mühlgrabner, Marianne Graninger, Stephan W. Aberle, Anna-Maria Husa, Lisa Ellen Shaw, Alexander Lercher, Pia Gattinger, Ricard Torralba-Gombau, Doris Trapin, Thomas Penz, Daniele Barreca, Ingrid Fae, Sabine Wenda, Marianna Traungott, Gernot Walder, Winfried F. Pick, Volker Thiel, Franz Allerberger, Elisabeth Puchhammer-Stöckl, Wolfgang Weninger, Gottfried Fischer, Wolfgang Hoepler, Erich Pawelka, Alexander Zoufaly, Rudolf Valenta, Christoph Bock, Wolfgang Paster, René Geyeregger, Matthias Farlik, Florian Halbritter, Johannes B. Huppa, Judith H. Aberle and Andreas Bergthaler

## Abstract:

CD8+ T cell immunity to SARS-CoV-2 has been implicated in COVID-19 severity and virus control. Here, we identified nonsynonymous mutations in MHC-I–restricted CD8+ T cell epitopes after deep sequencing of 747 SARS-CoV-2 virus isolates. Mutant peptides exhibited diminished or abrogated MHC-I binding in a cell-free in vitro assay. Reduced MHC-I binding of mutant peptides was associated with decreased proliferation, IFN-γ production, and cytotoxic activity of CD8+ T cells isolated from HLA-matched patients with COVID-19. Single-cell RNA sequencing of ex vivo expanded, tetramer-sorted CD8+ T cells from patients with COVID-19 further revealed qualitative differences in the transcriptional response to mutant peptides. Our findings highlight the capacity of SARS-CoV-2 to subvert CD8+ T cell surveillance through point mutations in MHC-I–restricted viral epitopes.


# Data preparation and loading

1. Download the GEO directory into the folder where all the Github scripts are located. cd to that directory.
2. Rename the GEO directory as GEOupload if it is not already called as such. 
3. In the command line, type 'bash prepdirs.sh'. this will create the following directories: 

sars042- with files ready for analysis 

sars060- with files ready for analysis

plots- to dump all figures

(Note for non linux/unix environment users: the prepdirs.sh script just separates the files from each patient into individual directories, and renames each file so that they are called exactly 'features.tsv.gz',, 'barcodes.tsv.gz' and 'matrix.mtx.gz'. This is a circumstancial requirement for the R package Seurat, and can easily do this manually without the script.)


4. The R code is not automated. Start R, make sure you have all the required packages described, and run the code chunk by chunk. 

5. For the VDJ analysis, input files need to be obtained from EGA manually with approval of the responsible Data Access Committee.


## Links:

* Gene expression omnibus (GEO): <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166651">https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166651</a>
* European Genome-Phenome Archive (EGA): <a href="https://ega-archive.org/studies/EGAS00001005060">https://ega-archive.org/studies/EGAS00001005060</a>
* Paper: <a href="https://doi.org/10.1126/sciimmunol.abg6461">https://doi.org/10.1126/sciimmunol.abg6461</a>
