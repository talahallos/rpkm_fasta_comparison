## Description: Comparing RPKM genomic file to fasta gene file
## Author: Talah
## Initiation date: April 26

install.packages("AnaCoDa")
download.packages(AnaCoDa, directory)
library("AnaCoDa")
#install.packages("stringr")
library(stringr)
library(dplyr)
library(tibble)
#Need a different way to read in the fasta file, using "head(readLines(filepath), n = 1000)" shows lines starting with quotation marks!
#install.packages("seqinr")
library(seqinr)
install.packages("AnnotationDbi")

# File Names:
# GSE75897_RiboZero_RPKMs.txt
# GSE75897_RiboZero_RPKMs.gz
# orf_genomic_1000.fasta

setwd("/Users/talahallos/Desktop/research")
# list.files() can be used to orient within a file - see what's there



# AnaCoDa Framework

genome <- initializeGenomeObject(file = "orf_genomic_1000.fasta")
# changed to gene.assignment from geneAssignment
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))
model <- initializeModelObject(parameter = parameter, model = "ROC")
#adjusted sphi, num.mixtures, and gene.assignment which are now random values
parameter <- initializeParameterObject(genome = genome, sphi = c(0.5, 2), num.mixtures = 2, gene.assignment = sample.int(2, length(genome), replace = T))
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)), init.sepsilon = sepsilon)
model <- initializeModelObject(parameter = parameter, model = "ROC",  fix.observation.noise = TRUE)



gz_path <- "/Users/talahallos/Desktop/research/GSE75897_RiboZero_RPKMs.gz"
 # Open a connection to the output .gz file using gzfile()
gz_conn <- gzfile(gz_path, "wb")

gz_conn <- gzfile(gz_path, "rb")
file_contents <- readLines(gz_conn)
close(gz_conn)

#.txt file to .gz file!!!!
text_content <- readLines("GSE75897_RiboZero_RPKMs.txt")
gzfile_conn <- gzfile("GSE75897_RiboZero_RPKMs.gz", "wb")
writeLines(text_content,gzfile_conn)

# READING .GZ FILES
gzfile_conn <- gzfile("GSE75897_RiboZero_RPKMs.gz", "rb")
text_content <- readLines(gzfile_conn)
close(gzfile_conn)
cat(text_content, sep = "\n")

# Defining terms to the files & data
rpkm_genome <- read.delim("GSE75897_RiboZero_RPKMs.txt", header = FALSE, sep = "\t")
rpkm_names <- rpkm_genome[,1]
fasta_genome <- head(readLines("/Users/talahallos/Desktop/research/orf_genomic_1000.fasta"), n = 1000)
f_genome <- initializeGenomeObject(file = "orf_genomic_1000.fasta")
fasta_names <- getNames(f_genome)

common_names <- intersect(rpkm_names, fasta_names)
length(common_names)
#*****^^^ (rm())

# Trying to sort out a smaller set of names, 100, to work with.
common_names_subset <- head(common_names, n=100)
# getNames(f_genome) being compared to txt.genome

# Defined Variables:
# fasta_names
# rpkm_names
# common_names

#now wanting to make a table where the matches and their position are highlighted.

# Empty matches table created:
fasta_match_pos_tbl <- tibble(name = NA, pos = NA)

#only the first 10 lines to test
for(name in fasta_names[1:10]) {
  pos <- str_which(rpkm_names, name)
  n_pos <- length(pos)
  if(n_pos !=1) pos <- NA
  fasta_match_pos_tbl <- add_row(fasta_match_pos_tbl, name = name, pos = pos)
}

#now functional, whole thing
#now giving error?
for(name in fasta_names) {
  pos <- str_which(rpkm_names, name)
  n_pos <- length(pos)
  if(n_pos !=1) pos <- NA
  fasta_match_pos_tbl <- add_row(fasta_match_pos_tbl, name = name, pos = pos)
}

# to display more lines, use print(fasta_match_pos_tbl, n = #)




print(matches)

for (name in fasta_names) {
  position <- str_which(common_names, "^Y*")
  if (length(position) > 0) {
    matches <- c(position, line)
  }
}


# continuing with AnaCoDa 
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
runMCMC(mcmc = mcmc, genome = genome, model = model)

# realized was running for "orf_genomic_1000.fasta" file instead of the matches file!!
# trying to now make the file a list of matches to run the AnaCoDa framework on the matches
# What I've tried and the errors recieved as a result:
genome <- initializeGenomeObject(common_names)
# Error: Expecting a single string value: [type=character; extent=5090]
genome <- initializeGenomeObject(fasta_match_pos_tbl)
# Error: Expecting a single string value: [type=list; extent=2].

#now trying to make fasta_match_pos_tbl into string to correct the error message
fasta_match_pos_tbl_as_strings <- as.data.frame(lapply(fasta_match_pos_tbl, as.character))





