## Description: 
## Author: Talah
## Initiation date: April 26

#install.packages(‘AnaCoDa’)
#download.packages(AnaCoDa, directory)
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
# GSE75897_RiboZero_RPKMs.txt
# orf_genomic_1000.fasta




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


con <- gzfile("/Users/talahallos/Desktop/research/GSE75897_RiboZero_RPKMs.gz", "r")



#.txt file to .gz file!!!!
text_content <- readLines("GSE75897_RiboZero_RPKMs.txt")
gzfile_conn <- gzfile("GSE75897_RiboZero_RPKMs.gz", "wb")
writeLines(text_content,gzfile_conn)

# READING .GZ FILES
setwd("path/to/folder")
gzfile_conn <- gzfile("example.txt.gz", "rb")
text_content <- readLines(gzfile_conn)
close(gzfile_conn)
cat(text_content, sep = "\n")

rpkm_genome <- read.delim("GSE75897_RiboZero_RPKMs.gz", header = FALSE, sep = "\t")
rpkm_names <-rpkm.genome[,1]


head(readLines(filepath), n = 500)
# assigned filepath <- orf_genomic_1000.fasta file!
# assigned filecontents <- head(readLines(filepath), n = 500)
cat(filecontents, sep = "\n")

# SHORTCUTS TO THE TWO FILES:
# txt.genome = read.delim("GSE75897_RiboZero_RPKMs.gz", header = TRUE, sep = "\t")
# fasta.genome = head(readLines(filepath), n = 1000)

txtgenome <- readLines("GSE75897_RiboZero_RPKMs.gz")
indices <- grep("^Y", txtgenome)
selected_lines <- txtgenome[indices]
print(selected_lines)

fastagenome <- readLines("orf_genomic_1000.fasta")
indices2 <- grep("^>", fastagenome)
length(indices2)
print(indices2)
selected_lines2 <- fastagenome[indices2]
length(selected_lines2)
print(selected_lines2)


file_fasta <-
sequences <- read.fasta("orf_genomic_1000.fasta")
my_sequence <- sequences[["orf_genomic_fasta"]]
print(my_sequence)
#didn't work

install.packages("AnnotationDbi")
#output gives link https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
av <- available.packages(filters=list())
av[av[, "Package"] == "AnnotationDbi", ]

?getNames 

sub("^>([^ ]+) .*", "\\1", selected_lines2) 
fasta_names <- sub("^>([^ ]+) .*", "\\1", selected_lines2)

genome <- initializeGenomeObject(file = "orf_genomic_1000.fasta")
fasta_names <- getNames(genome)

common_names <- intersect(rpkm_names, fasta_names)
length(common_names)
#*****^^^ (rm())

common_names_subset <- head(common_names, n=100)


getNames(genome)
#being compared to
txt.genome

# xxxxxx  common_lines <- intersect(txt.genome, fasta.genome)

fasta_names
rpkm_names
matches

for (line in fasta_names) {
  position <- str_which(common_names, )
}
# Error in str_which(common_names, ) : could not find function "str_which"


for (line in fasta_names) {
  position <- str_which(common_names, )
}
# Error in str_detect(string, pattern, negate = negate) : argument "pattern" is missing, with no default
# changing pattern to line (next to common_names)


#figure out things to this point--
#connection to github. commit it to github then start to work on and change things 

match.tbl <- tibble(name = fasta_names, position = NA)

for (name in fasta_names) {
  position <- str_which(common_names, name)
  print(position)
  if (length(position) > 0) {
    matches <- c(position, name)
  }
  print(matches)
}

print(matches)

for (name in fasta_names) {
  position <- str_which(common_names, "^Y*")
  if (length(position) > 0) {
    matches <- c(position, line)
  }
}




#trying to make the getNames(genome) list into its own file
mylist <- getNames(genome)
write.csv(as.data.frame(mylist), "/Users/talahallos/Desktop/research/mylist", row.names = FALSE)
my_list_from_file <- read.csv("/Users/talahallos/Desktop/research/mylist", header = TRUE, stringsAsFactors = FALSE)
my_list_from_file









