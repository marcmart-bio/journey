#Bioconductor basic elements 

### IRanges ###
library(IRanges)

###
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir2 <- IRanges(start = c(1,3,5), width = 3)

start(ir1) #retrieve vector
width(ir2) <- 1 # resize IRanges object

names(ir1) <- paste("A", 1:3, sep = "") # name IRanges object

c(ir1, ir2) #concatenate objects
resize(ir1, width = 1, fix = "center")

###
ir1 <- IRanges(start = c(1,3,5), width = 1)
ir2 <- IRanges(start = c(4,5,6), width = 1)

union(ir1, ir2) #equivalent to reduce(c(ir1,ir2))

intersect(ir1, ir2)

###
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3)

ir1
ir2 

ov <- findOverlaps(ir1, ir2)
ov

queryHits(ov)
subjectHits(ov)

unique(queryHits(ov)) #exact elements of the query that overlaps anything in the subject

args(findOverlaps)

countOverlaps(ir1, ir2) # simple count the number of overlaps of each element of the query

nearest(ir1, ir2) #which range in ir2 are close to ir1


### GenomicRanges - GRanges ###

#BiocManager::install("GenomicRanges")
library(GenomicRanges)
gr <- GRanges(seqnames = c("chr1"), strand = c("+", "-", "+"), ranges = IRanges(start = c(1,3,5), width = 3))
gr

flank(gr, 5)
promoters(gr)

seqinfo(gr) 

seqlengths(gr) <- c("chr1" <- 10)
seqinfo(gr)
seqlevels(gr)

gaps(gr) #stuff in chrs not covered by a range in GRanges

seqnames(gr) <- c("chr1", "chr2") #need to use seqlevels to make ths attribution

seqlevels(gr) <- c("chr1", "chr2") 
seqnames(gr) <- c("chr1", "chr2", "chr1")

seqlevels(gr) <- c("chr2", "chr1") 
sort(gr)


# set a genome may protect further computational biology - e.g.
genome(gr) <- "hg19" 
gr
seqinfo(gr)
gr2 <- gr 
genome(gr2) <- "hg18"
findOverlaps(gr, gr2)

###
ir = IRanges(start = 1:3, width = 2)
ir
df = DataFrame(ir = ir, score = rnorm(3)) # DataFrame() is an special data frame object in GRanges 
df
df[1,1]
df[1,2]
df$ir

gr <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = c("+", "-", "+"), ranges = IRanges(start = c(1,3,5), width = 3))
gr

values(gr) <- DataFrame(score = rnorm(3)) #including metadata with DataFrame in GRanges object
gr

gr2 <- GRanges(seqnames = c("chr1"), strand = "*", ranges = IRanges(start = c(1,3,5), width = 3))
gr2
gr
findOverlaps(gr, gr2)
findOverlaps(gr, gr2, ignore.strand = TRUE)
subsetByOverlaps(gr, gr2)

findOverlaps(gr2, gr)
subsetByOverlaps(gr2, gr)

###
df = data.frame(chr = "chr1", start = 1:3, end = 4:6, score = rnorm(3)) #classic R dataframe
df

makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

###

gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr

seqlevels(gr, force = TRUE) <- "chr1" #does not work in the current build

#alternatively
dropSeqlevels(gr, "chr2")
keepSeqlevels(gr, "chr1")

gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))
newStyle = mapSeqlevels(seqlevels(gr), "NCBI")
newStyle

gr = renameSeqlevels(gr, newStyle)
gr

### AnnotationHub ###
BiocManager::install("AnnotationHub")
library(AnnotationHub)

ah <- AnnotationHub()

ah[1] #look for snapshot
#ah[[1]] #retrieve the sequence online

#BiocManager::install("rtracklayer")

unique(ah$dataprovider)

ah2 <- subset(ah, species == "Homo sapiens")
ah2

query(ah, c("H3K4me3", "Gm12878"))

ah2 <- display(ah)
ah2

### querying Drosophila

ah[ah$species == "Drosophila willistoni"] #access
wil <- ah[["AH108769"]] #retrieve
wil

dros <- query(ah, c("Drosophila")) # 765 records

wil <- display(dros)
wil 

wilgr <- wil[["AH108769"]]

