install.packages(ATACseqQC)
## ----echo=FALSE, results="hide", warning=FALSE, message=FALSE-----------------
suppressPackageStartupMessages({
  library(ATACseqQC)
  library(ChIPpeakAnno)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(phastCons100way.UCSC.hg19)
  library(MotifDb)
  library(GenomicAlignments)
})
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ----eval=FALSE---------------------------------------------------------------
 library(BiocManager)
 BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
            "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene",
            "phastCons100way.UCSC.hg19"))
library(ATACseqQC)
library(ChIPpeakAnno)
library(MotifDb)
library(GenomicAlignments)
library(motifStack)
library(Rsamtools)
library(GenomicRanges)
library(txdbmaker)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
## -----------------------------------------------------------------------------

 
 ## load the library
library(ATACseqQC)
## input the bamFile from the ATACseqQC package 
bamfile <- system.file("extdata", "~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", 
                       package="ATACseqQC", mustWork=TRUE)
bamfile.labels <- gsub(".bam", "", basename(bamfile))

bamfile <- BamFile("path/to/your.bam", index = "path/to/your.bam.bai")
bamfile <- BamFile("~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", index = "~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam.bai")
seqinfo(bamfile)
quickBamFlagSummary(bamfile)
scanBamHeader(bamfile)
gal <- readGAlignments(bamfile)
yieldSize(bamfile) <- 100000
open(bamfile)
repeat {
  chunk <- scanBam(bamfile)
  if (length(chunk[[1]]$qname)==0) break
  # process chunk ...
}
close(bamfile)
yieldSize(bamfile) <- NA

readBAM <- function(bamfile){
  
  bam <- scanBam(bamfile)
  .unlist <- function (x){
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  bam_field <- names(bam[[1]])
  
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  
  #return a list that can be called as a data frame
  return(bam_df)
}

# Load the bam file
bam1 <- readBAM(bamfile)


bamfile <- system.file("~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", "bam", package="ATACseqQC", mustWork=TRUE)
bamfile.labels <- gsub(".shifted.bam", "", basename("~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam"))
bamfile.labels <- gsub("shifted.bam", "", basename("~/Desktop/Run149May2025/TR46/Af293BamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam"))
## ----eval=FALSE---------------------------------------------------------------
#source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))

## -----------------------------------------------------------------------------
bamQC(bam1, outPath=NULL)
res <- bamQC(bamfile = "~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", index = "~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam.bai")
res <- bamQC(bamfile = "~/Desktop/Run149May2025/TR46/Af293BamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", index = "~/Desktop/Run149May2025/TR46/Af293BamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam.bai")
estimateLibComplexity(readsDupFreq("~/Desktop/Run149May2025/TR46/Af293BamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam"))

## -----------------------------------------------------------------------------
## generate fragement size distribution
fragSize <- fragSizeDist("~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", bamfile.labels)
fragSize <- fragSizeDist("~/Desktop/Run149May2025/TR46/Af293BamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", bamfile.labels)

## -----------------------------------------------------------------------------
## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                    "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                  "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                  "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                  "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                  "U2"))
library(Rsamtools)
bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags
## files will be output into outPath
outPath <- "splited"
dir.create(outPath)
## shift the coordinates of 5'ends of alignments in the bam file
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## if you don't have an available TxDb, please refer
## GenomicFeatures::makeTxDbFromGFF to create one from gff3 or gtf file.
seqlev <- "chr1" ## subsample data for quick run
chromosomeinfo <- GCF_000002655.1_ASM265v1_genomic
seqinformation <- seqinfo(chromosomeinfo)
which <- as(seqinformation[seqlev], "GRanges")
which <- as(seqinformation, "GRanges")
gal <- readBamFile("~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path("~/Desktop/Run149May2025/TR46/ShiftedBamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam", "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)

## -----------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
gal1 <- readGAlignments("~/Desktop/Run149May2025/TR46/Af293BamFiles/149-N15_ATAC_TR46_cyp51A_hph__Rep2_S111_L004_R1_001_val_1.fq.gz.shifted.bam")
txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
gffFile <- system.file("extdata", "GFF3_files", "a.gff3", package="txdbmaker")
txdb <- makeTxDbFromGFF("~/Desktop/AfuGenomesTRUE/GCF_000002655.1_ASM265v1_genomic.gff", 
                        format=c("gff3"),
                        dataSource = NA,
                        organism = NA,
                        taxonomyId = NA,
                        chrominfo = NULL,
                        miRBaseBuild = NA,
                        metadata = NULL,
                        dbxrefTag = "transcript_id")

readGFF("~/Desktop/GCF_000002655.1_ASM265v1_genomic.gff")
seqinfo(txs)
seqinfo(gal1)
txs <- transcripts(txdb)
pt <- PTscore(gal1, txs)
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")

## -----------------------------------------------------------------------------
nfr <- NFRscore(gal1, txs)
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))

## -----------------------------------------------------------------------------
tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")

## -----------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
## run program for chromosome 1 only
txs <- txs[seqnames(txs) %in% "chr1"]
genome <- Hsapiens
## split the reads into NucleosomeFree, mononucleosome, 
## dinucleosome and trinucleosome.
## and save the binned alignments into bam files.
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath,
                              conservation=phastCons100way.UCSC.hg19)
## list the files generated by splitGAlignmentsByCut.
dir(outPath)

## ----eval=FALSE---------------------------------------------------------------
#  objs <- splitBam(bamfile, tags=tags, outPath=outPath,
#                   txs=txs, genome=genome,
#                   conservation=phastCons100way.UCSC.hg19)

## ----eval=FALSE---------------------------------------------------------------
#  ## split reads by fragment length
#  ## NOT RUN IN THIS example
#  objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)

## ----fig.height=4, fig.width=4------------------------------------------------
library(ChIPpeakAnno)
bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))
## Plot the cumulative percentage of tag allocation in nucleosome-free 
## and mononucleosome bam files.
cumulativePercentage(bamfiles[1:2], as(seqinformation["chr1"], "GRanges"))

## ----fig.height=8, fig.width=4------------------------------------------------
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
## estimate the library size for normalization
(librarySize <- estLibSize(bamfiles))
## calculate the signals around TSSs.
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
#plot heatmap
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)

## ----fig.show="hide"----------------------------------------------------------
## get signals normalized for nucleosome-free and nucleosome-bound regions.
out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")

## -----------------------------------------------------------------------------
## rescale the nucleosome-free and nucleosome signals to 0~1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")

## -----------------------------------------------------------------------------
## foot prints
library(MotifDb)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
print(CTCF[[1]], digits=2)
sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                         genome=genome, ## Don't have a genome? ask ?factorFootprints for help
                         min.score="90%", seqlev=seqlev,
                         upstream=100, downstream=100)

## ----fig.height=6, fig.width=6------------------------------------------------
featureAlignedHeatmap(sigs$signal, 
                      feature.gr=reCenterPeaks(sigs$bindingSites,
                                               width=200+width(sigs$bindingSites[1])), 
                      annoMcols="score",
                      sortBy="score",
                      n.tile=ncol(sigs$signal[[1]]))

sigs$spearman.correlation
sigs$Profile.segmentation

## -----------------------------------------------------------------------------
vp <- vPlot(shiftedBamfile, pfm=CTCF[[1]], 
            genome=genome, min.score="90%", seqlev=seqlev,
            upstream=200, downstream=200, 
            ylim=c(30, 250), bandwidth=c(2, 1))

distanceDyad(vp, pch=20, cex=.5)

## -----------------------------------------------------------------------------
path <- system.file("extdata", package="ATACseqQC", mustWork=TRUE)
bamfiles <- dir(path, "*.bam$", full.name=TRUE)
gals <- lapply(bamfiles, function(bamfile){
  readBamFile(bamFile=bamfile, tag=character(0), 
              which=GRanges("chr1", IRanges(1, 1e6)), 
              asMates=FALSE)
})
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicAlignments)
plotCorrelation(GAlignmentsList(gals), txs, seqlev="chr1")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
