
#Place your fastq files in C:\Users\00090473\AppData\Local\R\win-library\4.2\ShortRead\extdata\E-MTAB-1147
#https://compgenomr.github.io/book/filtering-and-trimming-reads.html

library(QuasR)
library(ShortRead)


# obtain a list of fastq file paths
fastqFiles <- system.file(package="ShortRead",
                          "extdata/E-MTAB-1147",
                          c("Lep14_1_L48BG_TGCCGGTCAG-CCTGATACAA_L001_R1.fastq.gz",
                            "Lep14_1_L48BG_TGCCGGTCAG-CCTGATACAA_L001_R2.fastq.gz")
)

# defined processed fastq file names
outfiles <- paste(tempfile(pattern=c("processed_1_",
                              "processed_2_")),".fastq",sep="")

# process fastq files
# remove reads that have more than 1 N, (nBases)
# trim 3 bases from the end of the reads (truncateEndBases)
# Remove ACCCGGGA patern if it occurs at the start (Lpattern)
# remove reads shorter than 40 base-pairs (minLength)
preprocessReads(fastqFiles, outfiles, 
                nBases=1,
                truncateEndBases=3,
                Lpattern="ACCCGGGA",
                minLength=40)





# obtain a list of fastq file paths
fastqFile <- system.file(package="ShortRead",
                          "extdata/E-MTAB-1147",
                          "Lep14_1_L48BG_TGCCGGTCAG-CCTGATACAA_L001_R1.fastq.gz")

# read fastq file
fq = readFastq(fastqFile)

# get quality scores per base as a matrix
qPerBase = as(quality(fq), "matrix")

# get number of bases per read that have quality score below 20
# we use this
qcount = rowSums( qPerBase <= 20) 

# Number of reads where all Phred scores >= 20
fq[qcount == 0]



# write out fastq file with only reads where all 
# quality scores per base are above 20
writeFastq(fq[qcount == 0], 
           paste(fastqFile, "Qfiltered", sep="_"))
