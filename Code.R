library(QuasR)
library(ShortRead)


# obtain a list of fastq file paths from the working directory
fastqFiles <- list.files(pattern="\\.fastq\\.gz")

# defined processed fastq file names
outfiles <- paste(tempfile(pattern=c("processed_1_", "processed_2_")), ".fastq", sep="")

# process fastq files
# remove reads that have more than 1 N, (nBases)
# trim 3 bases from the end of the reads (truncateEndBases)
# Remove ACCCGGGA pattern if it occurs at the start (Lpattern)
# remove reads shorter than 40 base-pairs (minLength)
preprocessReads(fastqFiles, outfiles, 
                nBases=1,
                truncateEndBases=3,
                Lpattern="ACCCGGGA",
                minLength=40)

# obtain a list of fastq file paths
fastqFile <- list.files(pattern="\\.fastq\\.gz")[1]

# read fastq file
fq = readFastq(fastqFile)

# get quality scores per base as a matrix
qPerBase = as(quality(fq), "matrix")

# get number of bases per read that have quality score below 20
qcount = rowSums(qPerBase <= 20) 

# Number of reads where all Phred scores >= 20
fq[qcount == 0]

# write out fastq file with only reads where all 
# quality scores per base are above 20
writeFastq(fq[qcount == 0], 
           paste(fastqFile, "Qfiltered", sep="_"))

# use this Python code to extract the .gz file:  https://github.com/Aria-Dolatabadian/Python-tips/blob/main/Extract%20fastq.gz%20to%20fastq
