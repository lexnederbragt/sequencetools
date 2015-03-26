# Ad hoc script to fix problem of legacy illumina files
# read 1 and read 2 are in the same file, 
# read 1 first, but in random order and with orphans
# the script pairs up the reads into two separate outputfiles
# restoring the order
# orphans go into a separate file
# run as python split_in_pairs.py file.fastq

import sys, os
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def write_pair(readID, read1, read2, fh_out_R1, fh_out_R2):
    """Write paired reads to two files"""
    fh_out_R1.write("@%s\n%s\n+\n%s\n" % (readID + "/1", read1[0], read1[1]))
    fh_out_R2.write("@%s\n%s\n+\n%s\n" % (readID + "/2", read2[0], read2[1]))
    pass

def write_orph(readID, read, fh_out_orphan):
    """Write orphan read to separate file"""
    fh_out_orphan.write("@%s\n%s\n+\n%s\n" % (readID, read[0], read[1]))
    pass


fastqIn = sys.argv[1]

# create output filenames
# Input: sometext.fq.gz or sometext.fq
# Out read 1: sometext.R1.fq
# Out read 1: sometext.R2.fq
# Out orphan reads: sometext.orhan.fq
fileName = os.path.basename(fastqIn.rstrip(".gz"))
fileOutname, ext = os.path.splitext(fileName)
fastqOutR1   = "%s.R1%s" %(fileOutname, ext)
fastqOutR2   = "%s.R2%s" %(fileOutname, ext)
fastqOutOrph = "%s.orphan%s" %(fileOutname, ext)

# determine whether gzipped
# open file appropriately
if fastqIn.endswith(".gz"):
    fh_in = gzip.open(fastqIn, 'r')
else:
    fh_in = open(fastqIn, 'r')


# dictionary to keep encountered reads
reads1 = {}
reads2 = {}

# open output files
fh_out_R1     = open(fastqOutR1, 'wb')
fh_out_R2     = open(fastqOutR2, 'wb')

# parse input fastq file
for title, seq, qual in FastqGeneralIterator(fh_in):
    pairID = title[-2:]
    readID = title[:-2]
    if pairID == "/1":
        # read 1
        # seen? then write as paired
        if readID in reads2:
            read1 = [seq, qual]
            read2 = reads2[readID]
            write_pair(readID, read1, read2, fh_out_R1, fh_out_R2)
            # remove from dictionary
            del reads2[readID]
        else:         
           # add to dictionary
           reads1[readID] = [seq, qual]
    elif pairID == "/2":
        # read 2
        # seen? then write as paired
        if readID in reads1:
            read1 = reads1[readID]
            read2 = [seq, qual]
            write_pair(readID, read1, read2, fh_out_R1, fh_out_R2)
            # remove from dictionary
            del reads1[readID]
        else:
           # add to dictionary
            reads2[readID] = [seq, qual]
    else:
        raise Exception("Read identifier does not end in /1 or /2: %s" %title)
fh_in.close()
fh_out_R1.close()
fh_out_R2.close()

# print out orphans, i.e. the remaining reads in the dictionaries
fh_out_orphan = open(fastqOutOrph, 'wb')
for readID in reads1:
    write_orph(readID + "/1", reads1[readID], fh_out_orphan)
for readID in reads2:
    write_orph(readID + "/2", reads2[readID], fh_out_orphan)
fh_out_orphan.close()
