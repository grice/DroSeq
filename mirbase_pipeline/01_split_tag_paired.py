# removes PCR 5 mer on p7 side prime end and appends molecular barcode from p5 (read1) side to read name
import argparse, sys, gzip
from itertools import izip

def parse_args():
    prs = argparse.ArgumentParser()

    prs.add_argument("fastq_R1", type=str, help="fastq file read 1")
    prs.add_argument("fastq_R2", type=str, help="fastq file read 2")
    prs.add_argument("out1", type=str, help="fastq out file read 1")
    prs.add_argument("out2", type=str, help="fastq out file read 2")
    prs.add_argument("--PCRnmer", type=int, default=5, help="PCR n-mer length")
    prs.add_argument("--molTagLen", type=int, default=10, help="RT molecular tag length")
    prs.add_argument("--R1trim", type=int, default=20, help="bp after R1 to trim")
    prs.add_argument("--R2trim", type=int, default=20, help="bp after R2 to trim")
    x = prs.parse_args()
    return x

if __name__ == "__main__":
    args = parse_args()

    with gzip.open(args.fastq_R1, "rb") as f1, gzip.open(args.fastq_R2, "rb") as f2, \
            gzip.open(args.out1, "wb") as w1, gzip.open(args.out2, "wb") as w2:

        count = 0
        for i,j in izip(f1,f2):
            if count%4 == 0: 
                name1, name2 = i, j
            if count%4 == 1:
                seq1,  seq2  = i, j
            if count%4 == 3:
                qual1, qual2 = i, j
                

                # name shuffling stuff, remove mol tag and nmer sequence
                # append mol tag to sequence names
                new_name1 = "@{0}:{1}".format(seq1[:args.molTagLen], name1.lstrip("@"))
                new_seq1  = seq1[args.molTagLen+args.R1trim:]
                new_qual1 = qual1[args.molTagLen+args.R1trim:]

                new_name2 = "@{0}:{1}".format(seq1[:args.molTagLen], name2.lstrip("@"))
                new_seq2  = seq2[args.PCRnmer+args.R2trim:]
                new_qual2 = qual2[args.PCRnmer+args.R2trim:]

                w1.write("{0}{1}+\n{2}".format(new_name1, new_seq1, new_qual1))
                w2.write("{0}{1}+\n{2}".format(new_name2, new_seq2, new_qual2))

            if count % 4000 == 0:
                print "Read: {0}".format(count/4000)
            count +=1
