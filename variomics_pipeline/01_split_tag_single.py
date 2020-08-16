# removes PCR 5 mer on p7 side prime end and appends molecular barcode from p5 (read1) side to read name
import argparse, sys, gzip
from itertools import izip

def parse_args():
    prs = argparse.ArgumentParser()

    prs.add_argument("fastq_R1", type=str, help="fastq file read 1")
    prs.add_argument("out1", type=str, help="fastq out file read 1")
    prs.add_argument("--PCRnmer", type=int, default=5, help="PCR n-mer length")
    prs.add_argument("--molTagLen", type=int, default=10, help="RT molecular tag length")
    prs.add_argument("--ltrim", type=int, default=20, help="bp after molTag to trim")
    prs.add_argument("--rtrim", type=int, default=20, help="bp after PCRnmer to trim")
    prs.add_argument("--maxLen", type=int, default=110, help="max fragment length after trimming")
    x = prs.parse_args()
    return x

if __name__ == "__main__":
    args = parse_args()

    with gzip.open(args.fastq_R1, "rb") as f1, gzip.open(args.out1, "wb") as w1:

        count = 0
        for i in f1:
            if count%4 == 0: 
                name1 = i
            if count%4 == 1:
                seq1 = i
            if count%4 == 3:
                qual1 = i
                

                # name shuffling stuff, remove mol tag and nmer sequence
                # append mol tag to sequence names
                new_name1 = "@{0}:{1}".format(seq1[:args.molTagLen], name1.lstrip("@"))
                new_seq1  = seq1[args.molTagLen+args.ltrim:-(args.PCRnmer+args.rtrim)]
                new_qual1 = qual1[args.molTagLen+args.ltrim:-(args.PCRnmer+args.rtrim)]

                if len(new_seq1) > args.maxLen:
                    count += 1
                    continue
                else:
                    w1.write("{0}{1}\n+\n{2}\n".format(new_name1, new_seq1, new_qual1))

            if count % 4000 == 0:
                print "Read: {0}".format(count/4)
            count +=1
