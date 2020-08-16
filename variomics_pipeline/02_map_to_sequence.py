from __future__ import print_function
import sys, pysam, re
import argparse as prs
#from countMutations import parseCigarString as pCS

def parseArgs():
    arg = prs.ArgumentParser()
    arg.add_argument('bamFile', type=str, help='input bamfile that contains aligned sequences')
    arg.add_argument('sequence', type=str, help='miR sequence to compare')
    arg.add_argument("refname", type=str, help='reference sequence name')
    arg.add_argument('outFile', type=str, help='out txt file')
    arg.add_argument('--maxMutations', type=int, default=6, help="max number of mutations in a read")
    o = arg.parse_args()
    return o

def mutstring2seq(mutstring, reference):
    # catch empty mutstrings
    if mutstring == None: return reference
    
    # use the number pattern to split the sequence
    pattern = "[0-9]+"

    # arrays of letters and numbers
    posmatch = re.findall(pattern, mutstring)
    letmatch = re.split(pattern, mutstring)[1:]

    #for pos, let in zip(posmatch, letmatch):
    #    mutDict[pos] = let
    #return letmatch, posmatch[:-1]
    sequence = ""
    # iterate through the reference sequence
    for n, base in enumerate(reference):
        # check if the base needs to be changed
        if str(n) in posmatch:
            # the base is in the list of changes, iterate through changes
            for k, pos in enumerate(posmatch):
                # if it is to be changed see what to change it to
                if pos == str(n):
                    mut = letmatch[k]

                    # insertions
                    if mut[0] == "^":
                        # sequence after carrot is what to insert
                        sequence += mut[1:]

                        # catch a list with last kind of change being an insertion
                        if k + 1 >= len(posmatch):
                            sequence += base
                        
                        # check if there is another kind of mutation at the same position
                        else:
                            # if there is don't add the base after the insertion
                            if posmatch[k+1] == str(n):
                                continue
                            # for a simple insertion add the base after
                            else:
                                sequence += base
                    
                    # substitutions
                    if mut in "AGCT":
                        sequence += mut
                    # deletions
                    if mut == "-":
                        continue
        else:
            sequence += base
    return sequence

def parseCigar2(read, ref_seq):
    #ref_seq      = read.get_reference_sequence()
    align_seq    = read.query_alignment_sequence
    align_tuples = read.get_aligned_pairs(with_seq=True)

    mutCount = 0

    mutstring = ""

    # catch missing bases in 5p end
    for i in range(align_tuples[0][1]):
        mutstring += str(i) + "-"
        mutCount += 1

    last_ref = 0
    for i in range(len(align_tuples)):
        query, ref, base = align_tuples[i]

        # if there is an insertion, use the location of the next base in the sequence
        if ref == None:
            if align_tuples[i-1][1] == None:
                mutstring += align_seq[query]
            else:
                mutstring += str(last_ref + 1) + "^" + str(align_seq[query])
            mutCount += 1

        # deletions
        elif query == None:
            mutstring += str(ref) + "-"
            mutCount += 1

        #substitutions
        elif base in "agct":
            mutstring += str(ref) + str(align_seq[query])
            mutCount += 1
        
        # track the last reference base for insertion mapping
        if ref != None:
            last_ref = ref

    
    # catch 3p end deletions
    for i in range(align_tuples[-1][1]+1, len(ref_seq)):
        mutstring += str(i) + "-"
        mutCount += 1

    if len(mutstring) == 0:
        mutstring = None

    return mutstring, mutCount


if __name__ == "__main__":
    arg = parseArgs()
    
    samfile = pysam.AlignmentFile(arg.bamFile, "rb")
    w       = open(arg.outFile, "w")

    nameStore = {}
    correct = 0

    for count, read in enumerate(samfile):
        moltag = read.query_name.split(":")[0]
        cigar  = read.cigarstring
        # start  = read.query_alignment_start

        seq    = read.query_alignment_sequence
        # quals  = read.qual


        # if there is no cigar string the read didn't map
        if cigar == None: 
            continue
        #print(count,"#"*30)
        # md_tag = read.get_tag("MD")

        #print(cigar, md_tag)
        #print(read.get_aligned_pairs(with_seq=True))

        #if count == 500: break

        # align, quals, shortName, mutCount, qualFlag  = pCS(cigar, seq, quals, arg.sequence, start)

        # if shortName != parseCigar2(read, arg.sequence):
        #     print("#"*30)
        #     print(count)
        #     print(arg.sequence)
        #     print(seq, cigar)
        #     print(shortName)
        #     print(parseCigar2(read, arg.sequence))
        #     print(read.get_aligned_pairs(with_seq=True))
        #print(parseCigar2(read))

        #print shortName, qualFlag, mutCount, read.mapping_quality
        shortName, mutCount = parseCigar2(read, arg.sequence)
        reverseSequence = mutstring2seq(shortName, arg.sequence)

        # catch mistagged sequences, e.g. those that were misidentified via miseq barcode
        if read.reference_name != arg.refname: continue

        if seq != reverseSequence:
            if read.reference_name != "hsa-mir-16-1": continue
            print("#"*30)
            print(read.reference_name)
            print(arg.sequence)
            print(seq)
            print(mutstring2seq(shortName, arg.sequence))
            print(shortName)
        else:
            correct += 1

        # filter very mismatched reads and those that have low quality
        if mutCount > arg.maxMutations: continue
        
        #line = "{0} {1}\n".format(moltag, shortName)
        #w.write(line)

        if shortName not in nameStore:
            nameStore[shortName] = [moltag]
        else:
            nameStore[shortName].append(moltag)
        
        #if count == 100000: break
        if count % 10000 == 0:
            print(count)

    # moltag collapse
    for name in nameStore:
        count = len(set(nameStore[name]))
        line = "{0} {1}\n".format(name, count)
        w.write(line)

    w.close()
    print(correct)

