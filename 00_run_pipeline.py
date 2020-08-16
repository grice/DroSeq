import sys, shlex, subprocess
import batchSubmit as batch


__threads__ = 2
__btDB__    = "bt2_index/select"
__qual__    = 30

skip = False

# read the manifest file and generate a list of jobs
manifest = []
for line in open(sys.argv[1]):
    # remove empty lines and commented out lines
    if len(line.rstrip()) == 0: continue
    if line.lstrip()[0] == "#": continue

    manifest.append( line.rstrip().split() )
#for i in manifest:
#    print len(i), i

#sys.exit()


# use sickle to trim the reads, only want high quality reads
# https://github.com/najoshi/sickle
jobs = []
manifest_next = []
for i,j,k,l,m,n in manifest:
    r1_out = "trim."+i 
    r2_out = "trim."+j
    manifest_next.append([r1_out, r2_out])
     
    # paired end mode, dont trim 5prime end output gzipped file
    line = "sickle pe -f {0} -r {1} -o {2} -p {3} -x -g -t sanger -q {5} -s {4} -l 100".format(i, j, r1_out, r2_out, "single."+i, __qual__) 
    jobs.append(shlex.split(line))

if not skip:
    batch.batchSubmit(jobs, __threads__)


# use flash to merge the reads
jobs = []
manifest_next = []
for i,j,k,l,m,n in manifest:
    r1_out = "merge."+i 
    r2_out = "merge."+j
    out = "merge." + k
    manifest_next.append([out + ".extendedFrags.fastq.gz", k,l])
    expOverlap = l 
    # paired end mode, dont trim 5prime end output gzipped file
    line = "flash {0} {1} -z -o {2} -t {3} -M {4} -m {5}".format(i,j,out,__threads__, int(expOverlap)+6, int(expOverlap)-6)

    print line
    if not skip:
        subprocess.call(shlex.split(line))

    # remove the unmapped reads
    line = "rm merge.{0}.notCombined_1.fastq.gz merge.{0}.notCombined_2.fastq.gz".format(k)
    if not skip:
        subprocess.call(shlex.split(line))




# append the molecular barcode and do some read trimming to remove primer binding sites
jobs = []
barcoded = []
for i,j,l in manifest_next:
    #r1_out = "bc." + i
    #r2_out = "bc." + j
    name     = j

    r1_out = "bc." + name + ".fastq.gz"
    barcoded.append([r1_out, name])

    # run the molecular tagging script
    line  = "python 01_split_tag_single.py {0} {1} --maxLen {2}".format(i,r1_out,int(l)-50)
    jobs.append(shlex.split(line))

print line
if not skip:
    batch.batchSubmit(jobs, __threads__)

# align the files in the manifest using bowtie2 in threaded mode
# output a bamfile
bamout = []
for i,j in zip(manifest, barcoded):
    r1,name = j[0], j[1]
    sam = "{0}.bam".format(name)
    bamout.append(sam)
    samfile = open(sam, "wb")

    # align using bowtie2
    line = "bowtie2 --norc -p {2} -x {0} -U {1}".format(__btDB__, r1, __threads__)

    # pull out sequences that are the exact length
    # !skip indels for now
    line2= 'grep -E "{0}M|@"'.format(int(i[3])-55)

    # pipe to a bam file
    line3= "samtools view -bS -"

    # run the commands with safe pipes
    p1 = subprocess.Popen(shlex.split(line), stdout=subprocess.PIPE)
    #p2 = subprocess.Popen(shlex.split(line2), stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(shlex.split(line3), stdin=p1.stdout, stdout=samfile) 
    p1.stdout.close()
    #p2.stdout.close()
    p1.wait()
    p3.wait()

txtout = []
for i,j in zip(manifest, bamout):
    name = i[2]
    sequence = i[4]
    refname = i[5]
    outname = name + "_mutfind.txt"
    txtout.append(outname)

    line = "python 02_map_to_sequence.py {0} {1} {3} {2}".format(j,sequence, outname, refname)
    subprocess.call(shlex.split(line))
    #jobs.append(shlex.split(line))

if not skip:
    batch.batchSubmit(jobs, __threads__)

sys.exit()
# data analysis portion
for i,j in zip(manifest, bamout):
    name=i[2]
    bamfile = j
    analysisName = name + "_analysis.txt"
    line = "python 02_count_and_deplex.py {0} {1}".format(bamfile, analysisName)
    subprocess.call(shlex.split(line))
