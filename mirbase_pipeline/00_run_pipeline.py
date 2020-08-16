import sys, shlex, subprocess
import batchSubmit as batch


__threads__ = 8
__btDB__    = "../bowtie_genomes/combo_rc_padded"
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
for i,j,k in manifest:
    r1_out = "trim."+i 
    r2_out = "trim."+j
    manifest_next.append([r1_out, r2_out, k])
     
    # paired end mode, dont trim 5prime end output gzipped file
    line = "sickle pe -f {0} -r {1} -o {2} -p {3} -x -g -t sanger -q {5} -s {4} -l 35".format(i, j, r1_out, r2_out, "single."+i, __qual__) 
    jobs.append(shlex.split(line))

if not skip:
    batch.batchSubmit(jobs, __threads__)


# append the molecular barcode and do some read trimming to remove primer binding sites
jobs = []
barcoded = []
for i,j,k in manifest_next:
    r1_out = "bc." + i
    r2_out = "bc." + j
    name     = k

    barcoded.append([r1_out, r2_out, name])

    # run the molecular tagging script
    line  = "python 01_split_tag_paired.py {0} {1} {2} {3}".format(i, j, r1_out,r2_out)
    jobs.append(shlex.split(line))

print(line)
if not skip:
    batch.batchSubmit(jobs, __threads__)

# align the files in the manifest using bowtie2 in threaded mode
# output a bamfile
bamout = []
for i,j in zip(manifest, barcoded):
    r1,r2,name = j[0], j[1], j[2]
    sam = "{0}.sam".format(name)
    bamout.append(sam)
    samfile = open(sam, "wb")

    # align using bowtie2
    line = "bowtie2 --norc -p {2} -x {0} -1 {1} -2 {3}".format(__btDB__, r1, __threads__, r2)

    # pipe to a bam file
    line3= "samtools view -bS -"

    # run the commands with safe pipes
    #p1 = subprocess.Popen(shlex.split(line), stdout=subprocess.PIPE)
    p1 = subprocess.Popen(shlex.split(line), stdout=samfile)
    #p3 = subprocess.Popen(shlex.split(line3), stdin=p1.stdout, stdout=samfile) 
    #p1.stdout.close()

    p1.wait()
    #p3.wait()


# data analysis portion
for i,j in zip(manifest, bamout):
    name=i[2]
    bamfile = j
    analysisName = name + "_analysis.txt"
    line = "python 02_map_to_sequence.py {0} {1}".format(bamfile, analysisName)
    subprocess.call(shlex.split(line))
