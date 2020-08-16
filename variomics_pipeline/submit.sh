#$ -N varmir_windels_21
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -l m_mem_free=2G
#$ -l h_rt=20000
#$ -j y

#python 00_run_pipeline.py manifest2.txt
python 00_run_pipeline.py manifest2_21.txt
