. /u/local/Modules/default/init/modules.sh
module load singularity

base_dir=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/PEG_experiment
indir=${base_dir}/Reference_genomes/
mapfile=${base_dir}/midas_db/peg.mapfile
outdir=${base_dir}/midas_db/

singularity exec $H2_CONTAINER_LOC/MIDAS-1.3.2.sif build_midas_db.py $indir $mapfile $outdir 