#!/bin/bash

#run tasks up until the GetContamination step (parallel process)
/mnt/storage/apps/anaconda3/bin/snakemake -s MitochondriaPipeline.snakefile --cluster-config /home/mi724/Tools/MitochondriaPipeline/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete

#run the GetContamination step alone without parallelism
/mnt/storage/apps/anaconda3/bin/snakemake -s GetContamination.snakefile --cluster-config /home/mi724/Tools/MitochondriaPipeline/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 1 --rerun-incomplete

#run everything after GetContamination step (parallel process)
/mnt/storage/apps/anaconda3/bin/snakemake -s FilterContamination.snakefile --cluster-config /home/mi724/Tools/MitochondriaPipeline/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete
