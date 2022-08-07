# Mitochondria Pipeline
Snakefile to call Mitochondria Short Variant Discovery. Converting the Terra [WDL scripts](https://app.terra.bio/#workspaces/help-gatk/Mitochondria-SNPs-Indels-hg38/workflows/help-gatk/1-MitochondriaPipeline) into sankefiles. Follows logic described on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-). Input fasta files and blacklist files for ChrM were downloaded from [here](gsutil.sh)

## How to run
/mnt/storage/apps/anaconda3/bin/snakemake -s /home/mi724/Tools/MitochondriaPipeline/MitochondriaPipeline.snakefile --cluster-config /home/mi724/Tools/MitochondriaPipeline/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete

## Steps in order
1. SubsetBamToChrM
2. RevertSam
3. AlignAndMarkDuplicates - On Not Shifted MT Reference Genome
4. AlignAndMarkDuplicates - On Shifted MT Reference Genome
5. CollectWgsMetrics
6. M2 as CallMt
7. M2 as CallShiftedMt


## Steps in each WDL file

### SubsetBamToChrM.wdl
1. SubsetBamToChrM
2. RevertSam
3. AlignAndCall.AlignAndCall as AlignAndCall
4. CoverageAtEveryBase 
5. SplitMultiAllelicSites
### AlignAndCall.wdl
1. AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt
2. AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt
3. CollectWgsMetrics
4. M2 as CallMt
5. M2 as CallShiftedMt
6. LiftoverAndCombineVcfs
7. MergeStats
8. Filter as InitialFilter
9. SplitMultiAllelicsAndRemoveNonPassSites
10. GetContamination
11. Filter as FilterContamination
### AlignmentPipeline.wdl
1. GetBwaVersion
2. AlignAndMarkDuplicates
