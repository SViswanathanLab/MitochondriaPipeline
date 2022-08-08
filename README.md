# Mitochondria Pipeline
Snakefile to call Mitochondria Short Variant Discovery. Converting the Terra [WDL scripts](https://app.terra.bio/#workspaces/help-gatk/Mitochondria-SNPs-Indels-hg38/workflows/help-gatk/1-MitochondriaPipeline) into sankefiles. Follows logic described on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-). Input fasta files and blacklist files for ChrM were downloaded from [here](gsutil.sh)

GetContamination rule needs the haplocheckCLI.jar file from the below repo.
```
git clone https://github.com/leklab/haplocheckCLI.git
```

## How to run
Snakemake uses the following format to run snakefiles on a high performance cluster.
```
/path/to/snakemake -s /path/to/snakefile --cluster-config  /path/to/cluster_qsub.yaml --cluster "qsub #variables inside the cluster_qsub.yaml file to add to qsub commands" --jobs #of_jobs_to_submit_at_a_time 
```

The below is how to run this snakefile (MitochondriaPipeline.snakefile) from inside my personal directory.
```
/mnt/storage/apps/anaconda3/bin/snakemake -s /home/mi724/Tools/MitochondriaPipeline/MitochondriaPipeline.snakefile --cluster-config /home/mi724/Tools/MitochondriaPipeline/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete
```

## Files
#### 1. config/config.yaml
The config file contains paths to scripts, tools and input files, aw well as parameters.
#### 2. config/samples.yaml
This file contains the path to the input bam files from TCGA. Follow the format provided in the file when adding more files.
#### 3. config/cluster_qsub.yaml
This file contains parameters for the qsub submissions. Default applies to all rules, but if a specific rule name is specified in the file, it will take on those parameters.
#### 4. code/coverage.R
This file contains the R commands that were used in the WDL script for the CollectWgsMetrics task.
#### 5. code/CoverageAtEveryBase.R
This file contains the R commands that were used in the WDL script for the CoverageAtEveryBase task.
#### 6. gsutil.sh
This file contains the copying commands from google storage to argos for the ChrMT input files.
#### 7. MitochondriaPipeline.snakefile
Snakefile to run the analysis.

## Steps in order inside Snakefile
1. SubsetBamToChrM
2. RevertSam
3. AlignAndMarkDuplicates - On Not Shifted MT Reference Genome
4. AlignAndMarkDuplicates - On Shifted MT Reference Genome
5. CollectWgsMetrics
6. M2 as CallMt
7. M2 as CallShiftedMt
8. LiftoverAndCombineVcfs
9. MergeStats
10. Filter as InitialFilter
11. SplitMultiAllelicsAndRemoveNonPassSites
12. GetContamination
13. Filter as FilterContamination
14. CoverageAtEveryBase 
15. SplitMultiAllelicSites

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
