# Mitochondria Pipeline
Snakefile to call Mitochondria Short Variant Discovery. Converting the Terra [WDL scripts](https://app.terra.bio/#workspaces/help-gatk/Mitochondria-SNPs-Indels-hg38/workflows/help-gatk/1-MitochondriaPipeline) into sankefiles. Follows logic described on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-).

## Steps

### MitochondriaPipeline.snakefile
1. SubsetBamToChrM
2. RevertSam
3. AlignAndCall.AlignAndCall as AlignAndCall
4. CoverageAtEveryBase 
5. SplitMultiAllelicSites
### AlignmentPipeline.snakefile
1. GetBwaVersion
2. AlignAndMarkDuplicates
### AlignAndCall.snakefile
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



