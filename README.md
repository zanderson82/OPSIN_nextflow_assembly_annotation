# Introduction
This repository contains a workflow that resolves the gene copy-number, order, and phasing for the opsin genes located at chromosome Xq28. The assembly and annotation steps of this workflow were used in (Anderson et al., 2026) Long-read sequencing with targeted assembly of the opsin locus accurately evaluates genes in expressed positions. https://www.medrxiv.org/content/10.64898/2026.03.17.26348636v1

## Usage
This workflow is written for nextflow version 26
```
nextflow run main.nf \
--bam_dir \  path to directory of bam files
--input_suffix \  input bam file suffix
--region_name \  name or identifier for samples in batch
--region \  coordinates in chr-start-end format 
--metadata_table \  list of samples
--output_dir \  publish location for outputs
--final_output_name \  final summary file name
--nested_bams \  flag that looks for bams in nested sub-directories ${bam_dir}/*/${sample_id}*${input_suffix} (default is FALSE; to run just add the flag to the nextflow command)
-resume
```

### Reference genome specific coordinates:
If your starting input bam file has been aligned to the GRChg38 reference genome, then you will use chrX:153121316-155216212 as your coordinates.

If your starting input bam file has been aligned to the T2T-CHM13 reference genome, then you will use chrX:151389254-153479422 as your coordinates.

### Note about input bam directory
- If the bams are nested, use the --nested_bams flag
## Dependencies and environments
All dependencies should be available through conda, docker or github
- samtools - version 1.22 or newer
- hifiasm
- minimap2 - version 2.28 or newer
- exonerate

### YAML files for samtools, exonerate, and minimap2 are located in ``./resources/``
This workflow is designed to use conda environments for the different steps or modules (found in ``./modules/``).
The conda environment that is to be used for each step is denoted in ``./nextflow.config``. Once you create your own environments with the YAML files, you will change the paths to each respective environment.

Example:
This is what the ``./nextflow.config`` currently looks like, starting at line 22.

```
process {
    withLabel: 'extract_reads' {
        conda = '/usr/share/millerlab/samtools-1.22'
        cpus = 10
    }

    withLabel: 'run_hifiasm_XY' {
        cpus = 10
    }

    withLabel: 'run_hifiasm_XX' {
        cpus = 10
    }

    withLabel: 'align_to_assembly' {
        conda = '/usr/share/millerlab/minimap-2.28'
        cpus = 10
    }

    withLabel: 'run_exonerate' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }

    withLabel: 'convert_gff_to_bed' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }

    withLabel: 'analyze_haplotype' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }
    withLabel: 'concatenate_results' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }
}
```

The paths to each conda environment must be changed to your machine-specific paths for the workflow to run.


## Input and output file formats
The starting file for this workflow must be an **aligned bam file**. The bam index is also needed. Note that the reference genome will effect the genomic coordinates that you use. 

The other file you will need is a tab-separated metadata file that has the sample identifier (Sample ID) and the Sex (XX or XY)

There are two summary output files:
1. A summary annotation file that each sample and haplotype is appended to. This file has the following columns: 

|Column name|Contents|
|-----------|--------|
|sample_id|Sample identifier from the original metadata file|
|sex|Sex of the sample (can be XX or XY)|
|haplotype|hap1 or hap2 for XX and primary for XY samples|
|structure|Order of genes annotated on a haplotype (e.g., L-M)|
|lw_count|Number of OPN1LW genes annotated on the haplotype|
|mw_count|Number of OPN1MW genes annotated on the haplotype|
|total_genes|Total number of genes annotated on the haplotype|
|lcr_count|Number of locus control region(s) (LCR) annotated on the haplotype|
|total_contigs|Number of contigs assigned to said haplotype with annotations|
|contigs_with_lcr|Number of contigs with LCR annotations|
|contigs_without_lcr|Number of contigs without an LCR annotation|
|orphan_genes|Genes annotated on contigs that don't have an LCR annotation|
|arrays_found|Number of LCR + L or M genes found|
|is_reverse|Indicator if the array was assembled in reverse|
|orientation_ambiguous|Indicates that there were an equal number of annotations on + and - strands|
|primary_contig|Name of contig that is marked as primary (LCR annotation + most gene annotations if there are multiple contigs with arrays)|
|primary_lcr_position|Coordinate of primary LCR annotation on its contig (most helpful if there are multiple LCR annotations)|
|primary_lcr_ratio|Ratio of mapq0 that map to the LCR annotation site to the total number of reads|
|primary_lcr_reads|Number of reads that map to the primary LCR annotation site|
|primary_lcr_mapq0|Number of reads that map to the primary LCR annotation site with a mapq score of 0|


