nextflow.enable.dsl=2

// Pull modules
include { extract_reads } from "./modules/01_extract_reads.nf"
include { run_hifiasm_XY } from "./modules/02.1_run_hifiasm_XY.nf"
include { run_hifiasm_XX } from "./modules/02.2_run_hifiasm_XX.nf"
include { cat_assemblies } from "./modules/02.3_cat_assemblies.nf"
include { align_to_assembly } from "./modules/03_align_to_assembly.nf"
include { run_exonerate } from "./modules/04_run_exonerate.nf"
include { convert_vulgar_to_gff } from "./modules/05_convert_vulgar_to_gff.nf"
include { convert_gff_to_bed } from "./modules/06_convert_gff_to_bed.nf"
include { add_alignment_stats } from "./modules/07_add_alignment_stats.nf"
include { analyze_haplotype } from "./modules/08_analyze_haplotype.nf"
include { concatenate_results } from "./modules/09_concatenate_haplotype_analysis.nf"


// Calculate hg_size from region param
def calc_hg_size(region) {
    def matcher = region =~ /[^:]+:(\d+)-(\d+)/
    if (!matcher) return ''
    def size = matcher[0][2].toLong() - matcher[0][1].toLong()
    if (size >= 1_000_000_000) return "${Math.round(size / 1_000_000_000)}g"
    if (size >= 1_000_000)     return "${Math.round(size / 1_000_000)}m"
    return "${Math.max(1, Math.round(size / 1000))}k"
}






// Define entry workflow and channels
workflow {
    def hg_size = calc_hg_size(params.region)
// This initial step is reading in the rows of the metadata file, splitting
    samples_ch = Channel
    .fromPath(params.metadata_table)
    .splitCsv(sep: '\t', header: ['sample_id', 'sex'], skip: 2)
    .map { row -> tuple(row.sample_id.toString().trim(), row.sex.toString().trim()) }


// Channel for finding the input aligned bam file
// Resolve globs explicitly to avoid empty/ambiguous file() behavior
    input_bam_ch = samples_ch
        .map { sample_id, sex ->
            def subdir = params.nested_bams ? "*/" : ""
            def bam_pattern = "${params.bam_dir}/${subdir}${sample_id}*${params.input_suffix}"
            def bai_pattern = "${bam_pattern}.bai"
            def bam_matches = files(bam_pattern)
            def bai_matches = files(bai_pattern)

            if (!bam_matches || bam_matches.size() == 0) {
                error "No BAM matched for ${sample_id} (${sex}) with pattern: ${bam_pattern}"
            }
            if (!bai_matches || bai_matches.size() == 0) {
                error "No BAI matched for ${sample_id} (${sex}) with pattern: ${bai_pattern}"
            }
            if (bam_matches.size() > 1) {
                error "Multiple BAM matches for ${sample_id} (${sex}): ${bam_matches}"
            }
            if (bai_matches.size() > 1) {
                error "Multiple BAI matches for ${sample_id} (${sex}): ${bai_matches}"
            }

            tuple(sample_id, sex, bam_matches[0], bai_matches[0])
        }


// Extract reads from the input aligned bam file
    extract_reads(input_bam_ch)

    // Branch based on sex
    extract_reads.out.fastq
        .map { sample_id, sex, fastq -> tuple(sample_id, sex, fastq, hg_size) }
        .filter { sample_id, sex, fastq, size -> sex == "XY" }
        .set { xy_fastq_ch }

    extract_reads.out.fastq
        .map { sample_id, sex, fastq -> tuple(sample_id, sex, fastq, hg_size) }
        .filter { sample_id, sex, fastq, size -> sex == "XX" }
        .set { xx_fastq_ch }
    

    run_hifiasm_XY(xy_fastq_ch)
    run_hifiasm_XX(xx_fastq_ch)


    // Combine the hap1 and hap2 fasta file into a channel to then be concatenated
    run_hifiasm_XX.out.hap1_fa
        .join(run_hifiasm_XX.out.hap2_fa, by: [0, 1])
        .set { xx_assemblies_ch}
    
    // concatenate hap1 and hap2 from run_hifiasm_XX
    cat_assemblies(xx_assemblies_ch)

    // add the fastq file to the tuple from extract reads for XX samples
    cat_assemblies.out.combined_diploid_fa
        .join(extract_reads.out.fastq, by: [0, 1])
        .set { xx_assembly_and_fastq_for_alignment_ch }
    
    // add the fastq file to the tuple from extract reads for XY samples
    run_hifiasm_XY.out.primary_fa
        .join(extract_reads.out.fastq, by: [0, 1])
        .set { xy_assembly_and_fastq_for_alignmet_ch }

    // combine XX and XY samples into one channel for alignment
    // new tuple will be sample_id, sex, assembly_fa, fastq
    xx_assembly_and_fastq_for_alignment_ch
        .mix(xy_assembly_and_fastq_for_alignmet_ch)
        .set { joined_assemblies_and_fastqs_for_alignment_ch }


    // outputs are (sample_id, sex, bam_path) as reads_to_assembly_bam
    // and (sample_id, sex, bam_index) as reads_to_assembly_bai
    align_to_assembly(joined_assemblies_and_fastqs_for_alignment_ch)

    // set up channels for exonerate
    // add XY primary hap label
    run_hifiasm_XY.out.primary_fa
        .map { sample_id, sex, primary_fa -> tuple(sample_id, sex, primary_fa, "primary") }
        .set { xy_assembly_for_exonerate_ch }
    
    run_hifiasm_XX.out.hap1_fa
        .map { sample_id, sex, hap1_fa -> tuple(sample_id, sex, hap1_fa, "hap1") }
        .set { xx_hap1_for_exonerate_ch }

    run_hifiasm_XX.out.hap2_fa
        .map { sample_id, sex, hap2_fa -> tuple(sample_id, sex, hap2_fa, "hap2") }
        .set { xx_hap2_for_exonerate_ch }

    // tuple structure is now sample_id, sex, assembly_fa, hap_name
    xy_assembly_for_exonerate_ch
        .mix(xx_hap1_for_exonerate_ch)
        .mix(xx_hap2_for_exonerate_ch)
        .set { all_assemblies_for_exonerate_ch }

    // Run exonerate on tuple of sample_id, sex, assembly_fa and hap_name
    // Output a tupe of sample_id, sex, hap_name, and exonerate_output_vulgar
    run_exonerate(all_assemblies_for_exonerate_ch)

    // run convert_vulgar_to_gff and output a tuple of sample_id, sex, hap_name and the vulgar_to_gff_output
    convert_vulgar_to_gff(run_exonerate.out.exonerate_output_vulgar)

    // run convert gff to bed and output a tuple of sample_id, sex, hap_name, and the gff_to_bed_output
    convert_gff_to_bed(convert_vulgar_to_gff.out.vulgar_to_gff_output)

    // Create new channel with sample_id, sex, hap_name, gff_to_bed_output, reads_to_assembly_bam, reads_to_assembly_bai
    convert_gff_to_bed.out.gff_to_bed_output
        .combine(align_to_assembly.out.reads_to_assembly_bam, by: [0,1])
        .combine(align_to_assembly.out.reads_to_assembly_bai, by: [0,1])
        .set { beds_and_bams_for_adding_stats_ch }
    
    // add alignment stats and output a tuple of sample_id, sex, hap_name, and stats_gff
    add_alignment_stats(beds_and_bams_for_adding_stats_ch)

    analyze_haplotype(add_alignment_stats.out.stats_gff)


    

   analyze_haplotype.out.haplotype_analysis
        .map { sample_id, sex, hap_name, tsv -> tsv }
        .set { new_analysis_tsvs_ch }

    existing_analysis_tsvs = Channel.fromPath("${params.output_dir}/*.analysis.tsv")

    new_analysis_tsvs_ch
        .mix(existing_analysis_tsvs)
        .unique { it.name }
        .collect()
        .set { all_tsvs_ch }

    concatenate_results(all_tsvs_ch)

    }