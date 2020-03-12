#!/usr/bin/env nextflow
/*
========================================================================================
                                   RNAseq
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/rnaseq-nf
 @#### Authors
 Sarai Varona <s.varona@isciii.es>
----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Preprocessing
   - 1.1: Build aligner index
     - 1.1.1: Build STAR index
	   - 1.1.2: Build HISAT2 splice sites file
	   - 1.1.3: Build HISAT2 index
	 - 1.2: Convert GFF3 to GTF
	 - 1.3: Build BED12 file
	 - 1.4: FastQC for raw sequencing reads quality control
   - 1.5: Trimmomatic
 - 2. : Alignment
   - 2.1: Align with STAR
	 - 2.2: Align with HISAT2
	 - 2.3: preseq analysis
 - 3. : Remove duplicates
	 - 3.1: Mark duplicates
	 - 3.2: dupRadar
 - 4. : Feature counts
   - 4.1: Feature counts
	 - 4.2: Merge feature counts
 - 5. : Assembly
   - 5.1: stringtie FPKM
 - 6. : Differential expression
   - 6.1: edgeR (3 or mor samples needed)
	 - 6.2: DESeq2
 - 7. : Quality control
   - 7.1: RSeQC analysis
 - 8. : Stats
   - 8.1 : MultiQC
 - 9. : Output Description HTML
 ----------------------------------------------------------------------------------------


*/
/*
 * This is a message
 */
def helpMessage() {
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
    Options:
      --genome                      Name of iGenomes reference
      --singleEnd                   Specifies that the input is single end reads

    Strandedness:
      --forward_stranded            The library is forward stranded
      --reverse_stranded            The library is reverse stranded
      --unstranded                  The default behavior

    Aligner:
      --aligner                     Choose aligner for RNA data between 'star' or 'hisat2'. By default STAR aligner.

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index                  Path to STAR index
      --hisat2_index                Path to HiSAT2 index
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF file
      --gff                         Path to GFF3 file
      --bed12                       Path to bed12 file
      --saveReference               Save the generated reference files the the Results directory.
      --saveAlignedIntermediates    Save the BAM files from the Alignment step  - not done by default

    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_adapters_parameters Trimming parameters for adapters. <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. Default 2:30:10
      --trimmomatic_window_length   Window size. Default 4
      --trimmomatic_window_value    Window average quality required. Default 20
      --trimmomatic_mininum_length  Minimum length of reads

    Presets:
      --pico                        Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forward_stranded --clip_r1 3 --three_prime_clip_r2 3
      --fcExtraAttributes           Define which extra parameters should also be included in featureCounts (default: gene_names)
      --fcGroupFeatures             Define the attribute type used to group features. (default: 'gene_name')
      --fcGroupFeaturesType         Define the type attribute used to group features based on the group attribute (default: 'gene_name')

    Other options:
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      --sampleLevel                 Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --seqCenter                   Add sequencing center in @RG line of output BAM header

    QC options:
      --skip_qc                     Skip all QC steps apart from MultiQC
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_genebody_coverage      Skip calculating genebody coverage
      --skip_preseq                 Skip Preseq
      --skip_dupradar               Skip dupRadar (and Picard MarkDups)
      --skip_edger                  Skip edgeR MDS plot and heatmap
      --skip_multiqc                Skip MultiQC

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}
//Other
params.help = false
params.name = false
//References
params.genome = false
params.genomes = false
params.readPaths = false
//Library data
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false
params.singleEnd = false
params.pico = false
//Alignment
params.aligner = 'star'
params.seqCenter = false
//Saving files
params.saveAlignedIntermediates = false
params.saveReference = false
params.saveTrimmed = false
//Feature Counts
params.fcGroupFeatures = 'gene_name'
params.fcGroupFeaturesType = 'gene_name'
params.fcExtraAttributes = 'gene_name'
//Skip steps
params.skip_qc = false
params.skip_fastqc = false
params.skip_genebody_coverage = false
params.skip_preseq = false
params.skip_dupradar = false
params.skip_edger = false
params.skip_multiqc = false
params.skip_rseqc = false
// Defaults
sampleLevel = false
hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
subsampFilesizeThreshold = 10000000000 // Don't subsample BAMs for RSeQC gene_body_coverage if less than this
readPaths = null
star_memory = false // Cluster specific param required for hebbe

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false


ch_mdsplot_header = Channel.fromPath("$baseDir/assets/mdsplot_header.txt")
ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt")
ch_biotypes_header = Channel.fromPath("$baseDir/assets/biotypes_header.txt")
Channel.fromPath("$baseDir/assets/where_are_my_files.txt")
       .into{ch_where_trim_galore; ch_where_star; ch_where_hisat2; ch_where_hisat2_sort}

// Trimming default
params.notrim = false

// Default trimming options
params.trimmomatic_adapters_file = "\$TRIMMOMATIC_PATH/adapters/NexteraPE-PE.fa"
params.trimmomatic_adapters_parameters = "2:30:10"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "50"

forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Validate inputs
if (params.aligner != 'star' && params.aligner != 'hisat2'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}
if( params.star_index && params.aligner == 'star' ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}
else if ( params.hisat2_index && params.aligner == 'hisat2' ){
    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*")
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
}
else if ( params.fasta ){
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
           .into { ch_fasta_for_star_index; ch_fasta_for_hisat_index}
}
else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_makeBED12;
              gtf_star; gtf_dupradar; gtf_featureCounts; gtf_stringtieFPKM }
} else if( params.gff ){
  gffFile = Channel.fromPath(params.gff)
                   .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
}

if( params.bed12 ){
    bed12 = Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .into {bed_rseqc; bed_genebody_coverage}
}
if( params.aligner == 'hisat2' && params.splicesites ){
    Channel
        .fromPath(params.splicesites)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: ${params.splicesites}" }
        .into { indexing_splicesites; alignment_splicesites }
}
if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
params.multiqc_config = "${baseDir}/conf/multiqc_config.yaml"

if (params.multiqc_config){
    multiqc_config = file(params.multiqc_config)
}
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimming }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_trimming }
}


// Header log info
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
if(params.genome) summary['Genome'] = params.genome
if(params.pico) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
summary['Strandedness']     = ( unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None' )
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trimmomatic adapters file'] = params.trimmomatic_adapters_file
    summary['Trimmomatic adapters parameters'] = params.trimmomatic_adapters_parameters
    summary["Trimmomatic window length"] = params.trimmomatic_window_length
    summary["Trimmomatic window value"] = params.trimmomatic_window_value
    summary["Trimmomatic minimum length"] = params.trimmomatic_mininum_length
}
if(params.aligner == 'star'){
    summary['Aligner'] = "STAR"
    if(params.star_index)          summary['STAR Index']   = params.star_index
    else if(params.fasta)          summary['Fasta Ref']    = params.fasta
} else if(params.aligner == 'hisat2') {
    summary['Aligner'] = "HISAT2"
    if(params.hisat2_index)        summary['HISAT2 Index'] = params.hisat2_index
    else if(params.fasta)          summary['Fasta Ref']    = params.fasta
    if(params.splicesites)         summary['Splice Sites'] = params.splicesites
}
if(params.gtf)                 summary['GTF Annotation']  = params.gtf
if(params.gff)                 summary['GFF3 Annotation']  = params.gff
if(params.bed12)               summary['BED Annotation']  = params.bed12
summary['Save prefs']     = "Ref Genome: "+(params.saveReference ? 'Yes' : 'No')+" / Trimmed FastQ: "+(params.saveTrimmed ? 'Yes' : 'No')+" / Alignment intermediates: "+(params.saveAlignedIntermediates ? 'Yes' : 'No')
//summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * PREPROCESSING - Build STAR index
 */
if(params.aligner == 'star' && !params.star_index && params.fasta){
    process makeSTARindex {
        label 'high_memory'
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/../REFERENCES/star_index" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
		cpus '20'
		penv 'openmp'

        input:
        file fasta from ch_fasta_for_star_index
        file gtf from gtf_makeSTARindex

        output:
        file "star" into star_index


        script:
        def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''

        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN 20 \\
            --sjdbGTFfile $gtf \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            $avail_mem
        """
    }
}
/*
 * PREPROCESSING - Build HISAT2 splice sites file
 */
if(params.aligner == 'hisat2' && !params.splicesites){
    process makeHisatSplicesites {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/../REFERENCES/hisat_index" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeHisatSplicesites

        output:
        file "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites, alignment_splicesites

        script:
        """
        hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
        """
    }
}
/*
 * PREPROCESSING - Build HISAT2 index
 */
if(params.aligner == 'hisat2' && !params.hisat2_index && params.fasta){
    process makeHISATindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/../REFERENCES/hisat_index" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        cpus '20'
		penv 'openmp'

        input:
        file fasta from ch_fasta_for_hisat_index
        file indexing_splicesites from indexing_splicesites
        file gtf from gtf_makeHISATindex

        output:
        file "${fasta.baseName}.*.ht2*" into hs2_indices

        script:
        log.info "[HISAT2 index build] Over ${params.hisatBuildMemory} GB available, so using splice sites and exons in HISAT2 index"
        extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
        ss = "--ss $indexing_splicesites"

        """
        $extract_exons
        hisat2-build -p 20 $ss $exon $fasta ${fasta.baseName}.hisat2_index
        """
    }
}
/*
 * PREPROCESSING - Convert GFF3 to GTF
 */
if(params.gff){
  process convertGFFtoGTF {
      tag "$gff"

      input:
      file gff from gffFile

      output:
      file "${gff.baseName}.gtf" into gtf_makeSTARindex, gtf_makeHisatSplicesites, gtf_makeHISATindex, gtf_makeBED12,
            gtf_star, gtf_dupradar, gtf_featureCounts, gtf_stringtieFPKM

      script:
      """
      gffread $gff -T -o ${gff.baseName}.gtf
      """
  }
}

/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/../REFERENCES/bed12" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/01-fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trimmomatic
 */
process trimming {
    label 'low_memory'
    tag "$prefix"
    publishDir "${params.outdir}/02-preprocessing", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf(".log") > 0) "logs/$filename"
                else if (params.saveTrimmed && filename.indexOf(".fastq.gz")) "trimmed/$filename"
                else null
        }

    input:
    set val(name), file(reads) from raw_reads_trimming
    file wherearemyfiles from ch_where_trim_galore.collect()

    output:
    file '*_filtered_*.fastq.gz' into trimmed_reads,trimmed_paired_reads,trimmed_paired_reads_bwa
    file '*_unpaired_*.fastq.gz' into trimmed_unpaired_reads, trimmed_unpaired_reads_picard
    file '*_fastqc.{zip,html}' into trimmomatic_fastqc_reports, trimmomatic_fastqc_reports_picard
    file '*.log' into trimmomatic_results, trimmomatic_results_picard


    script:
    prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    java -jar $TRIMMOMATIC_PATH/trimmomatic-0.33.jar PE -threads 1 -phred33 $reads $prefix"_filtered_R1.fastq" $prefix"_unpaired_R1.fastq" $prefix"_filtered_R2.fastq" $prefix"_unpaired_R2.fastq" ILLUMINACLIP:${params.trimmomatic_adapters_file}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${name}.log
    gzip *.fastq
    fastqc -q *_filtered_*.fastq.gz
    """
}


/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}
if(params.aligner == 'star'){
    hisat_stdout = Channel.from(false)
    process star {
        label 'high_memory'
        tag "$prefix"
        publishDir "${params.outdir}/03-alignment", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") == -1) "logs/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
            }
        cpus '20'
		penv 'openmp'

        input:
        file reads from trimmed_reads
        file index from star_index.collect()
        file gtf from gtf_star.collect()
        file wherearemyfiles from ch_where_star.collect()

        output:
        set file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log
        file "where_are_my_files.txt"
        file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_filtered_)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        def star_mem = task.memory ?: params.star_memory ?: false
        def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
        seqCenter = params.seqCenter ? "--outSAMattrRGline ID:$prefix 'CN:$params.seqCenter'" : ''
        """
        STAR --genomeDir $index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN 20 \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate $avail_mem \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix $prefix $seqCenter

        samtools index ${prefix}Aligned.sortedByCoord.out.bam
        """
    }
    // Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM; bam_forSubsamp; bam_skipSubsamp  }
}


/*
 * STEP 3 - align with HISAT2
 */
if(params.aligner == 'hisat2'){
    star_log = Channel.from(false)
    process hisat2Align {
        label 'high_memory'
        tag "$prefix"
        publishDir "${params.outdir}/03-alignment", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
            }
        cpus '20'
		penv 'openmp'

        input:
        file reads from trimmed_reads
        file hs2_indices from hs2_indices.collect()
        file alignment_splicesites from alignment_splicesites.collect()
        file wherearemyfiles from ch_where_hisat2.collect()

        output:
        file "${prefix}.bam" into hisat2_bam
        file "${prefix}.hisat2_summary.txt" into alignment_logs
        file "where_are_my_files.txt"

        script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2l?/
        prefix = reads[0].toString() - ~/(_R1)?(_filtered_)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
//        seqCenter = params.seqCenter ? "--rg-id ${prefix} --rg CN:${params.seqCenter.replaceAll('\\s','_')}" : ''
//                   --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
//                   --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
        def rnastrandness = ''
        if (forward_stranded && !unstranded){
            rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
        } else if (reverse_stranded && !unstranded){
            rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
        }
        if (params.singleEnd) {
            """
            hisat2 -x $index_base \\
                   -U $reads \\
                   $rnastrandness \\
                   --known-splicesite-infile $alignment_splicesites \\
                   -p 20 \\
                   --met-stderr \\
                   --new-summary \\
                   --summary-file ${prefix}.hisat2_summary.txt \\
                   | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
            """
        } else {
            """
            hisat2 -x $index_base \\
                   -1 ${reads[0]} \\
                   -2 ${reads[1]} \\
                   $rnastrandness \\
                   --known-splicesite-infile $alignment_splicesites \\
                   --no-mixed \\
                   --no-discordant \\
                   -p 20 \\
                   --met-stderr \\
                   --new-summary \\
                   --summary-file ${prefix}.hisat2_summary.txt \\
                   | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
            """
        }
    }

    process hisat2_sortOutput {
        label 'mid_memory'
        tag "${hisat2_bam.baseName}"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: { filename ->
                if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") "aligned_sorted/$filename"
                else null
            }
        cpus '20'
		pènv 'openmp'

        input:
        file hisat2_bam
        file wherearemyfiles from ch_where_hisat2_sort.collect()

        output:
        file "${hisat2_bam.baseName}.sorted.bam" into bam_count, bam_rseqc, bam_preseq, bam_markduplicates, bam_featurecounts, bam_stringtieFPKM,bam_forSubsamp, bam_skipSubsamp
        file "${hisat2_bam.baseName}.sorted.bam.bai" into bam_index_rseqc, bam_index_genebody
        file "where_are_my_files.txt"

        script:
        def suff_mem = ("${(task.memory.toBytes() - 6000000000) / 20}" > 2000000000) ? 'true' : 'false'
        def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / 20}" : ''
        """
        samtools sort \\
            $hisat2_bam \\
            -@ 20 ${avail_mem} \\
            -o ${hisat2_bam.baseName}.sorted.bam
        samtools index ${hisat2_bam.baseName}.sorted.bam
        """
    }
}

/*
 * STEP 4 - RSeQC analysis
 */
process rseqc {
    label 'mid_memory'
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/04-rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("clipping_profile.r") > 0)                "clipping_profile/data/$filename"
            else if (filename.indexOf("clipping_profile.xls") > 0)              "clipping_profile/data/$filename"
            else if (filename.indexOf("clipping_profile.R1.pdf") > 0)           "clipping_profile/plots/$filename"
            else if (filename.indexOf("clipping_profile.R2.pdf") > 0)           "clipping_profile/plots/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction_annotation_log.txt") > 0)       "junction_annotation/logs/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_rseqc
    file index from bam_index_rseqc
    file bed12 from bed_rseqc.collect()

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    prefix = bam_rseqc.baseName - '_filteredAligned.sortedByCoord.out'
    if (params.singleEnd) {
        paired="SE"
    } else {
        paired="PE"
    }
    """
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${prefix}.infer_experiment.txt
    clipping_profile.py -i $bam_rseqc -s $paired -o ${prefix}
    junction_annotation.py -i $bam_rseqc -o ${prefix} -r $bed12
    bam_stat.py -i $bam_rseqc > ${prefix}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${prefix} -r $bed12 2> ${prefix}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${prefix} -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${prefix}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${prefix}.read_duplication
    """
}

/*
 * Step 4.1 Subsample the BAM files if necessary

bam_forSubsamp
    .filter { it.size() > params.subsampFilesizeThreshold }
    .map { [it, params.subsampFilesizeThreshold / it.size() ] }
    .set{ bam_forSubsampFiltered }
bam_skipSubsamp
    .filter { it.size() <= params.subsampFilesizeThreshold }
    .set{ bam_skipSubsampFiltered }

process bam_subsample {
    tag "${bam.baseName - '.sorted'}"

    input:
    set file(bam), val(fraction) from bam_forSubsampFiltered

    output:
    file "*_subsamp.bam" into bam_subsampled

    script:
    """
    samtools view -s $fraction -b $bam | samtools sort -o ${bam.baseName}_subsamp.bam
    """
}
 */

/*
 * Step 4.2 Rseqc genebody_coverage
 */
process genebody_coverage {
    label 'mid_memory'
    tag "${bam.baseName - '.sorted'}"
       publishDir "${params.outdir}/04-rseqc" , mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    when:
    !params.skip_qc && !params.skip_genebody_coverage

    input:
    file bam from bam_forSubsamp
    file bed12 from bed_genebody_coverage.collect()

    output:
    file "*.{txt,pdf,r}" into genebody_coverage_results

    script:
    prefix = bam.baseName - '_filteredAligned.sortedByCoord.out'
    """
    samtools index $bam
    geneBody_coverage.py \\
        -i $bam \\
        -o ${prefix} \\
        -r $bed12
    mv log.txt ${prefix}.log.txt
    """
}

/*
 * STEP 5 - preseq analysis
 */
process preseq {
    tag "${bam_preseq.baseName - '.sorted'}"
    publishDir "${params.outdir}/05-preseq", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_preseq

    input:
    file bam_preseq

    output:
    file "${prefix}.ccurve.txt" into preseq_results

    script:
    prefix = bam_preseq.baseName - '_filteredAligned.sortedByCoord.out'
    """
    preseq lc_extrap -v -B $bam_preseq -o ${prefix}.ccurve.txt
    """
}


/*
 * STEP 6 Mark duplicates
 */
process markDuplicates {
    tag "${bam.baseName - '.sorted'}"
    publishDir "${params.outdir}/06-removeDuplicates/picard", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam from bam_markduplicates

    output:
    file "${prefix}.markDups.bam" into bam_md
    file "${prefix}.markDups_metrics.txt" into picard_results
    file "${prefix}.markDups.bam.bai"

    script:
    prefix = bam.baseName - '_filteredAligned.sortedByCoord.out'
//    markdup_java_options = (task.memory.toGiga() > 8) ? ${params.markdup_java_options} : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""

    """
    picard MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${prefix}.markDups.bam \\
        METRICS_FILE=${prefix}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${prefix}.markDups.bam
    """
}


/*
 * STEP 7 - dupRadar
 */
process dupradar {
    label 'low_memory'
    tag "${bam_md.baseName - '.sorted.markDups'}"
    publishDir "${params.outdir}/06-removeDuplicates", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve_mqc.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }
    cpus '20'
	penv 'openmp'

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam_md
    file gtf from gtf_dupradar.collect()

    output:
    file "*.{pdf,txt}" into dupradar_results

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    def dupradar_direction = 0
    if (forward_stranded && !unstranded) {
        dupradar_direction = 1
    } else if (reverse_stranded && !unstranded){
        dupradar_direction = 2
    }
    def paired = params.singleEnd ? 'single' :  'paired'
    """
    dupRadar.r $bam_md $gtf $dupradar_direction $paired 20
    """
}


/*
 * STEP 8 Feature counts
 */
process featureCounts {
    label 'low_memory'
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/07-featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()
    file biotypes_header from ch_biotypes_header.collect()

    output:
    file "${sample_name}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${sample_name}_gene.featureCounts.txt.summary" into featureCounts_logs
    file "${sample_name}_biotype_counts*mqc.{txt,tsv}" into featureCounts_biotype

    script:
    def featureCounts_direction = 0
    def extraAttributes = params.fcExtraAttributes ? "--extraAttributes ${params.fcExtraAttributes}" : ''
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    // Try to get real sample name
    sample_name = bam_featurecounts.baseName - '_filteredAligned.sortedByCoord.out'
    """
    featureCounts -a $gtf -g ${params.fcGroupFeatures} -o ${sample_name}_gene.featureCounts.txt $extraAttributes -p -s $featureCounts_direction $bam_featurecounts
    featureCounts -a $gtf -g ${params.fcGroupFeaturesType} -o ${sample_name}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    cut -f 1,7 ${sample_name}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${sample_name}_biotype_counts_mqc.txt
    mqc_features_stat.py ${sample_name}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${sample_name}_biotype_counts_gs_mqc.tsv
    """
}

/*
 * STEP 9 - Merge featurecounts
 */
process merge_featureCounts {
    tag "${input_files[0].baseName - '.sorted'}"
    publishDir "${params.outdir}/07-featureCounts", mode: 'copy'

    input:
    file input_files from featureCounts_to_merge.collect()

    output:
    file 'merged_gene_counts.txt'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand,gene_name"'
    """
    $merge $input_files | csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | sed 's/.markDups.bam//g' > merged_gene_counts.txt
    """
}


/*
 * STEP 10 - stringtie FPKM
 */
process stringtieFPKM {
    tag "${bam_stringtieFPKM.baseName - '.sorted'}"
    publishDir "${params.outdir}/08-stringtieFPKM", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
            else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
            else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
            else "$filename"
        }

    input:
    file bam_stringtieFPKM
    file gtf from gtf_stringtieFPKM.collect()

    output:
    file "${prefix}_transcripts.gtf"
    file "${prefix}.gene_abund.txt"
    file "${prefix}.cov_refs.gtf"
    file ".command.log" into stringtie_log
    file "${prefix}_ballgown"

    script:
    prefix = bam_stringtieFPKM.baseName - '_filteredAligned.sortedByCoord.out'
    def st_direction = ''
    if (forward_stranded && !unstranded){
        st_direction = "--fr"
    } else if (reverse_stranded && !unstranded){
        st_direction = "--rf"
    }
    """
    stringtie $bam_stringtieFPKM \\
        $st_direction \\
        -o ${prefix}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${prefix}.gene_abund.txt \\
        -C ${prefix}.cov_refs.gtf \\
        -e \\
        -b ${prefix}_ballgown
    """
}

/*
 * STEP 11 - edgeR MDS and heatmap
 */
process sample_correlation {
    label 'low_memory'
    tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
    publishDir "${params.outdir}/09-sample_correlation", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_edger

    input:
    file input_files from geneCounts.collect()
    val num_bams from bam_count.count()
    file mdsplot_header from ch_mdsplot_header
    file heatmap_header from ch_heatmap_header

    output:
    file "*.{txt,pdf,csv}" into sample_correlation_results

    when:
    num_bams > 2 && (!params.sampleLevel)

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    """
    edgeR_heatmap_MDS.r $input_files
    cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
    mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
    cat $heatmap_header log2CPM_sample_distances_mqc.csv >> tmp_file
    mv tmp_file log2CPM_sample_distances_mqc.csv
    """
}

/*
 * STEP 12 MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/99-stats/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file multiqc_config from multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimommatic/*') from trimmomatic_results.collect()
    file ('trimommatic/*') from trimmomatic_fastqc_reports.collect()
    file ('alignment/*') from alignment_logs.collect()
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('rseqc/*') from genebody_coverage_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect()
    file ('featureCounts_biotype/*') from featureCounts_biotype.collect()
    file ('stringtie/stringtie_log*') from stringtie_log.collect()
    file ('sample_correlation_results/*') from sample_correlation_results.collect().ifEmpty([]) // If the Edge-R is not run create an Empty array

    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'

    """
    multiqc -d . --config $multiqc_config
    """
}


/*
 * STEP 13 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/99-stats/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


workflow.onComplete {
    log.info "BU-ISCIII - Pipeline complete"
}
