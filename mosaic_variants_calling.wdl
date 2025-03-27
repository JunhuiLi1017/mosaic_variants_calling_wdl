version 1.0

struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
    Int boot_disk_size
}

workflow sortBamWorkflow {
    input {
        File input_bam_cram
        File input_bai
        File ref_fasta
        File ref_fai
        File ref_dict
        File dbsnp138
        File dbsnp138_idx
        File known_indels
        File known_indels_idx
        File mills_1000G
        File mills_1000G_idx
        Int threads = 8          # Shared across tasks
        Boolean pcr_based = true  # PCR-based deduplication flag
        # Task-specific memory overrides
        Int mem_gb_CramToBam = 8
        Int mem_gb_SortBam = 12
        Int mem_gb_MarkDup = 16
        Int mem_gb_BaseRecal = 16
        Int mem_gb_ApplyBQSR = 12

        # mutect2
        File intervals
        File pon
        File pon_idx
        Int scatter_count
        File gnomad
        File gnomad_idx
        String? m2_extra_args
        String? m2_extra_filtering_args
        String? split_intervals_extra_args
        Boolean? compress_vcfs

        # mosaic
        Float ? low_af
        Float ? high_af
        File segdup_bed #SegDup_and_clustered.GRCh38.bed
        File allrepeats_bed #allrepeats_forindel.hg38.bed
        File umap_file
        File feature_extraction_script
        File avail_acess_1000g
        String cov
        File prediction_script
        File model_snv
        File model_ins
        File model_del
        Int? emergency_extra_disk
        File humandb_refGene
        File humandb_refGeneMrna
        File humandb_dbnsfp42a
        File gnomad_2_1_1_0_001
        File gnomad_2_1_1_0


        # depth of coverate
        Int mem_gb_m2 = 16
        Int mem_gb_depthcov = 16
        Int mem_gb_feature = 16
        Int mem_gb_mosaic_tier = 16

        # runtime
        String gatk_docker
        File? gatk_override
        String basic_bash_docker = "ubuntu:16.04"

        Int? preemptible
        Int? max_retries
        Int small_task_cpu = 2
        Int small_task_mem = 4
        Int small_task_disk = 20
        Int boot_disk_size = 12
        Int learn_read_orientation_mem = 8000
        Int filter_alignment_artifacts_mem = 9000

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk

        # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
        # Large is for Bams/WGS vcfs
        # Small is for metrics/other vcfs
        Float large_input_to_output_multiplier = 3
        Float small_input_to_output_multiplier = 2.0
        Float cram_to_bam_multiplier = 7.0
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])
    Int emergency_extra_disk_or_default = select_first([emergency_extra_disk,10])
    Boolean compress = select_first([compress_vcfs, false])

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
    Int cram_size = ceil(size(input_bam_cram, "GB") + size(input_bai, "GB"))
    Int gnomad_vcf_size = if defined(gnomad) then ceil(size(gnomad, "GB")) else 0
    Int bound_size = ceil(size(dbsnp138, "GB") + size(dbsnp138_idx, "GB") + size(known_indels, "GB") + size(known_indels_idx, "GB") + size(mills_1000G, "GB") + size(mills_1000G_idx, "GB"))
    Int humandb_size = ceil(size(humandb_refGene, "GB") + size(humandb_refGeneMrna, "GB") + size(humandb_dbnsfp42a, "GB"))
    Int gnomad_size = ceil(size(gnomad_2_1_1_0_001, "GB") + size(gnomad_2_1_1_0, "GB"))
    
    # If no tar is provided, the task downloads one from broads ftp server
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + gatk_override_size + emergency_extra_disk_or_default

    # logic about output file names -- these are the names *without* .vcf extensions
    String output_basename = basename(basename(input_bam_cram, ".bam"),".cram")  #hacky way to strip either .bam or .cram
    String unfiltered_name = output_basename + "-unfiltered"
    String filtered_name = output_basename + "-filtered"

    String output_vcf_name = output_basename + ".vcf"

    Int cram_to_bam_disk = ceil(cram_size * cram_to_bam_multiplier + ref_size)

    Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries_or_default, "preemptible": preemptible_or_default, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": small_task_disk + disk_pad, "boot_disk_size": boot_disk_size}

    # Step 1: Convert CRAM to BAM
    call CramToBam {
        input:
            input_bam_cram = input_bam_cram,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            mem_gb = mem_gb_CramToBam,
            threads = threads,
            disk_space_gb = cram_to_bam_disk
    }

    # Calculate disk space for sorting
    Int sort_bam_disk = ceil(size(CramToBam.output_bam, "GB") * 3) + emergency_extra_disk_or_default

    # Step 2: Sort BAM with Sambamba
    call SortBamWithSambamba {
        input:
            unsort_coor_bam = CramToBam.output_bam,
            mem_gb = mem_gb_SortBam,
            threads = threads,
            disk_space_gb = sort_bam_disk
    }

    # Calculate disk space for removing duplicates
    Int mark_dup_bam_disk = ceil(size(SortBamWithSambamba.sort_bam, "GB") * 3) + emergency_extra_disk_or_default

    # Step 3: Remove Duplicates
    call MarkDup {
        input:
            input_sort_bam = SortBamWithSambamba.sort_bam,
            pcr_based = pcr_based,
            mem_gb = mem_gb_MarkDup,
            threads = threads,
            disk_space_gb = mark_dup_bam_disk
    }
    # Calculate disk space for base recalibration
    Int recal_disk = ceil(size(MarkDup.sort_mark_bam, "GB") * 3) + ref_size + bound_size + emergency_extra_disk_or_default

    # Step 4: Base Recalibration
    call BaseRecalibrator {
        input:
            mark_dup_bam = MarkDup.sort_mark_bam,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            dbsnp138 = dbsnp138,
            dbsnp138_idx = dbsnp138_idx,
            known_indels = known_indels,
            known_indels_idx = known_indels_idx,
            mills_1000G = mills_1000G,
            mills_1000G_idx = mills_1000G_idx,
            mem_gb = mem_gb_BaseRecal,
            threads = threads,
            disk_space_gb = recal_disk
    }

    # Calculate disk space for apply bqsr 
    Int apply_bqsr_disk = ceil(size(MarkDup.sort_mark_bam, "GB") * 3) + ref_size + emergency_extra_disk_or_default

    # Step 5: Apply BQSR
    call ApplyBQSR {
        input:
            mark_dup_bam = MarkDup.sort_mark_bam,
            recal_table = BaseRecalibrator.recal_table,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            mem_gb = mem_gb_ApplyBQSR,
            threads = threads,
            disk_space_gb = apply_bqsr_disk
    }

    # Step 6: split intervals
    call SplitIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = standard_runtime
    }
    
    #TODO: do we need to change this disk size now that NIO is always going to happen (for the google backend only)
    Int bqsr_bam_size = ceil(size(ApplyBQSR.recal_bam, "GB")) + ceil(size(ApplyBQSR.recal_bai, "GB"))
    Int m2_output_size = bqsr_bam_size / scatter_count
    Int m2_per_scatter_size = bqsr_bam_size + ref_size + gnomad_vcf_size + m2_output_size + disk_pad

    # Step 7: mutect2
    scatter (subintervals in SplitIntervals.interval_files ) {
        call M2 {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                bqsr_bam = ApplyBQSR.recal_bam,
                bqsr_bai = ApplyBQSR.recal_bai,
                pon = pon,
                pon_idx = pon_idx,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                preemptible = preemptible,
                max_retries = max_retries,
                m2_extra_args = m2_extra_args,
                compress = compress,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                disk_space = m2_per_scatter_size,
                mem_gb=mem_gb_m2
        }
    }

    call MergeVCFs {
        input:
            input_vcfs = M2.unfiltered_vcf,
            input_vcf_indices = M2.unfiltered_vcf_idx,
            output_name = unfiltered_name,
            compress = compress,
            runtime_params = standard_runtime
    }

    call MergeStats { input: stats = M2.stats, runtime_params = standard_runtime }
    
    call Filter {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            intervals = intervals,
            unfiltered_vcf = MergeVCFs.merged_vcf,
            unfiltered_vcf_idx = MergeVCFs.merged_vcf_idx,
            output_name = filtered_name,
            compress = compress,
            mutect_stats = MergeStats.merged_stats,
            m2_extra_filtering_args = m2_extra_filtering_args,
            runtime_params = standard_runtime,
            disk_space = ceil(size(MergeVCFs.merged_vcf, "GB") * small_input_to_output_multiplier) + disk_pad
    }

    # Conditionally set low_af based on whether pcr_based is true or not
    Float default_low_af = select_first([low_af,0.02])
    Float final_low_af = if (defined(pcr_based) && pcr_based) then default_low_af else 0.03
    Float final_high_af = select_first([high_af, 0.4])
    call MT2InitialFilter {
        input:
            input_vcf = Filter.filtered_vcf,
            low_af = final_low_af,
            high_af = final_high_af,
            runtime_params = standard_runtime,
            disk_space = ceil(size(Filter.filtered_vcf, "GB") * small_input_to_output_multiplier) + 5
    }

    call repeat_filter {
        input:
            mt2pon_bed = MT2InitialFilter.mt2_initial_filter_bed,
            segdup_bed = segdup_bed,
            runtime_params = standard_runtime,
            disk_space = ceil((size(MT2InitialFilter.mt2_initial_filter_bed, "GB") + size(segdup_bed, "GB")) * small_input_to_output_multiplier) + 5
    }

    call annovar_formatter {
        input:
            input_bed = repeat_filter.filter_segdup_bed,
            runtime_params = standard_runtime,
            disk_space = ceil(size(repeat_filter.filter_segdup_bed, "GB") * small_input_to_output_multiplier) + 5
    }

    call MAF0_extraction_SNV {
        input:
            input_list = annovar_formatter.output_variant,
            runtime_params = standard_runtime,
            disk_space = ceil(size(annovar_formatter.output_variant, "GB") * small_input_to_output_multiplier) + 5
    }

    call MAF0_extraction_INS {
        input:
            input_list = annovar_formatter.output_variant,
            allrepeats_bed = allrepeats_bed,
            runtime_params = standard_runtime,
            disk_space = ceil((size(annovar_formatter.output_variant, "GB") + size(allrepeats_bed, "GB"))* small_input_to_output_multiplier) + 5
    }

    call MAF0_extraction_DEL {
        input:
            input_list = annovar_formatter.output_variant,
            allrepeats_bed = allrepeats_bed,
            runtime_params = standard_runtime,
            disk_space = ceil((size(annovar_formatter.output_variant, "GB") + size(allrepeats_bed, "GB"))* small_input_to_output_multiplier) + 5
    }

    call FeatureExtraction {
        input:
            bed_file_snv = MAF0_extraction_SNV.output_bed_snv,
            bed_file_ins = MAF0_extraction_INS.output_bed_ins,
            bed_file_del = MAF0_extraction_DEL.output_bed_del,
            umap_file = umap_file,
            bqsr_bam = ApplyBQSR.recal_bam,
            bqsr_bai = ApplyBQSR.recal_bai,
            reference_fasta = ref_fasta,
            feature_extraction_script = feature_extraction_script,
            runtime_params = standard_runtime,
            mem_gb = mem_gb_feature,
            disk_space = ceil((size(MAF0_extraction_SNV.output_bed_snv, "GB") * 3 + size(umap_file, "GB") + size(ApplyBQSR.recal_bam, "GB") + ref_size) * small_input_to_output_multiplier) + 5
    }

    call Filter1000G {
        input:
            feature_snv = FeatureExtraction.output_feature_snv,
            feature_ins = FeatureExtraction.output_feature_ins,
            feature_del = FeatureExtraction.output_feature_del,
            avail_acess_1000g = avail_acess_1000g,
            runtime_params = standard_runtime,
            disk_space = ceil((size(FeatureExtraction.output_feature_snv, "GB") * 3 + size(avail_acess_1000g, "GB")) * small_input_to_output_multiplier) + 5
    }

    call PredictionMosaic {
        input:
            cov=cov,
            prediction_script = prediction_script,
            model_snv = model_snv,
            model_del = model_del,
            model_ins = model_ins,
            features_snv = Filter1000G.output_filtered_snv,
            features_del = Filter1000G.output_filtered_del,
            features_ins = Filter1000G.output_filtered_ins,
            runtime_params = standard_runtime,
            disk_space = ceil((size(Filter1000G.output_filtered_snv, "GB") * 3) * small_input_to_output_multiplier) + 5
    }

    # Set multipliers and padding

    call AnnotationMosaic {
        input:
            predictions_snv = PredictionMosaic.predictions_snv,
            predictions_ins = PredictionMosaic.predictions_ins,
            predictions_del = PredictionMosaic.predictions_del,
            humandb_refGene = humandb_refGene,
            humandb_refGeneMrna = humandb_refGeneMrna,
            humandb_dbnsfp42a = humandb_dbnsfp42a,
            runtime_params = standard_runtime,
            disk_space = ceil((size(PredictionMosaic.predictions_snv, "GB") * 3 + humandb_size) * small_input_to_output_multiplier) + 5
    }

    call ExtractSubvcf {
        input:
            vcf_gz = Filter.filtered_vcf,
            predictions_snv = AnnotationMosaic.input_predictions_snv,
            predictions_ins = AnnotationMosaic.input_predictions_ins,
            predictions_del = AnnotationMosaic.input_predictions_del,
            runtime_params = standard_runtime,
            disk_space = ceil((size(AnnotationMosaic.input_predictions_snv, "GB") * 3 + size(Filter.filtered_vcf, "GB")) * small_input_to_output_multiplier) + 5
    }

    call MosaicTier {
        input:
            vcf = ExtractSubvcf.sub_vcf_snv,
            anno = AnnotationMosaic.output_predictions_snv,
            gnomad_2_1_1_0_001 = gnomad_2_1_1_0_001,
            gnomad_2_1_1_0 = gnomad_2_1_1_0,
            runtime_params = standard_runtime,
            mem_gb = mem_gb_mosaic_tier,
            disk_space = ceil(gnomad_size * small_input_to_output_multiplier) + 5
    }

    scatter(subintervals in SplitIntervals.interval_files) {
        call DepthOfCoverage {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                bqsr_bam = ApplyBQSR.recal_bam,
                bqsr_bam_idx = ApplyBQSR.recal_bai,
                gatk_docker = gatk_docker,
                disk_space = m2_per_scatter_size,
                mem_gb=mem_gb_depthcov
        }
    }

    call MergeCovStatistics {
        input:
            statistics = DepthOfCoverage.depth_of_cov_statistics,
            runtime_params = standard_runtime,
            disk_space = ceil(size(DepthOfCoverage.depth_of_cov_statistics, "GB") * small_input_to_output_multiplier) + 5
    }

    output {
        File recal_bam_flag = ApplyBQSR.recal_bamflagstat
        File merged_vcf = MergeVCFs.merged_vcf
        File merged_vcf_idx = MergeVCFs.merged_vcf_idx
        File filtered_vcf = Filter.filtered_vcf
        File filtering_stats = Filter.filtering_stats
        File mutect_stats = MergeStats.merged_stats

        File mt2_initial_filter_bed = MT2InitialFilter.mt2_initial_filter_bed
        File mt2_segdup_bed = repeat_filter.filter_segdup_bed
        File var_anno= annovar_formatter.output_variant
        File extraction_snv = MAF0_extraction_SNV.output_bed_snv
        File extraction_ins = MAF0_extraction_INS.output_bed_ins
        File extraction_del = MAF0_extraction_DEL.output_bed_del
        File feature_snv = FeatureExtraction.output_feature_snv
        File feature_ins = FeatureExtraction.output_feature_ins
        File feature_del = FeatureExtraction.output_feature_del
        File filter_1000g_snv = Filter1000G.output_filtered_snv
        File filter_1000g_ins = Filter1000G.output_filtered_ins
        File filter_1000g_del = Filter1000G.output_filtered_del
        File mosaic_snv = PredictionMosaic.predictions_snv
        File mosaic_ins = PredictionMosaic.predictions_ins
        File mosaic_del = PredictionMosaic.predictions_del
        File in_anno_snv = AnnotationMosaic.input_predictions_snv
        File in_anno_ins = AnnotationMosaic.input_predictions_ins
        File in_anno_del = AnnotationMosaic.input_predictions_del
        File anno_snv = AnnotationMosaic.output_predictions_snv
        File anno_ins = AnnotationMosaic.output_predictions_ins
        File anno_del = AnnotationMosaic.output_predictions_del
        File snv_sub_vcf = ExtractSubvcf.sub_vcf_snv
        File snv_subvcf_bed = ExtractSubvcf.subvcf_bed_snv
        File snv_geno_dp_af = ExtractSubvcf.geno_dp_af_snv
        File ins_sub_vcf = ExtractSubvcf.sub_vcf_ins
        File ins_subvcf_bed = ExtractSubvcf.subvcf_bed_ins
        File ins_geno_dp_af = ExtractSubvcf.geno_dp_af_ins
        File del_sub_vcf = ExtractSubvcf.sub_vcf_del
        File del_subvcf_bed = ExtractSubvcf.subvcf_bed_del
        File del_geno_dp_af = ExtractSubvcf.geno_dp_af_del
        File vcf3 = MosaicTier.vcf3
        File vcf4 = MosaicTier.vcf4
        File anno3 = MosaicTier.anno3
        File anno4 = MosaicTier.anno4
        File final_tier = MosaicTier.tier
        File stat_all = MergeCovStatistics.stat_all
    }
}

# Task Definitions (unchanged except for fixes)
task CramToBam {
    input {
        File input_bam_cram
        File ref_fasta
        File ref_fai
        Int disk_space_gb = 50
        Int mem_gb = 8
        Int threads = 8
    }

    String sample = sub(basename(input_bam_cram), "\\..*$", "")
    Int mem_per_thread_gb = floor(mem_gb / threads)

    command <<<
        samtools view \
            -b \
            -@ ~{threads} \
            -m ~{mem_per_thread_gb}G \
            -T ~{ref_fasta} \
            -o ~{sample}.bam \
            ~{input_bam_cram}
    >>>

    output {
        File output_bam = "~{sample}.bam"
    }

    runtime {
        docker: "junhuili/terra_samtools:1.20"
        memory: "~{mem_gb} GB"
        cpu: "~{threads}"
        disks: "local-disk ~{disk_space_gb} HDD"
    }
}

task SortBamWithSambamba {
    input {
        File unsort_coor_bam
        Int mem_gb = 16
        Int threads = 8
        Int disk_space_gb = 100
    }

    String sample = sub(basename(unsort_coor_bam), "\\..*$", "")
    Int mem_per_thread_gb = floor(mem_gb / threads)

    command <<<
        sambamba sort \
            -t ~{threads} \
            -m ~{mem_per_thread_gb}G \
            -o ~{sample}.sortbycoor.bam \
            ~{unsort_coor_bam}
    >>>

    output {
        File sort_bam = "~{sample}.sortbycoor.bam"
        File sort_bai = "~{sample}.sortbycoor.bam.bai"
    }

    runtime {
        docker: "junhuili/sambamba:1.0.1"
        memory: "~{mem_gb} GB"
        cpu: "~{threads}"
        disks: "local-disk ~{disk_space_gb} HDD"
    }
}

task MarkDup {
    input {
        File input_sort_bam
        Boolean pcr_based = true
        Int threads = 8
        Int mem_gb = 96
        Int disk_space_gb = 100
    }

    String sample = sub(basename(input_sort_bam), "\\..*$", "")
    String dedup_params = if pcr_based then "--REMOVE_DUPLICATES true" else "--REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"

    command <<<
        java -jar /usr/local/bin/picard.jar MarkDuplicates \
            --INPUT ~{input_sort_bam} \
            --METRICS_FILE ~{sample}.sort.rmdup.matrix \
            --OUTPUT ~{sample}.sort.rmdup.bam \
            ~{dedup_params}
        samtools index ~{sample}.sort.rmdup.bam
    >>>

    output {
        File sort_mark_bam = "~{sample}.sort.rmdup.bam"
        File sort_mark_bai = "~{sample}.sort.rmdup.bam.bai"
        File rmdupMatrix = "~{sample}.sort.rmdup.matrix"
    }

    runtime {
        docker: "junhuili/picard-samtools:2.23.3"
        memory: "~{mem_gb} GB"
        cpu: "~{threads}"
        disks: "local-disk ~{disk_space_gb} HDD"
    }
}

task BaseRecalibrator {
    input {
        File mark_dup_bam
        File ref_fasta
        File ref_fai
        File ref_dict
        File dbsnp138
        File dbsnp138_idx
        File known_indels
        File known_indels_idx
        File mills_1000G
        File mills_1000G_idx
        Int threads = 8
        Int mem_gb = 16
        Int disk_space_gb = 5
    }

    String sample = sub(basename(mark_dup_bam), "\\..*$", "")

    command <<<
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xms10g" BaseRecalibrator \
            -I ~{mark_dup_bam} \
            -O ~{sample}.recal_data.table \
            -R ~{ref_fasta} \
            --known-sites ~{dbsnp138} \
            --known-sites ~{known_indels} \
            --known-sites ~{mills_1000G}
    >>>

    output {
        File recal_table = "~{sample}.recal_data.table"
    }

    runtime {
        docker: "broadinstitute/gatk:4.1.8.1"  # Fixed "broadinstitute/gatk:4.1.8.1" placeholder
        memory: "~{mem_gb} GB"
        cpu: "~{threads}"
        disks: "local-disk ~{disk_space_gb} HDD"
    }
}

task ApplyBQSR {
    input {
        File mark_dup_bam
        File recal_table
        File ref_fasta
        File ref_fai
        File ref_dict
        Int threads = 8
        Int mem_gb = 12
        Int disk_space_gb = 200
    }

    String sample = sub(basename(mark_dup_bam), "\\..*$", "")

    command <<<
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xms10g -XX:ParallelGCThreads=~{threads}" ApplyBQSR \
            -I ~{mark_dup_bam} \
            -O ~{sample}.bam \
            -R ~{ref_fasta} \
            --bqsr-recal-file ~{recal_table}
        samtools flagstat ~{sample}.bam > ~{sample}.txt
    >>>

    output {
        File recal_bam = "~{sample}.bam"
        File recal_bai = "~{sample}.bai"
        File recal_bamflagstat = "~{sample}.txt"
    }

    runtime {
        docker: "broadinstitute/gatk:4.1.8.1"  # Fixed "broadinstitute/gatk:4.1.8.1" placeholder
        memory: "~{mem_gb} GB"
        cpu: "~{threads}"
        disks: "local-disk ~{disk_space_gb} HDD"
    }
}

task SplitIntervals {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count
      String? split_intervals_extra_args

      # runtime
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"

        mkdir interval-files
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" SplitIntervals \
            -R ~{ref_fasta} \
            ~{"-L " + intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File bqsr_bam
      File bqsr_bai
      File? pon
      File? pon_idx
      File? gnomad
      File? gnomad_idx
      String? m2_extra_args
      Boolean compress

      File? gatk_override

      # runtime
      String gatk_docker
      Int? mem_gb
      Int? preemptible
      Int? max_retries
      Int? disk_space
      Int? cpu
      Boolean use_ssd = false
    }

    String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1000 else 8500
    Int command_mem = machine_mem - 500

    parameter_meta{
      intervals: {localization_optional: true}
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      bqsr_bam: {localization_optional: true}
      bqsr_bai: {localization_optional: true}
      pon: {localization_optional: true}
      pon_idx: {localization_optional: true}
      gnomad: {localization_optional: true}
      gnomad_idx: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        echo "" > normal_name.txt

        gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{bqsr_bam} -O tumor_name.txt -encode
        tumor_command_line="-I ~{bqsr_bam} -tumor `cat tumor_name.txt`"

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
            -R ~{ref_fasta} \
            $tumor_command_line \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            -O "~{output_vcf}" \
            ~{m2_extra_args}
        m2_exit_code=$?
        set +e

        # the script only fails if Mutect2 itself fails
        exit $m2_exit_code
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        File stats = "~{output_stats}"
    }
}

task MergeVCFs {
    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      String output_name
      Boolean compress
      Runtime runtime_params
    }
    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task MergeStats {
    input {
      Array[File]+ stats
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}


        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_stats = "merged.stats"
    }
}

task Filter {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File unfiltered_vcf
      File unfiltered_vcf_idx
      Boolean compress
      File? mutect_stats
      String? m2_extra_filtering_args
      String output_name
      Runtime runtime_params
      Int? disk_space
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{ref_fasta} \
            -O ~{output_vcf} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats filtering.stats \
            ~{m2_extra_filtering_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "filtering.stats"
    }
}

task MT2InitialFilter {
    input {
        File input_vcf
        Float low_af 
        Float high_af
        Runtime runtime_params
        Int? disk_space
    }

    String sample = basename(basename(input_vcf,".gz"), "-filtered.vcf")
    String output_bed = sample + ".init.bed"

    command <<<
        set -e

        echo "Processing input VCF file..."
        echo ~{input_vcf}
        echo ~{low_af}
        echo ~{high_af}

        input_vcf=~{input_vcf}
        if [[ $input_vcf == *.gz ]]; then
            cat_cmd="zcat"
        else
            cat_cmd="cat"
        fi

        echo "Running filtering command..."

        $cat_cmd $input_vcf | grep -v "^#" | grep PASS | \
        gawk '{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){print $0}}' | \
        cut -f1,2,4,5,10 | sed "s/:/\t/g" | sed "s/,/\t/g" | \
        awk '$8>= ~{low_af} && $8< ~{high_af}' | grep -v '0|1' | grep -v '1|0' > temp1.bed

        $cat_cmd $input_vcf | grep -v "^#" | grep PASS | \
        gawk '{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){print $0}}' | \
        cut -f1,2,4,5,10 | sed "s/:/\t/g" | sed "s/,/\t/g" | \
        awk '$8>= 0.02 && $8< ~{high_af}' | grep -E "0\|1|1\|0" > temp2.bed

        cat temp1.bed temp2.bed | \
        cut -f 1-4,6-8 | awk -v sample="~{sample}" '{OFS="\t";print $1,$2-1,$2,$3,$4,sample,$5,$6,$7}' \
        > ~{output_bed}

        echo "Checking output file..."
        if [ ! -s ~{output_bed} ]; then
            echo "Error: Output file is empty. Check the previous commands for errors." >&2
            exit 1
        fi

        echo "Filtering complete. Output file: ~{output_bed}"
        set +e
    >>>

    output {
        File mt2_initial_filter_bed = "~{output_bed}"
    }

    runtime {
        docker: "junhuili/terra_init:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task repeat_filter {
    input {
        File mt2pon_bed
        File segdup_bed
        Runtime runtime_params
        Int? disk_space
    }

    # Extract sample name from the input file name
    String output_filename = sub(basename(mt2pon_bed), "init\\.bed$", "noSegDup\\.bed")

    command {
        subtractBed -a ${mt2pon_bed} -b ${segdup_bed} > ${output_filename}
    }

    output {
        File filter_segdup_bed = "${output_filename}"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task annovar_formatter {
    input {
        File input_bed
        Runtime runtime_params
        Int? disk_space
    }
    String output_bed = sub(basename(input_bed),"noSegDup\\.bed$","anno\\.list")

    command <<<
        cat ~{input_bed} | awk '{OFS="\t";len=length($4)-length($5);if(len<=0){print $1,$3,$3,$4,$5,$6}if(len>0){print $1,$3,$3+len,$4,$5,$6}}' > ~{output_bed}
        sed -i 's/chr//g' ~{output_bed}
    >>>

    output {
        File output_variant = "~{output_bed}"
    }

    runtime {
        docker: "ubuntu:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task MAF0_extraction_SNV {
    input {
        File input_list
        Runtime runtime_params
        Int? disk_space
    }
    String out_bed = sub(basename(input_list),"anno\\.list", "MAF0\\.SNV\\.chr\\.bed")

    command <<<
        cat ~{input_list} | awk '{OFS="\t";print $1,$2-1,$2,$4,$5,$6}' | awk 'length($4)==1 && length($5)==1' | awk 'BEGIN{OFS="\t"} $1="chr"$1' > ~{out_bed}
    >>>

    output {
        File output_bed_snv = "~{out_bed}"
    }

    runtime {
        docker: "ubuntu:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task MAF0_extraction_INS {
    input {
        File input_list
        File allrepeats_bed
        Runtime runtime_params
        Int? disk_space
    }
    String out_bed = sub(basename(input_list),"anno\\.list", "MAF0\\.INS\\.chr\\.bed")

    command <<<
        cat ~{input_list} | awk '{OFS="\t";print $1,$2-1,$2,$4,$5,$6}' | awk 'length($4) < length($5)' > temp_ins.bed
        bedtools subtract -a temp_ins.bed -b ~{allrepeats_bed} | awk 'BEGIN{OFS="\t"} $1="chr"$1' > ~{out_bed}
    >>>

    output {
        File output_bed_ins = "~{out_bed}"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task MAF0_extraction_DEL {
    input {
        File input_list
        File allrepeats_bed
        Runtime runtime_params
        Int? disk_space
    }
    String out_bed = sub(basename(input_list),"anno\\.list", "MAF0\\.DEL\\.chr\\.bed")

    command <<<
        cat ~{input_list} | awk '{OFS="\t";print $1,$2-1,$2,$4,$5,$6}' | awk 'length($4) > length($5)' > temp_del.bed
        bedtools subtract -a temp_del.bed -b ~{allrepeats_bed} | awk 'BEGIN{OFS="\t"} $1="chr"$1' > ~{out_bed}
    >>>

    output {
        File output_bed_del = "~{out_bed}"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task FeatureExtraction {
    input {
        File bed_file_snv
        File bed_file_ins
        File bed_file_del
        File umap_file
        File bqsr_bam
        File bqsr_bai
        File reference_fasta
        File feature_extraction_script
        Runtime runtime_params
        Int? disk_space
        Int? mem_gb
    }
    String sample_name=sub(basename(bed_file_snv),"\\.MAF0\\.SNV\\.chr\\.bed","")
    Int machine_mem = if defined(mem_gb) then mem_gb * 1000 else 8000

    command <<<
        set -e
        export PATH=/opt/conda/bin:$PATH
        source activate py3.7.1 || conda activate py3.7.1
        export PYTHONWARNINGS="ignore"
        
        # Create a directory for BAM files
        bam_dir=$(dirname ~{bqsr_bam})
        
        python ~{feature_extraction_script} \
            ~{bed_file_snv} \
            ~{sample_name}.SNV.features \
            $bam_dir \
            ~{reference_fasta} \
            ~{umap_file} \
            1 bam > ~{sample_name}.SNV.features.log 2>&1 &

        python ~{feature_extraction_script} \
            ~{bed_file_ins} \
            ~{sample_name}.INS.features \
            $bam_dir \
            ~{reference_fasta} \
            ~{umap_file} \
            1 bam > ~{sample_name}.INS.features.log 2>&1 &

        python ~{feature_extraction_script} \
            ~{bed_file_del} \
            ~{sample_name}.DEL.features \
            $bam_dir \
            ~{reference_fasta} \
            ~{umap_file} \
            1 bam > ~{sample_name}.DEL.features.log 2>&1 &
        wait
        
        conda deactivate
        set +e
    >>>

    output {
        File output_feature_snv = "~{sample_name}.SNV.features"
        File output_log_snv = "~{sample_name}.SNV.features.log"
        File output_feature_ins = "~{sample_name}.INS.features"
        File output_log_ins = "~{sample_name}.INS.features.log"
        File output_feature_del = "~{sample_name}.DEL.features"
        File output_log_del = "~{sample_name}.DEL.features.log"
    }

    runtime {
        docker: "junhuili/terra_py_tools:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + "MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task Filter1000G {
    input {
        File feature_snv
        File feature_ins
        File feature_del
        File avail_acess_1000g
        Runtime runtime_params
        Int? disk_space
    }
    String sample_name = sub(basename(feature_snv),"\\.SNV\\.features","")
    command <<<
        cat <(head -n 1 ~{feature_snv}) <(awk 'NR==FNR{c[$1,$2]=$0}NR!=FNR{if(c[$1,$2]){print $0}}' \
            <(bedtools intersect -a <(sed 's/~/\t/g' ~{feature_snv} | awk '{print $2"\t"$3"\t"$3}' | grep -v "conflict_num") \
            -b ~{avail_acess_1000g} -wa) \
            <(paste <(sed 's/~/\t/g' ~{feature_snv} | cut -f 2-3) ~{feature_snv}) | cut -f 3-) > ~{sample_name}.SNV.no1000g.features &

        cat <(head -n 1 ~{feature_ins}) <(awk 'NR==FNR{c[$1,$2]=$0}NR!=FNF{if(c[$1,$2]){print $0}}' \
            <(bedtools intersect -a <(sed 's/~/\t/g' ~{feature_ins} | awk '{print $2"\t"$3"\t"$3}' | grep -v "conflict_num") \
            -b ~{avail_acess_1000g} -wa) \
            <(paste <(sed 's/~/\t/g' ~{feature_ins} | cut -f 2-3) ~{feature_ins}) | cut -f 3-) > ~{sample_name}.INS.no1000g.features &

        cat <(head -n 1 ~{feature_del}) <(awk 'NR==FNR{c[$1,$2]=$0}NR!=FNF{if(c[$1,$2]){print $0}}' \
            <(bedtools intersect -a <(sed 's/~/\t/g' ~{feature_del} | awk '{print $2"\t"$3"\t"$3}' | grep -v "conflict_num") \
            -b ~{avail_acess_1000g} -wa) \
            <(paste <(sed 's/~/\t/g' ~{feature_del} | cut -f 2-3) ~{feature_del}) | cut -f 3-) > ~{sample_name}.DEL.no1000g.features &
        wait
    >>>

    output {
        File output_filtered_snv= "~{sample_name}.SNV.no1000g.features"
        File output_filtered_ins= "~{sample_name}.INS.no1000g.features"
        File output_filtered_del= "~{sample_name}.DEL.no1000g.features"
    }

    runtime {
        docker: "junhuili/terra_r:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
    
}

task PredictionMosaic {
    input {
        String cov="50X"
        File prediction_script
        File model_snv
        File model_ins
        File model_del
        File features_snv
        File features_del
        File features_ins
        Runtime runtime_params
        Int? disk_space
    }
    String sample_name = sub(basename(features_snv),"\\.SNV\\.no1000g\\.features","")

    command {
        Rscript ${prediction_script} \
        ${features_snv} \
        ${model_snv} \
        Refine \
        ${sample_name}.${cov}.SNV.predictions &

        Rscript ${prediction_script} \
        ${features_ins} \
        ${model_ins} \
        Refine \
        ${sample_name}.${cov}.INS.predictions &

        Rscript ${prediction_script} \
        ${features_del} \
        ${model_del} \
        Refine \
        ${sample_name}.${cov}.DEL.predictions &
        wait
    }

    output {
        File predictions_snv = "${sample_name}.${cov}.SNV.predictions"
        File predictions_ins = "${sample_name}.${cov}.INS.predictions"
        File predictions_del = "${sample_name}.${cov}.DEL.predictions"
    }

    runtime {
        docker: "junhuili/terra_r:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task AnnotationMosaic {
    input {
        File predictions_snv
        File predictions_ins
        File predictions_del
        File humandb_refGene
        File humandb_refGeneMrna
        File humandb_dbnsfp42a
        Runtime runtime_params
        Int? disk_space
    }

    String sample_name=sub(basename(predictions_snv),"\\.SNV\\.predictions","")
    command <<<
        set -e
        humandb=$(dirname ~{humandb_refGene})
        grep "mosaic" ~{predictions_snv} | cut -f 1,35 | awk 'BEGIN { FS = "~" } ; { print $2"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6 }' | tail -n +2 | awk '$6=="mosaic" {print $0}' > ~{sample_name}.SNV.input.predictions.txt
        cat <(grep "hap=3" ~{predictions_ins}) | awk 'BEGIN { FS = "~" } ; { print $2"\t"$3"\t"$3"\t"$4"\t"$5 }' | tail -n +2 > ~{sample_name}.INS.input.predictions.txt
        cat <(grep "hap=3" ~{predictions_del}) | awk 'BEGIN { FS = "~" } ; { print $2"\t"$3"\t"$3+length($4)-1"\t"$4"\t"$5 }' | tail -n +2 > ~{sample_name}.DEL.input.predictions.txt
    
        # Run ANNOVAR
        perl /usr/src/app/table_annovar.pl ~{sample_name}.SNV.input.predictions.txt $humandb -buildver hg38 -out ~{sample_name}.SNV.output.predictions -remove -protocol refGene,dbnsfp42a -operation g,f -nastring .
        perl /usr/src/app/table_annovar.pl ~{sample_name}.INS.input.predictions.txt $humandb -buildver hg38 -out ~{sample_name}.INS.output.predictions -remove -protocol refGene,dbnsfp42a -operation g,f -nastring .
        perl /usr/src/app/table_annovar.pl ~{sample_name}.DEL.input.predictions.txt $humandb -buildver hg38 -out ~{sample_name}.DEL.output.predictions -remove -protocol refGene,dbnsfp42a -operation g,f -nastring .
        echo "after anno"
        set +e
    >>>

    output {
        File input_predictions_snv = "${sample_name}.SNV.input.predictions.txt"
        File output_predictions_snv = "${sample_name}.SNV.output.predictions.hg38_multianno.txt"

        File input_predictions_ins = "${sample_name}.INS.input.predictions.txt"
        File output_predictions_ins = "${sample_name}.INS.output.predictions.hg38_multianno.txt"

        File input_predictions_del = "${sample_name}.DEL.input.predictions.txt"
        File output_predictions_del = "${sample_name}.DEL.output.predictions.hg38_multianno.txt"
    }

    runtime {
        docker: "junhuili/terra_perl_anno:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task ExtractSubvcf {
    input {
        File vcf_gz
        File predictions_snv
        File predictions_ins
        File predictions_del
        Runtime runtime_params
        Int? disk_space
    }
    String sample_name_snv = sub(basename(predictions_snv),"\\.input\\.predictions\\.txt","")
    String sample_name_ins = sub(basename(predictions_ins),"\\.input\\.predictions\\.txt","")
    String sample_name_del = sub(basename(predictions_del),"\\.input\\.predictions\\.txt","")
    command <<<
        cat ~{predictions_snv} | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' | sort -k1,1 -k2,2n > ~{sample_name_snv}.subvcf.bed
        if [[ "~{vcf_gz}" == *.gz ]]; then
            tabix ~{vcf_gz} -R ~{sample_name_snv}.subvcf.bed > ~{sample_name_snv}.sub.vcf.tmp
            cat <(zcat ~{vcf_gz} | grep "^#") <(cat ~{sample_name_snv}.sub.vcf.tmp) > ~{sample_name_snv}.sub.vcf
        else
            bedtools intersect -a ~{vcf_gz} -b ~{sample_name_snv}.subvcf.bed > ~{sample_name_snv}.sub.vcf.tmp
            cat <(cat ~{vcf_gz} | grep "^#") <(cat ~{sample_name_snv}.sub.vcf.tmp) > ~{sample_name_snv}.sub.vcf
        fi
        cut -f 1,2,4,5,10 ~{sample_name_snv}.sub.vcf.tmp | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}1' > ~{sample_name_snv}.Geno.DP.AF.txt

        cat ~{predictions_ins} | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' | sort -k1,1 -k2,2n > ~{sample_name_ins}.subvcf.bed
        if [[ "~{vcf_gz}" == *.gz ]]; then
            tabix ~{vcf_gz} -R ~{sample_name_ins}.subvcf.bed > ~{sample_name_ins}.sub.vcf.tmp
            cat <(zcat ~{vcf_gz} | grep "^#") <(cat ~{sample_name_ins}.sub.vcf.tmp) > ~{sample_name_ins}.sub.vcf
        else
            bedtools intersect -a ~{vcf_gz} -b ~{sample_name_ins}.subvcf.bed > ~{sample_name_ins}.sub.vcf.tmp
            cat <(cat ~{vcf_gz} | grep "^#") <(cat ~{sample_name_ins}.sub.vcf.tmp) > ~{sample_name_ins}.sub.vcf
        fi
        cut -f 1,2,4,5,10 ~{sample_name_ins}.sub.vcf.tmp | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}1' > ~{sample_name_ins}.Geno.DP.AF.txt

        cat ~{predictions_del} | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' | sort -k1,1 -k2,2n > ~{sample_name_del}.subvcf.bed
        if [[ "~{vcf_gz}" == *.gz ]]; then
            tabix ~{vcf_gz} -R ~{sample_name_del}.subvcf.bed > ~{sample_name_del}.sub.vcf.tmp
            cat <(zcat ~{vcf_gz} | grep "^#") <(cat ~{sample_name_del}.sub.vcf.tmp) > ~{sample_name_del}.sub.vcf
        else
            bedtools intersect -a ~{vcf_gz} -b ~{sample_name_del}.subvcf.bed > ~{sample_name_del}.sub.vcf.tmp
            cat <(cat ~{vcf_gz} | grep "^#") <(cat ~{sample_name_del}.sub.vcf.tmp) > ~{sample_name_del}.sub.vcf
        fi
        cut -f 1,2,4,5,10 ~{sample_name_del}.sub.vcf.tmp | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}1' > ~{sample_name_del}.Geno.DP.AF.txt
    >>>

    output {
        File sub_vcf_snv = "~{sample_name_snv}.sub.vcf"
        File subvcf_bed_snv = "~{sample_name_snv}.subvcf.bed"
        File geno_dp_af_snv = "~{sample_name_snv}.Geno.DP.AF.txt"
        File sub_vcf_ins = "~{sample_name_ins}.sub.vcf"
        File subvcf_bed_ins = "~{sample_name_ins}.subvcf.bed"
        File geno_dp_af_ins = "~{sample_name_ins}.Geno.DP.AF.txt"
        File sub_vcf_del = "~{sample_name_del}.sub.vcf"
        File subvcf_bed_del = "~{sample_name_del}.subvcf.bed"
        File geno_dp_af_del = "~{sample_name_del}.Geno.DP.AF.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task MosaicTier {
    input {
        File vcf
        File anno
        File gnomad_2_1_1_0_001
        File gnomad_2_1_1_0
        Runtime runtime_params
        Int? disk_space
        Int? mem_gb
    }
    String sample_name = sub(basename(vcf),"\\.vcf","")
    Int machine_mem = if defined(mem_gb) then mem_gb * 1000 else 10500

    command <<<
        awk 'NR==FNR{c[$1,$2,$4,$5]++;next} !(c[$1,$2,$4,$5]) {print $0}' ~{gnomad_2_1_1_0_001} ~{vcf} > ~{sample_name}.MFv4.vcf
        awk 'NR==FNR{c[$1,$2,$4,$5]++;next} !(c[$1,$2,$4,$5]) {print $0}' ~{gnomad_2_1_1_0} ~{vcf} > ~{sample_name}.MFv3.vcf
        cat <(head -n 1 ~{anno}) <(awk 'NR==FNR{{c[$1,$2,$4,$5]++;next}} c[$1,$2,$4,$5]>0{{print $0}}' ~{sample_name}.MFv4.vcf ~{anno}) > ~{sample_name}.MFv4.txt
        cat <(head -n 1 ~{anno}) <(awk 'NR==FNR{{c[$1,$2,$4,$5]++;next}} c[$1,$2,$4,$5]>0{{print $0}}' ~{sample_name}.MFv3.vcf ~{anno}) > ~{sample_name}.MFv3.txt
        awk 'BEGIN{OFS="\t"}NR==FNR{c[$1,$2,$4,$5]=$0} NR!=FNR{if (c[$2,$3,$5,$6]) {print "tier2\t"$0} else {print "tier3\t"$0}}' ~{sample_name}.MFv4.txt \
        <(awk 'BEGIN{OFS="\t"}NR==FNR{c[$1,$2,$4,$5]=$0} NR!=FNR{if (c[$1,$2,$4,$5]) {print "tier1\t"$0} else {print "tier3\t"$0}}' ~{sample_name}.MFv3.txt \
        ~{anno}) | sed 's/tier2\ttier1/tier1/g' | sed 's/tier2\ttier3/tier2/g' | sed 's/tier3\ttier3/tier3/g' |sed 's/tier3\ttier1/tier1/g' | sed 's/tier1\tChr/Tier\tChr/' \
        > ~{sample_name}.tier.txt
    >>>

    output {
        File vcf4 = "~{sample_name}.MFv4.vcf"
        File vcf3 = "~{sample_name}.MFv3.vcf"
        File anno4 = "~{sample_name}.MFv4.txt"
        File anno3 = "~{sample_name}.MFv3.txt"
        File tier = "~{sample_name}.tier.txt"
    }

    runtime {
        docker: "ubuntu:latest"
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task DepthOfCoverage {
    input {
        File intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File bqsr_bam
        File bqsr_bam_idx
        # runtime
        String gatk_docker
        Int? mem_gb
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Boolean use_ssd = false
        File? gatk_override
    }
    String name = sub(basename(bqsr_bam),"\\.bam","")
    Int machine_mem = if defined(mem_gb) then mem_gb * 1000 else 16500
    Int command_mem = machine_mem - 500

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx~{command_mem}m" DepthOfCoverage \
        -R ~{ref_fasta} \
        -I ~{bqsr_bam} \
        -L ~{intervals} \
        -O ~{name}
        set +e
    }
    output {
        File depth_of_cov_statistics = "~{name}.sample_statistics"
    }
    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 3])
        cpu: select_first([cpu, 2])
    }
}

task MergeCovStatistics {
        input {
            Array[File] statistics
            Runtime runtime_params
            Int? disk_space
        }
        
        command <<<
            file_list=$(printf "%s " ~{sep=' ' statistics})

            awk -F"," 'NR==1 && FNR==1 {header=$0; next} {for(i=2; i<=NF; i++) sum[i] += $i} 
            END {
                print header
                printf "%s", $1
                for(i=2; i<=NF; i++) printf ",%d", sum[i]
                print ""
            }' $file_list > output.sample_statistics.txt
        >>>
        output {
            File stat_all = "output.sample_statistics.txt"
        }

        runtime {
            docker: "ubuntu:latest"
            bootDiskSizeGb: runtime_params.boot_disk_size
            memory: runtime_params.machine_mem + " MB"
            disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
            preemptible: runtime_params.preemptible
            maxRetries: runtime_params.max_retries
            cpu: runtime_params.cpu
        }
}