# mosaic_variants_calling_wdl

Mosaic Variants Calling Pipeline Description
The mosaic_variants_calling pipeline is a scalable and reproducible bioinformatics workflow implemented in WDL (Workflow Description Language) for detecting mosaic variants from whole-genome sequencing (WGS) data. Designed to process paired-end FASTQ files aligned to a reference genome (e.g., hg38), this pipeline automates the identification of low-frequency somatic mutations, such as single nucleotide variants (SNVs), with high accuracy and efficiency. It is hosted on GitHub at github.com/JunhuiLi1017/mosaic_variants_calling_wdl.

The workflow integrates several industry-standard tools into a cohesive pipeline:

Read Alignment: Utilizes BWA-MEM to align paired-end FASTQ reads (R1 and R2) to a user-specified reference genome, producing a BAM file.
Duplicate Marking: Employs Picard MarkDuplicates to identify and flag PCR duplicates, ensuring accurate variant calling by reducing artifacts.
Base Quality Score Recalibration (BQSR): Applies GATK’s BaseRecalibrator and ApplyBQSR to improve base quality scores using known variant sites (e.g., dbSNP, Mills, and 1000G datasets), enhancing call reliability.
Variant Calling: Leverages GATK’s Mutect2 in tumor-only mode to detect mosaic SNVs, followed by filtering with FilterMutectCalls to refine variant calls based on statistical models.
Output: Generates a filtered VCF file containing high-confidence mosaic variants, alongside intermediate BAM files and quality metrics for downstream analysis.
Key features include:

Modularity: Structured as a WDL workflow with configurable inputs (e.g., FASTQ files, reference genome, known variant databases), allowing adaptation to diverse datasets.
Scalability: Designed for execution on cloud-based or high-performance computing platforms (e.g., Cromwell), supporting large-scale genomic studies.
Reproducibility: Encapsulates all dependencies and parameters in a single workflow, ensuring consistent results across runs.
This pipeline reflects my expertise in developing end-to-end bioinformatics solutions for genomic variant analysis, with applications in cancer research, rare disease studies, and population genetics.
