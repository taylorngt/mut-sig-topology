INDEX_SUFFIXES = ["1", "2", "3", "4", "5", "6", "7", "8"]
NUMBERED_CHROMOSOMES = ["1", "2", "3", "4", "5"]
COMPRESSION_LEVEL = "8"
ALL_SAMPLEIDs = ["Sin1", "Son1", "San8", "San11", "Sin1_2"]
TSTV_PARAMETERS = "w1000000af1"
CtoT_PARAMETERS = "w1000000af1"
SIG11_PARAMETERS = "w2000000af0"
GENE_PANELS = ["nicotinamide_gene_list"]


rule hisat_reference_indexing:
    input:
        reference_genome="reference_genome/a_thaliana.fa"
    output:
        ref_indeces=expand("reference_indeces/a_thaliana.{index_num}.ht2", index_num=INDEX_SUFFIXES)
    shell:
        "hisat2-build "
        "-f {input.reference_genome} reference_indeces/a_thaliana "

rule samtools_reference_indexing:
    input:
        reference_genome="reference_genome/a_thaliana.fa"
    output:
        reference_index="reference_genome/a_thaliana.fa.fai"
    shell:
        "samtools faidx {input.reference_genome}"

rule samtools_reference_dictionary_construction:
    input:
        reference_genome="reference_genome/a_thaliana.fa"
    output:
        reference_dictionary="reference_genome/a_thaliana.dict"
    shell:
        'samtools dict {input.reference_genome} -o {output.reference_dictionary}'

rule alignment_and_compression:
    input:
        reference_genome="reference_genome/a_thaliana.fa",
        ref_indeces=expand("reference_indeces/a_thaliana.{index_num}.ht2", index_num=INDEX_SUFFIXES),
        forward_reads="sample_reads/{sampleID}/{sampleID}_1.fastq",
        reverse_reads="sample_reads/{sampleID}/{sampleID}_2.fastq"
    output:
        compressed_alignment="compressed_alignments/{sampleID}.bam"
    threads: workflow.cores
    shell:
        "hisat2 -p {threads} -x reference_indeces/a_thaliana "
        "-q -1 {input.forward_reads} -2 {input.reverse_reads} "
        "--rg-id {wildcards.sampleID} --rg SM:{wildcards.sampleID} "
        "| "
        "samtools view -@ {threads} -b "
        "- > {output.compressed_alignment}"



rule sort_alignment:
    input:
        compressed_alignment="compressed_alignments/{sampleID}.bam"
    output:
        sorted_alignment="sorted_alignments/{sampleID}_sorted.bam"
    threads: workflow.cores
    params:
        compress = COMPRESSION_LEVEL
    shell:
        "samtools sort -@ {threads} -l {params.compress} {input.compressed_alignment} "
        "-o {output.sorted_alignment}"



rule alignment_indexing:
    input:
        sorted_alignment="sorted_alignments/{sampleID}_sorted.bam"
    output:
        alignment_index="sorted_alignments/{sampleID}_sorted.bam.bai"
    threads: workflow.cores
    shell:
        "samtools index --threads {threads} -b {input.sorted_alignment}"



rule mutect_calling:
    input:
        sample_alignment="sorted_alignments/{sampleID}_sorted.bam",
        sample_alignment_index="sorted_alignments/{sampleID}_sorted.bam.bai",
        wildtype_alignment="sorted_alignments/Ws2_sorted.bam",
        wildtype_alignment_index="sorted_alignments/Ws2_sorted.bam.bai",
        reference_genome="reference_genome/a_thaliana.fa",
        reference_index="reference_genome/a_thaliana.fa.fai",
        reference_dictionary="reference_genome/a_thaliana.dict"
    output:
        mutect_variant_calls="mutect_variant_calls/{sampleID}.vcf.gz"
    shell:
        "gatk Mutect2 "
        "-R {input.reference_genome} "
        "-I {input.sample_alignment} "
        "-I {input.wildtype_alignment} "
        "--tumor-sample {wildcards.sampleID} "
        "--normal-sample Ws2 "
        "-O {output.mutect_variant_calls}"



rule variant_artifact_labeling:
    input:
        mutect_variant_calls="mutect_variant_calls/{sampleID}.vcf.gz",
        reference_genome="reference_genome/a_thaliana.fa"
    output:
        artifact_labeled_variants="labeled_mutect_variant_calls/{sampleID}.vcf.gz"
    shell:
        "gatk FilterMutectCalls "
        "-R {input.reference_genome} "
        "-V {input.mutect_variant_calls} "
        "-O {output.artifact_labeled_variants}"



rule artifact_filtering_and_formatting:
    input:
        artifact_labeled_variants="labeled_mutect_variant_calls/{sampleID}.vcf.gz",
        chromosome_name_translation="reference_genome/TAIR10_chromosome_name_translation.txt"
    output:
        artifact_filtered_variants="filtered_mutect_variant_calls/{sampleID}.vcf.gz"
    shell:
        "bcftools view "
        "{input.artifact_labeled_variants} "
        "-f PASS | "
        "bcftools annotate "
        "--rename-chrs {input.chromosome_name_translation} "
        "-Oz -o {output.artifact_filtered_variants}"

rule variant_call_indexing:
    input:
        artifact_filtered_variants="filtered_mutect_variant_calls/{sampleID}.vcf.gz"
    output:
        variant_call_index = "filtered_mutect_variant_calls/{sampleID}.vcf.gz.csi"
    shell:
        'bcftools index {input.artifact_filtered_variants}'


rule filtered_variant_unzipping:
    input:
        compressed_variants="filtered_mutect_variant_calls/{sampleID}.vcf.gz"
    output:
        unzipped_variants=temp("filtered_mutect_variant_calls/{sampleID}.vcf")
    shell:
        "gunzip -k {input.compressed_variants}"
    


rule topology_generation:
    input: 
        artifact_filtered_variants = "filtered_mutect_variant_calls/{sampleID}.vcf",
        reference_genome = "reference_genome/a_thaliana.fa"
    params:
        chromosomes = NUMBERED_CHROMOSOMES,
        sampleIDs = ALL_SAMPLEIDs
    wildcard_constraints:
        afModifier = "(0|1)",
        method = "(CtoT|TSTV|Sig11)"
    output:
        topology_plot = "topology_plots/{sampleID}/{sampleID}{method}_w{windowSize}af{afModifier}_stair.png",
        distribution_plot = "topology_plots/{sampleID}/{sampleID}{method}_w{windowSize}af{afModifier}_hist.png",
        peak_regions = "topology_peak_regions/{sampleID}/{sampleID}{method}_w{windowSize}af{afModifier}_peaks.txt"
    script:
        'scripts/snpTopology.py'




rule peak_region_concatenation:
    input:
        TSTV_peaks= expand("topology_peak_regions/{{sampleID}}/{{sampleID}}TSTV_{TSTV_PARAMETERS}_peaks.txt", TSTV_PARAMETERS=TSTV_PARAMETERS),
        CtoT_peaks= expand("topology_peak_regions/{{sampleID}}/{{sampleID}}CtoT_{CtoT_PARAMETERS}_peaks.txt", CtoT_PARAMETERS=CtoT_PARAMETERS),
        Sig11_peaks= expand("topology_peak_regions/{{sampleID}}/{{sampleID}}Sig11_{SIG11_PARAMETERS}_peaks.txt", SIG11_PARAMETERS=SIG11_PARAMETERS)
    output:
        combined_peaks= "topology_peak_regions/{sampleID}/{sampleID}_combined_peaks.txt"
    shell:
        'cat {input.TSTV_peaks} {input.CtoT_peaks} {input.Sig11_peaks} '
        '> {output.combined_peaks}'



rule variant_region_filtering:
    input:
        combined_peaks= "topology_peak_regions/{sampleID}/{sampleID}_combined_peaks.txt",
        mutect_variant_calls = "filtered_mutect_variant_calls/{sampleID}.vcf.gz",
        mutect_variant_call_index = "filtered_mutect_variant_calls/{sampleID}.vcf.gz.csi"
    output:
        region_filtered_variant_calls= "region_filtered_variant_calls/{sampleID}.vcf.gz"
    shell:
        "bcftools view -R {input.combined_peaks} "
        "-Oz -o {output.region_filtered_variant_calls} "
        "{input.mutect_variant_calls}"



rule VCF_header_parsing:
    input:
        region_filtered_variant_calls= "region_filtered_variant_calls/{sampleID}.vcf.gz"
    output:
        temp_header = temp('temp_headers/{sampleID}_header.vcf')
    shell:
        'bcftools head {input.region_filtered_variant_calls} | '
        'head -n -1 > {output.temp_header}'

rule region_filtered_variant_unzipping:
    input:
        compressed_variants="region_filtered_variant_calls/{sampleID}.vcf.gz"
    output:
        unzipped_variants=temp("region_filtered_variant_calls/{sampleID}.vcf")
    shell:
        "gunzip -k {input.compressed_variants}"


rule variant_annotation:
    input:
        region_filtered_variants= "region_filtered_variant_calls/{sampleID}.vcf",
        nicotinamide_gene_panel= "gene_panels/nicotinamide_gene_panel.txt",      
        circadian_gene_panel = "gene_panels/circadian_gene_panel.txt"
    output:
        headless_tiered_variants= temp("tiered_variant_calls/{sampleID}_headerless.vcf")
        # gene_list= "gene_lists/{sampleID}_gene_list.txt"
    script:
        'scripts/VEP_API.py'



rule VCF_header_replacement:
    input:
        temp_header = 'temp_headers/{sampleID}_header.vcf',
        headerless_tiered_variants= "tiered_variant_calls/{sampleID}_headerless.vcf"
    output:
        tiered_variants= "tiered_variant_calls/{sampleID}_tiered.vcf"
    shell:
        'cat {input.temp_header} {input.headerless_tiered_variants} '
        '> {output.tiered_variants}'








#Off-shoot direcitves not currently used in main pipeline
rule mutation_matrix_generation:
    input:
        mutect_variant_calls = "filtered_mutect_variant_calls/{sampleID}.vcf",
        reference_genome = "reference_genome/a_thaliana.fa"   
    output:
        matrix_plot="mutation_count_matrices/allSamplesCatalogPlots/{sampleID}_mut_cat.png",
        mutation_matrix="mutation_count_matrices/{sampleID}_mut_matrix.csv"
    script:
        'scripts/MutationMatrixGenerator.py'



