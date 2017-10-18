from pathlib import Path
GATK_JAR_PTH='/gscmnt/gc7201/dinglab/medseq/lwang/conda_root/envs/gatk_wxs_cnv/share/gatk4-4.0b5-0/gatk-package-4.beta.5-local.jar'
GENOME_REF_FA = '/gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170530.fa'
EXOME_TARGET_TSV = '/gscmnt/gc2737/ding/sample_info/SureSelectV5UTR/S04380219_Padded.merged.new.bed.gatkTarget.tsv'

WXS_NORMAL_BAMS_PTH = '/gscmnt/gc2737/ding/sample_info/wxs_normal_abspath.txt'
WXS_TUMOR_BAMS_PTH = '/gscmnt/gc2737/ding/sample_info/wxs_tumor_abspath.txt'

# Read all the bam paths from the txt, they are not really the BAM files in use
_WXS_NORMAL_BAMS = []
_WXS_TUMOR_BAMS = []
with open(WXS_NORMAL_BAMS_PTH) as f:
    _WXS_NORMAL_BAMS.extend(Path(l.strip()) for l in f)
with open(WXS_TUMOR_BAMS_PTH) as f:
    _WXS_TUMOR_BAMS.extend(Path(l.strip()) for l in f)

# Extract the SRR ids for all tumor/normal bam files
NORMAL_SRR_IDS = [pth.name.split('.', 1)[0] for pth in _WXS_NORMAL_BAMS]
TUMOR_SRR_IDS = [pth.name.split('.', 1)[0] for pth in _WXS_TUMOR_BAMS]

# Create a SRRxxxx id mapping to the bam file
# TODO: clean up this part of ugly code
srr_id_to_bam_pth = {}
for bam_pth in _WXS_NORMAL_BAMS:
    srr_id = bam_pth.name.split('.', 1)[0]
    # if the alternative bam file exists, we use them
    alt_bam_pth = bam_pth.with_name(f'{bam_pth.name}.2.bam')
    if alt_bam_pth.exists():
        srr_id_to_bam_pth[srr_id] = alt_bam_pth
    elif bam_pth.exists():
        srr_id_to_bam_pth[srr_id] = bam_pth
    else:
        print(f"Normal BAM {bam_pth} does not exist. Skipped!")
        NORMAL_SRR_IDS.remove(srr_id)

for bam_pth in _WXS_TUMOR_BAMS:
    srr_id = bam_pth.name.split('.', 1)[0]
    # if the alternative bam file exists, we use them
    alt_bam_pth = bam_pth.with_name(f'{bam_pth.name}.2.bam')
    if alt_bam_pth.exists():
        srr_id_to_bam_pth[srr_id] = alt_bam_pth
    elif bam_pth.exists():
        srr_id_to_bam_pth[srr_id] = bam_pth
    else:
        print(f"Tumor BAM {bam_pth} does not exist. Skipped!")
        TUMOR_SRR_IDS.remove(srr_id)


def find_bam_pth(wildcards):
    return str(srr_id_to_bam_pth[wildcards.srr_id])


rule calculate_target_coverage:
    """Caluculate the coverage within the exome seq target region."""
    input:
        bam=find_bam_pth,
        exome_target_tsv=EXOME_TARGET_TSV,
    output: "coverage/{srr_id}.pcov"
    resources: mem=4000  # unit: MB
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'CalculateTargetCoverage '
        '-I {input.bam} -T {input.exome_target_tsv} '
        '-transform PCOV -groupBy SAMPLE -targetInfo FULL '
        '-O {output}'


rule make_panel_of_normal_pcov_list:
    input: expand('coverage/{srr_id}.pcov', srr_id=NORMAL_SRR_IDS)
    output: 'panel_of_normal/normal_pcov_files.list'
    run:
        with open(output[0], 'w') as f:
            for srr_pth in input:
                print(srr_pth, file=f)


rule combine_normal_samples_coverage:
    input: 'panel_of_normals/normal_pcov_files.list'
    output: 'panel_of_normals/combined_normal.pcov'
    resources: mem=4000  # 4000MB
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'CombineReadCounts -inputList {input} -O {output}'


rule create_panel_of_normals:
    input: 'panel_of_normals/combined_normal.pcov'
    output: 'panel_of_normals/normals.pon'
    resources: mem=32000  # 32GB
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'CreatePanelOfNormals -I {input} -O {output} '
        '-noQC --disableSpark '
        '--minimumTargetFactorPercentileThreshold 5'


rule normalize_somatic_read_counts:
    input:
        tumor_pcov="coverage/{srr_id}.pcov",
        pon="panel_of_normals/normals.pon",
    output:
        ptn="somatic_CNV/{srr_id}/processed/{srr_id}.ptn.tsv",
        tn="somatic_CNV/{srr_id}/processed/{srr_id}.tn.tsv",
    resources:
        # So the memory will start with 8GB, 16GB, and 32GB
        # For each attempt
        mem=lambda wildcards, attempt: 4000 * 2**attempt
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'NormalizeSomaticReadCounts '
        '-I {input.tumor_pcov} -PON {input.pon} '
        '-PTN {output.ptn} -TN {output.tn}'


rule perform_segmentation:
    input: tn="somatic_CNV/{srr_id}/processed/{srr_id}.tn.tsv"
    resources: mem=4000
    output: seg="somatic_CNV/{srr_id}/processed/{srr_id}.seg"
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'PerformSegmentation -TN {input.tn} -O {output.seg} -LOG'


rule plot_segmented_copy_ratio:
    input:
        ptn="somatic_CNV/{srr_id}/processed/{srr_id}.ptn.tsv",
        tn="somatic_CNV/{srr_id}/processed/{srr_id}.tn.tsv",
        seg="somatic_CNV/{srr_id}/processed/{srr_id}.seg",
        seq_dict=GENOME_REF_FA + '.dict',
    params:
        output_dir="somatic_CNV/{srr_id}/plots",
        sample_prefix="{srr_id}",
    resources: mem=4000
    output:
        expand(
            "somatic_CNV/{{srr_id}}/plots/"
            "{{srr_id}}{filename}",
            filename=[
                "_FullGenome.png", "_Before_After.png", "_Before_After_CR_Lim_4.png",
                "_preQc.txt", "_postQc.txt", "_dQc.txt", "_scaled_dQc.txt"
            ])
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'PlotSegmentedCopyRatio -PTN {input.ptn} -TN {input.tn} -S {input.seg} '
        '-SD {input.seq_dict} -O {params.output_dir} -pre {params.sample_prefix} -LOG'


rule call_segments:
    input:
        tn="somatic_CNV/{srr_id}/processed/{srr_id}.tn.tsv",
        seg="somatic_CNV/{srr_id}/processed/{srr_id}.seg",
    output:
        cnv_call="somatic_CNV/{srr_id}/{srr_id}.cnv"
    resources: mem=4000
    shell:
        'java -jar -Xmx{resources.mem}m {GATK_JAR_PTH} '
        'CallSegments -TN {input.tn} -S {input.seg} -O {output.cnv_call}'


rule all:
    input:
        expand("somatic_CNV/{srr_id}/{srr_id}.cnv", srr_id=TUMOR_SRR_IDS),
        expand("somatic_CNV/{srr_id}/plots/{srr_id}_FullGenome.png", srr_id=TUMOR_SRR_IDS),
    params: LSF="-N"
    shell:
        'python --version'
