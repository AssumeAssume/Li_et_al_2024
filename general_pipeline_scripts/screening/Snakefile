#vim: set syntax=python

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

"""
A CRISPR/Cas9 analysis workflow using MAGeCK, FastQC and VISPR.
"""


configfile: "config.yaml"


import sys
import yaml
from mageck_vispr import (postprocess_config, vispr_config, get_fastq,
                          annotation_available, get_counts, design_available,
                          need_annotate_bed_with_lfc,get_sample_name,
                          need_run_rra_in_mle,rra_treatment_string,rra_control_string,
                          get_norm_method, COMBAT_SCRIPT_PATH)


postprocess_config(config)

rule all:
    input:
        expand("results/{experiment}.vispr.yaml", experiment=config["experiments"]),
        (expand("annotation/{experiment}.rra.{nsample}_vs_{day0}.sgrnas.bed",experiment=config["experiments"],nsample=get_sample_name(config),day0=config["day0label"]) if ( ("day0label" in config) and annotation_available(config)) else [])

if "samples" in config:
    rule fastqc:
        input:
            lambda wildcards: config["replicates"][wildcards.replicate]
        output:
            directory("results/qc/{replicate}")
        log:
            "logs/fastqc/{replicate}.log"
        shell:
            "mkdir -p {output}; rm -rf {output}/*; "
            "fastqc -f fastq --extract -o {output} {input} 2> {log}"


    if "paired" in config:
        rule fastqc_paired:
            input:
                lambda wildcards: config["paired_rep"][wildcards.replicate]
            output:
                directory("results/qc/{replicate}_R2")
            log:
                "logs/fastqc/{replicate}_R2.log"
            shell:
                "mkdir -p {output}; rm -rf {output}/*; "
                "fastqc -f fastq --extract -o {output} {input} 2> {log}"


    if "adapter" in config["sgrnas"]:
        rule cutadapt:
            input:
                lambda wildcards: config["replicates"][wildcards.replicate]
            output:
                "results/trimmed_reads/{replicate}.fastq"
            shell:
                "cutadapt -a {config[sgrnas][adapter]} {input} > {output}"


    rule mageck_count:
        input:
            fastqs=[get_fastq(rep, config) for rep in config["replicates"]],
            library=config["library"]
        output:
            "results/count/all.count.txt",
            "results/count/all.count_normalized.txt",
            "results/count/all.countsummary.txt"
        params:
            labels=",".join(config["samples"].keys()),
            norm=get_norm_method(config),
            fastqs=" ".join(
                ",".join(get_fastq(rep, config) for rep in replicates)
                for replicates in config["samples"].values()),
            pairedfastqs=(
                "" if not "paired" in config
                else "--fastq-2 "+(" ".join(
                ",".join(config["paired_rep"][rep] for rep in replicates)
                for replicates in config["paired"].values()))),
            countpair=str(
                "" if not "countpair" in config
                else "--count-pair "+str(config["countpair"])),
            prefix="results/count/all",
            day0=(
                "" if not "day0label" in config
                else "--day0-label "+config["day0label"]),
            controlsg=(
                "" if not "control_sgrna" in config
                else "--control-sgrna "+config["control_sgrna"])
        log:
            "logs/mageck/count/all.log"
        shell:
            "mageck count --output-prefix {params.prefix} "
            "--norm-method {params.norm} "
            "--list-seq {input.library} "
            "--fastq {params.fastqs} --sample-label {params.labels} "
            "{params.pairedfastqs} "
            "{params.countpair} "
            "{params.controlsg} "
            "{params.day0} --trim-5 {config[sgrnas][trim-5]} 2> {log}"


if "counts" in config:
    rule mageck_qc:
        input:
            counts=get_counts(config),
        output:
            "results/count/all.count_normalized.txt",
            "results/count/all.countsummary.txt",
            "results/count/all_countsummary.R",
            "results/count/all_countsummary.Rnw"
        params:
            prefix="results/count/all",
            norm=get_norm_method(config),
            day0=(
                "" if not "day0label" in config
                else "--day0-label "+config["day0label"]),
            controlsg=(
                "" if not "control_sgrna" in config
                else "--control-sgrna "+config["control_sgrna"])
        log:
            "logs/mageck/count/all.log"
        shell:
            "mageck count --output-prefix {params.prefix} {params.day0} "
            "--norm-method {params.norm} "
            "{params.controlsg} "
            "--count-table {input.counts} 2> {log}"


if "library" in config:
    rule annotate_sgrnas:
        input:
            config["library"]
        output:
            "annotation/sgrnas.bed"
        params:
            annotation_file=("--annotation-table "+config["sgrnas"]["annotation-sgrna-file"] if ("annotation-sgrna-file" in config["sgrnas"] ) else " "),
            annotation_folder=("--annotation-table-folder "+config["sgrnas"]["annotation-sgrna-folder"] if ("annotation-sgrna-folder" in config["sgrnas"] ) else " ")
        log:
            "logs/annotation/sgrnas.log"
        shell:
            "mkdir -p annotation; "
            "mageck-vispr annotate-library {input} "
            "{params.annotation_file} "
            "{params.annotation_folder} "
            "--sgrna-len {config[sgrnas][len]} --assembly {config[assembly]} "
            "> {output} 2> {log}"


if "batchmatrix" in config:
    rule remove_batch:
        input:
            counts=config.get("counts", "results/count/all.count_normalized.txt"),
            batchmatrix=config["batchmatrix"]
        output:
            "results/count/all.count.batchcorrected.txt"
        log:
            "logs/combat.log"
        script:
            COMBAT_SCRIPT_PATH


ruleorder: mageck_mle > mageck_rra

rule mageck_rra:
    input:
        counts=get_counts(config)
    output:
        genesummary="results/test/{experiment}.gene_summary.txt",
        sgrnasummary="results/test/{experiment}.sgrna_summary.txt",
        indivoutput=(expand("results/test/{{experiment}}.{nsample}_vs_{day0}.sgrna_summary.txt",nsample=get_sample_name(config),day0=config["day0label"]) if "day0label" in config else [])
    params:
        prefix="results/test/{experiment}",
        #treatment=lambda wildcards: ",".join(config["experiments"][wildcards.experiment]["treatment"]),
        #control=lambda wildcards: ",".join(config["experiments"][wildcards.experiment]["control"]),
        treatment=lambda wildcards: rra_treatment_string(wildcards,config),
        control=lambda wildcards: rra_control_string(wildcards, config),
        day0=(
            "" if not "day0label" in config
            else "--day0-label "+config["day0label"]),
        norm=get_norm_method(config),
        controlsg=(
            "" if not "control_sgrna" in config
            else "--control-sgrna "+config["control_sgrna"]),
        cnv_correct=(
            "" if not config["correct_cnv"] 
            else "--cnv-norm "+config["cnv_norm"]+" --cell-line "+config["cnv_cell_line"]),
        additionalparameter=(
            "" if not "additional_rra_parameter" in config
            else " "+config["additional_rra_parameter"])
    log:
        "logs/mageck/test/{experiment}.log"
    shell:
        "mageck test --norm-method {params.norm} "
        "--output-prefix {params.prefix} "
        "--count-table {input} "
        "{params.controlsg} "
        #"--treatment-id {params.treatment} "
        #"--control-id {params.control} "
        "{params.control} "
        "{params.treatment} "
        "{params.day0} "
        "{params.cnv_correct} "
        "{params.additionalparameter} "
        "2> {log} "

if need_annotate_bed_with_lfc(config):
    rule annotate_sgrna_after_rra:
        input:
            library=config["library"],
            annotation="annotation/sgrnas.bed",
            input_sgrna_summary=(
                "results/test/{experiment}.rra.{nsample}_vs_{day0}.sgrna_summary.txt" if "day0label" in config 
                else "results/test/{experiment}.sgrna_summary.txt")
        output:
            ("annotation/{experiment}.rra.{nsample}_vs_{day0}.sgrnas.bed" if "day0label" in config 
            else "annotation/{experiment}.sgrnas.bed")
        params:
            input_column_string="LFC"
        log:
            ("logs/annotation/{experiment}.rra.{nsample}_vs_{day0}.sgrnas.log" if "day0label" in config 
            else "logs/annotation/{experiment}.sgrnas.log")
        shell:
            # mageck-vispr annotate-library yusa_library.csv --sgrna-len 19 --assembly mm10 --annotation-table /src/exome_scan/sgrna_annotation_mm10_exome_19bp.txt.bz2 --bedvalue results/test/ESC-mle.rra.esc1_vs_plasmid.sgrna_summary.txt --bedvalue-column KFC
            "mageck-vispr annotate-library {input.library} "
            "--sgrna-len {config[sgrnas][len]} --assembly {config[assembly]} "
            "--bedvalue {input.input_sgrna_summary} "
            "--bedvalue-column {params.input_column_string} "
            "> {output} 2> {log}"

rule mageck_mle:
    input:
        counts=get_counts(config),
        has_designmatrix=lambda wildcards: config["experiments"][wildcards.experiment]["designmatrix"],
        annotation="annotation/sgrnas.bed" if annotation_available(config) else []
        #cnv_profile=config["cnv_norm"] if config["correct_cnv"] else []
    output:
        "results/test/{experiment}.gene_summary.txt",
        "results/test/{experiment}.sgrna_summary.txt"
    params:
        prefix="results/test/{experiment}",
        efficiency=(
            "" if not annotation_available(config) or not config["sgrnas"]["annotate-sgrna-efficiency"]
            else "--sgrna-eff-name-column 3 --sgrna-eff-score-column 4 --sgrna-efficiency annotation/sgrnas.bed"),
        update_efficiency=(
            "" if not config["sgrnas"].get("update-efficiency", False)
            else "--update-efficiency"),
        norm=get_norm_method(config),
        designmatrix=(lambda wildcards: "" if not design_available(config)
            else "--design-matrix " +  config["experiments"][wildcards.experiment]["designmatrix"]),
        day0=(
            "" if not "day0label" in config
            else "--day0-label "+config["day0label"]),
        controlsg=(
            "" if not "control_sgrna" in config
            else "--control-sgrna "+config["control_sgrna"]),
        threads=(
            "" if not "threads" in config
            else "--threads "+str(config["threads"])),
        cnv_correct=(
            "" if not config["correct_cnv"] 
            else "--cnv-norm "+config["cnv_norm"]),
        additionalparameter=(
            "" if not "additional_mle_parameter" in config
            else " "+config["additional_mle_parameter"])
    log:
        "logs/mageck/test/{experiment}.log"
    shell:
        "mageck mle --norm-method {params.norm} "
        "--output-prefix {params.prefix} {params.efficiency} --genes-var 0 "
        "{params.update_efficiency} --count-table {input.counts} "
        "{params.cnv_correct} "
        "{params.threads} {params.controlsg} {params.designmatrix} {params.day0} "
        "{params.additionalparameter} "
        "2> {log}"


rule vispr:
    input:
        "annotation/sgrnas.bed" if annotation_available(config) else [],
        # lfcbed="annotation/{experiment}.sgrnas.bed" if need_annotate_bed_with_lfc(config) else [],
        results="results/test/{experiment}.gene_summary.txt",
        results2=(lambda wildcards: "results/test/{experiment}.rra.gene_summary.txt" if need_run_rra_in_mle(wildcards,config)  else []),
        #results2="results/test/{experiment}.rra.gene_summary.txt",
        sgrna_results="results/test/{experiment}.sgrna_summary.txt",
        counts=(get_counts(config, normalized=True) if "samples" in config else "results/count/all.count_normalized.txt" ),
        mapstats=("results/count/all.countsummary.txt" if "samples" in config else []),
        fastqc=(expand("results/qc/{replicate}", replicate=config["replicates"]) if "samples" in config else []),
        pairedfastqc=(expand("results/qc/{replicate}_R2", replicate=config["paired_rep"]) if "paired" in config else [])
    output:
        "results/{experiment}.vispr.yaml"
    run:
        vispr_config(input, output, wildcards, config)

