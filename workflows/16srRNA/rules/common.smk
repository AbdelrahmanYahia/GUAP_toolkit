
common_rules = config["common_rules"]
include: f'{common_rules}/common.smk'



EXTT = config["ext"]
EXT = EXTT.replace(".gz","")
DOWNSTREAM = config["downstream"]
USE_THREADS_QIIME = ALL_THREADS
USE_MEM = int((MEM/2)*1000)


if ALL_THREADS <= 4:
    USE_THREADS = 1
else:
    USE_THREADS = 4

def get_dada_RAW_input(wildcards):
    inputs = []
    inputs.extend(expand(
        f"{PATH}/{{sample}}_{RS}{{R}}{TAIL}.{EXT}",
        R = [1, 2],
        sample = samples_names
    ))
    return inputs

def get_final_output(wildcards):
    final_input = []

    final_input.extend([
        "QIIME2/visualization/table.qzv",
        "QIIME2/visualization/bar_plot.qzv",
        "QIIME2/visualization/rep-seqs.qzv",
        "QIIME2/classify/taxonomy.qza",
        "QIIME2/visualization/taxonomy.qzv",
        "QIIME2/unrooted-tree.nwk",
        "QIIME2/rooted-tree.nwk"
        ]) 

    if config["deblur"]:
        final_input.extend([
            "QIIME2/deblur/deblur-stats.qza",
            "QIIME2/visualization/deblur-stats.qzv",
            "QIIME2/rep-seqs.qza",
            "QIIME2/table.qza"
        ]) 

    elif config["use_QIIME2"]:
        final_input.extend([
            "QIIME2/table-stats.qza",
            "QIIME2/rep-seqs.qza",
            "QIIME2/table.qza"
        ])     
    else:
        final_input.extend([
            "DADA2/rep-seqs.fna",
            "DADA2/seqtab-nochim.txt",
            "QIIME2/rep-seqs.qza",
            "QIIME2/table.qza"
        ])

    if config["export_figs"]:
        final_input.extend([
            directory("Phyloseq_Figures/")
        ])

    if config["downstream"]:
        final_input.extend([
            "QIIME2/visualization/alpha-rarefaction.qzv",
            "QIIME2/visualization/shannon_vector-group-significance.qzv",
            directory("QIIME2/core-metrics-results/"),
            "ps.rds",
            "QIIME2/visualization/faith-pd-group-significance.qzv",
            "QIIME2/visualization/shannon_vector-group-significance.qzv",
            "QIIME2/visualization/evenness-group-significance.qzv",
            "QIIME2/visualization/unweighted-unifrac-condition-significance.qzv",
            "QIIME2/visualization/unweighted-unifrac-emperor.qzv",
            "QIIME2/visualization/bray-curtis-emperor.qzv"
        ])
    return final_input


def get_R_pattern(n):
    if config["trimmomatic"] or config["remove_primers"]:
        return(f"_{RS}{n}.{EXT}")
    else:
        if config["naming_pattern"] == "illumina":
            return(f"*S[0-9]+{lane}_{RS}{n}{TAIL}.{EXT}")
        else:
            return(f"*{lane}_{RS}{n}{TAIL}.{EXT}")


def get_analysis_input_dir():
    if config["trimmomatic"]:
        if config["remove_primers"]:
            return("cutadapt")
        else:
            return("trimmomatic")
    elif config["remove_primers"]: 
        return("cutadapt")
    else:
        return(PATH)


def get_rename_prefix():
    if config["deblur"]:
        return("deblur")
    elif config["use_QIIME2"]: 
        return("dada")
    else:
        return("r")


def get_analysis_input():
    if config["trimmomatic"]:
        if config["remove_primers"]:
            return("multiqc/multiqc_report.html")
        else:
            return("multiqc/multiqc_report.html")
    elif config["remove_primers"]: 
        return("multiqc/multiqc_report.html")
    
    elif config["skip_QC"] :
        return(unpack(get_dada_RAW_input))
    
    else:
        return("multiqc/multiqc_report.html")


def get_classifeir_input(wildcards):
    if config["classifier"] == "dada":
        return("DADA2/seqtab.nochim.RDS")
    else:
        return("QIIME2/rep-seqs.qza")


def get_multiqc_input(wildcards):
    final_input = []
    final_input.extend(expand(
        f"QC/{{sample}}_{RS}{{R}}_fastqc.{{ext}}",
        ext = ["zip", "html"],
        R = [1, 2],
        sample = samples
    ))

    if config["trimmomatic"]:
        final_input.extend(expand(
                f"trimmomatic/{{sample}}_{RS}{{R}}.{EXT}",
                R = [1, 2],
                sample = samples
        ))
    
    if config["remove_primers"]:
        final_input.extend(expand(
                f"cutadapt/{{sample}}_{RS}{{R}}.{EXT}",
                R = [1, 2],
                sample = samples
        ))    
    
    return final_input

AID = get_analysis_input_dir()
R1_pattern = get_R_pattern(1)
R2_pattern = get_R_pattern(2)


include: "ASV.smk"
include: "downstream_QIIME2.smk"
include: "handling.smk"
include: "classify.smk"
include: "downstream_R.smk"
include: "import_fx.smk"
include: "QC.smk"

