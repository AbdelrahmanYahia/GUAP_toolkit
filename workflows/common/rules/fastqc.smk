rule Fastqc:
    input:
        get_qc_input
    log:
        f"logs/QC/QC_{{sample}}_{RS}{{R}}.log"
    output:
        zip=f"QC/{{sample}}_R{{R}}_fastqc.zip",
        html=f"QC/{{sample}}_R{{R}}_fastqc.html"
        
    benchmark: f"benchamrks/QC/{{sample}}_{RS}{{R}}fastqc.txt"
    threads: 4
    params:
        path="QC"
    
    shell:
        "fastqc {input} --threads {threads} -o {params.path} > {log} 2>&1"
