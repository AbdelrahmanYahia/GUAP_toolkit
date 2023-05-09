
rule DE:
    input:
        "Counts/{aligner}_{quantifier}.counts"
    output:
        temp("Downstream/{aligner}_{quantifier}_tmp.txt")
    params: 
        source_dir = config["GUAP_DIR"],
        clinical = config["metadata"],
        samplecol = config["sample_name"],
        controlcol = config["control_name"],
        dirstr = config["dirstr"],
        re = config["ID_extension_ptrn"],
        samplezeros = config["sample_zeros"],
        controlzeros = config["control_zeros"],
        name = config["name"],
        outdir = "Downstream"
    shell:
        """Rscript {params.source_dir}workflows/DE/GDEv1.2.R -i {input}\
         -o {params.outdir} --tab -c {params.clinical} -S {params.samplecol} \
         -C {params.controlcol} \
         -z {params.samplezeros} -x {params.controlzeros} -n {params.name}

         touch {output}
        """