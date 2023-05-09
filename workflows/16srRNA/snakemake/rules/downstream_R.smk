
rule export_ps_figs:
    input:
        "ps.rds"
    output:
        directory("Phyloseq_Figures/")
    params:
        con=config["condition_name"],
        name=config["name"],
        mydir=f"{GUAP_FOLDER}/bin/16s"

    shell:
        """
        mkdir -p Phyloseq_Figures
        Rscript {params.mydir}/scripts/Export_phyloseq_figs.R \
            -i "ps.rds" \
            -n {params.name} \
            -c {params.con} \
            -w "./Phyloseq_Figures"
        """

