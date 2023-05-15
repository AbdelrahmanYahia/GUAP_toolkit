
rule conda_env_choise:
    shell:
        """
        ## check current set flags
        echo "$-"
        ## switch conda env
        source ~/anaconda3/etc/profile.d/conda.sh && conda activate r-reticulate
        ## Confirm that set flags are same as prior to conda activate command
        echo "$-"

        ## switch conda env again
        conda activate dev
        echo "$-"
        which R
        samtools --version

        ## revert to previous: r-reticulate
        conda deactivate
        """

rule gunzip_vcf:
    input: "GATK/{sample}_{aligner}_picard.vcf.gz"
    output: "GATK/{sample}_{aligner}_picard.vcf"
    benchmark: "benchamrks/{sample}_{aligner}_GATK_gunzip.txt"

    shell: "gunzip {input}"

