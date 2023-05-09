
rule featureCounts:
    input:
        expand("{aligner}/{sample}/{bamfilename}.bam", 
                aligner = aligner,
                sample = samples, 
                bamfilename = bamfilename)

    output: 
        "Counts/{aligner}_featurecounts.precounts"

    log: 
        "logs/counts/{aligner}_featurecounts.log"

    threads: config["quant_t"]
    benchmark: "benchamrks/counts/{aligner}_fc.txt"
    params:
        gtf = ref_gtf,
        feature = config["feature"],
        extra = config['quantifier_extra_args']
        
    shell:
        """
        featureCounts -T {threads} -g {params.feature} \
            -o {output} -p \
            -a {params.gtf} \
            {input} \
            > {log} 2>&1
        """

rule modify_Fcounts:
    input:
        "Counts/{aligner}_featurecounts.precounts"
    output:
        "Counts/{aligner}_featurecounts.counts"
    params:
        GUAP_DIR = config["GUAP_DIR"],
        dirstr = dirstr,
        re = repattern
    shell:
        """tail -n +2 {input} | cut -f 1,7- > {output}
        Rscript {params.GUAP_DIR}/workflows/RNAseq/scripts/remove_pattern_in_name.R -i {output} --tab -d {params.dirstr} -r {params.re}
        """


rule htseq:
    input:
        expand("{aligner}/{sample}/{bamfilename}.bam", 
                sample = samples,
                bamfilename = bamfilename,
                aligner = aligner)

    output: 
        "Counts/{aligner}_htseq.counts"
    benchmark: "benchamrks/counts/{aligner}_htseq.txt"
    threads: config["quant_t"]

    params:
        gtf = ref_gtf,
        feature = config["feature"],
        extra = config['quantifier_extra_args']

    log: 
        "logs/counts/{aligner}_htseq.log"

    shell:
        f"""
        string='{samples}'; string=$(echo $string | tr -d "[]' "); string=$(echo $string | tr "," "\t")
        echo -e "geneid $string" > {{output}}
        htseq-count -t transcript {{params.extra}}\
                    -f bam -i {{params.feature}} \
                    -r pos -n {{threads}} {{input}} \
                    {{params.gtf}} >> {{output}} 2> {{log}}
        tail -n 5 {{output}} > Counts/htseq.stats
        head -n -5 {{output}} > Counts/temp && mv Counts/temp {{output}}
        """


        
