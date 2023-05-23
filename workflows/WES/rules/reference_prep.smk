rule index_ref:
    input: "{reference}"
    output: "{reference}.fai"
    shell: "samtools faidx {input}"

rule refrence_dict:
    input: "{reference}"
    output: "{refernce.dict}"
    shell: "picard CreateSequenceDictionary -R {input}"
