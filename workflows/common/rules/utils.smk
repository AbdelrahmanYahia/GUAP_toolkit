def get_sample_number(sample):
    return SampleTable.loc[SampleTable['sample_id'] == sample, 'sample_number'].values[0]

def get_qc_input(wildcards):
    inputs = []
    inputs.extend(expand(
        f"{PATH}/{{sample}}_{RS}{{R}}{TAIL}.{EXT}",
        R = [1, 2],
        sample = samples_names
    ))
    if config["trimmomatic"] is True:
        inputs.extend(expand(
                f"trimmomatic/{{sample}}_{RS}{{R}}.{EXT}",
                R = [1, 2],
                sample = samples
        ))
    return inputs

def get_decompress_input(wildcards):
    inputs = []
    inputs.extend(expand(
        f"{PATH}/{{sample}}_{RS}{{R}}{TAIL}.{EXTT}",
        R = [1, 2],
        sample = samples_names
    ))
    return inputs

def get_align_input(wildcards):
    sample = wildcards.sample
    units = SampleTable.loc[sample]
    if config["trimmomatic"] is True:
            mydict = dict(
        zip(
            ["R1", "R2"],
                [
                    f"trimmomatic/{sample}_{RS}1.{EXT}",
                    f"trimmomatic/{sample}_{RS}2.{EXT}",
                ],
        )
    )
    else:
        source = PATH

        mydict = dict(
            zip(
                ["R1", "R2"],
                [
                    f"{PATH}/{sample}_{units.sample_number}{lane}_{RS}1{TAIL}.{EXT}",
                    f"{PATH}/{sample}_{units.sample_number}{lane}_{RS}2{TAIL}.{EXT}",
                ],
            )
        )
    return mydict

def get_raw_fasta(wildcards):
    sample = wildcards.sample
    units = SampleTable.loc[sample]

    source = PATH

    mydict = dict(
        zip(
            ["R1", "R2"],
                [
                    f"{PATH}/{sample}_{units.sample_number}{lane}_{RS}1{TAIL}.{EXT}",
                    f"{PATH}/{sample}_{units.sample_number}{lane}_{RS}2{TAIL}.{EXT}",
                ]
        )
    )
    return mydict
