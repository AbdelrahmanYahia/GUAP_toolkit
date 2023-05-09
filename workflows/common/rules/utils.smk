def get_sample_number(sample):
    return SampleTable.loc[SampleTable['sample_id'] == sample, 'sample_number'].values[0]

def get_qc_input(wildcards):
    inputs = []
    inputs.extend(expand(
        f"{config['input']}/{{sample}}_{RS}{{R}}{TAIL}.{EXT}",
        ext = ["zip", "html"],
        R = [1, 2],
        sample = samples_names
    ))
    if config["trimmomatic"] is True:
        inputs.extend(expand(
                f"trimmomatic/{{sample}}_{RS}{{R}}.{EXT}",
                ext = ["zip", "html"],
                R = [1, 2],
                sample = samples
        ))
    return inputs

def get_decompress_input(wildcards):
    inputs = []
    inputs.extend(expand(
        f"{config['input']}/{{sample}}_{RS}{{R}}{TAIL}.{EXTT}",
        ext = ["zip", "html"],
        R = [1, 2],
        sample = samples_names
    ))
    if config["trimmomatic"] is True:
        inputs.extend(expand(
                f"trimmomatic/{{sample}}_{RS}{{R}}.{EXTT}",
                ext = ["zip", "html"],
                R = [1, 2],
                sample = samples
        ))
    return inputs

def get_raw_fasta(wildcards):
    sample = wildcards.sample
    units = SampleTable.loc[sample]

    source = PATH

    mydict = dict(
        zip(
            ["R1", "R2"],
                [
                    f"{config['input']}/{sample}_{units.sample_number}{lane}_{RS}1{TAIL}.{EXT}",
                    f"{config['input']}/{sample}_{units.sample_number}{lane}_{RS}2{TAIL}.{EXT}",
                ]
        )
    )
    return mydict
