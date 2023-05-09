import os
import sys
import argparse
from ..utils.globals import *
from ..utils.parse_input import parse_input_args


def smk_cmd(all_args):
    smk_cmd = ""

    if all_args['dry_run'] is True:
        smk_cmd += " -n"
        glogger.prnt_warning("performing dry run")
        if all_args['export_dag'] is True:
            smk_cmd += f" --rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'"
            glogger.prnt_warning("exporting dag")
        else:
            pass
    elif all_args['export_dag'] is True and all_args['dry_run'] != True:
        smk_cmd += f"--rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'"

    else:
        smk_cmd = smk_cmd + f"{all_args['smk_extra_args']}"


    snakemake_cmd = f"snakemake --snakefile '{GUAP_DIR}/workflows/RNAseq/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} {smk_cmd}"
    return snakemake_cmd

# Create a custom action to preserve quoted strings
class StoreQuotedString(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If the value is surrounded by quotes, remove them
        if isinstance(values, str) and values.startswith('"') and values.endswith('"'):
            values = values[1:-1]
        setattr(namespace, self.dest, values)

# base class for subcommands (workflows)
class WorkflowCli:
    def __init__(self,subparser):
        self.parser = subparser.add_parser(self.name, help=self.help, usage=self.usage)
        self.add_arguments(self.parser)

    def add_arguments(self, parser):
        pass
    
    def run(self, args):
        pass
