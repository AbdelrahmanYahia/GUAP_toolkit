import argparse
from . import workflows
from ..utils.globals import *
from ..gui import main

# main cli class
class Cli:
    def __init__(self):
        # create argument parser
        self.parser = argparse.ArgumentParser(description="GUAPtoolkit for Genomics analysis",
                                              prog="guap", 
                                              epilog = '''This tool is still under development''',
                                              usage=f"guap [WES|RNA|16s|WGS] [--input <dir>] [--output <dir>] [...]")

        # setup workflows parsers
        subparsers = self.parser.add_subparsers(title=f"{NC_}Workflows{NC}", dest='workflow',
                                                description="choose workflow", 
                                                metavar="")
        GUI_subparser = subparsers.add_parser("GUI", help="launch GUI GUAPtoolkit")
        rRNAworkflowCli = workflows.rRNA(subparsers)
        WESworkflowCli = workflows.WES(subparsers)
        RNAworkflowCli = workflows.RNA(subparsers)
        WGSworkflowCli = workflows.WGS(subparsers)
        
        # Set up top-level arguments
        self.parser.add_argument('--full-help', action='store_true', help='print help of all workflows')

        args = self.parser.parse_args()

        if args.full_help:
            print('''this is the help for the full workflows:\nyou can use guap workflow --help to see specific ones''')
            exit(1)

        if args.workflow == None:
            self.parser.print_help()
            exit(0)

        elif args.workflow == 'GUI':
            GUI = main.GUAP_GUI()
            GUI.start()

        elif args.workflow == '16s':
            rRNAworkflowCli.run(args)

        elif args.workflow == 'WES':
            WESworkflowCli.run(args)

        elif args.workflow == 'RNA':
            RNAworkflowCli.run(args)

        elif args.workflow == 'WGS':
            WGSworkflowCli.run(args)


        else:
            print("unknown arg!")


if __name__ == '__main__':
    cli = Cli()
