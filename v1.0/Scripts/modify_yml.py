import yaml
import os
import readline

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

fname = "config.yaml"

stream = open(fname, 'r')
data = yaml.safe_load(stream)

readline.set_completer_delims(' \t\n=')
readline.parse_and_bind("tab: complete")

inclass = input("path to input classifier: ")
if os.path.isfile(inclass):
    data['indexes']["classifier"] = inclass
    with open(fname, 'w') as yaml_file:
        yaml_file.write( yaml.dump(data, default_flow_style=False))
else:
    print(f"{bcolors.FAIL}Error : No file !{bcolors.ENDC}")

