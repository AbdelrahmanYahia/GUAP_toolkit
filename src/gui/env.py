#!/usr/bin/python3

# import packages 
import multiprocessing
import customtkinter
import sys
import logging
import psutil
from datetime import datetime
from datetime import date
import platform
from .analysis import parse_input
# set theme in Customtkinter
customtkinter.set_appearance_mode("Light")  # Modes: "System" (standard), "Dark", "Light"
# customtkinter.set_default_color_theme("/home/abdelrahman/Desktop/GUAP_GUI/bin/custom_theme.json")

all_threads = multiprocessing.cpu_count()
platform_ = str(platform.system())

if all_threads < 4:
    logging.error(f"\033[;31;1mYour System doesn't have enough threads to run the analyis")
    exit(1)
all_mem =int( (psutil.virtual_memory().total ) / 1000000000 )
if all_mem < 4:
    logging.error(f"\033[;31;1mYour System doesn't have enough memory to run the analyis")
    exit(1)

today = date.today()
d1 = today.strftime("%d-%b-%Y")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")

system_info = f"""{d1}    GUAP GUI started at   {current_time}

System platform:   {platform_}
total threads:      {all_threads}
total memory:     {all_mem}

"""

# stores all vars inside a dict
frames_atr = {}
VARS = {}
INPUTS = {
    'bash_continue': True,
    'snakemake': True,
    'verbose': False,
    'snakemake_dry_run': False,
    'snakemake_dag': False,
}

# stores system std.err/out
old_stderr = sys.stderr
old_stdout = sys.stdout    

def print_env_input():
    data = {}
    objects = {}
    
    for key, widget in INPUTS.items():
        try:
            data[key] = widget.get()
            print(f"updated: {key}   :   {widget.get()}")
        except:
            if str(widget).startswith(".!"):
                objects[key] = widget
            else:
                data[key] = widget
                print(f"{key}   :   {widget}")
    # yaml.safe_dump(data)
    for key, vlaue in objects.items():
        print(f"objects: {key}   :   {vlaue}")

    # parse_input.validate_sample_input(data)

