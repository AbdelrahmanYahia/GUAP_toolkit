#!/bin/env python

import datetime 
now = datetime.datetime.now()
start_time = now.strftime("%Y-%m-%d %H:%M:%S")
import os
os.environ['start_time'] = start_time

from src.cli import cli
from src.utils.globals import *

if __name__ == '__main__':
    print(f"""
 ________________________________________________
|                                                 |
|                                                 |
|                                                 |
|    ▄████     █    ██     ▄▄▄          ██▓███    |
|   ██▒ ▀█▒    ██  ▓██▒   ▒████▄       ▓██░  ██▒  |
|  ▒██░▄▄▄░   ▓██  ▒██░   ▒██  ▀█▄     ▓██░ ██▓▒  |
|  ░▓█  ██▓   ▓▓█  ░██░   ░██▄▄▄▄██    ▒██▄█▓▒ ▒  |
|  ░▒▓███▀▒   ▒▒█████▓     ▓█   ▓██▒   ▒██▒ ░  ░  |
|   ░▒   ▒    ░▒▓▒ ▒ ▒     ▒▒   ▓▒█░   ▒▓▒░ ░  ░  |
|    ░   ░    ░░▒░ ░ ░      ▒   ▒▒ ░   ░▒ ░       |
|  ░ ░   ░     ░░░ ░ ░      ░   ▒      ░░         |
|        ░       ░              ░  ░              |
|                                                 |
|                                                 |
|       {GRE}GUAP toolkit for Genomics analysis{NC}        |
| ________________________________________________|
    """)
    cli = cli.Cli()