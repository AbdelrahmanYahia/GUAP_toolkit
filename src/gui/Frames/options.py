#!/usr/bin/python3
import customtkinter
from .. import env
from ..analysis import rRNA
from .. import main


class options_frame(customtkinter.CTkFrame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)        

    def grid(self, sticky=("nswe"), **kwargs):
        super().grid(sticky=sticky, **kwargs)
        
