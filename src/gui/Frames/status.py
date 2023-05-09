#!/usr/bin/python3
import customtkinter
import tkinter
from .. import env

class statusbar(customtkinter.CTkFrame):
    status_strvar = ""
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        statusbar.status_strvar = tkinter.StringVar()
        statusbar.status_strvar.set("GUAP toolkit initiated")
        env.VARS["status_var"] = customtkinter.CTkLabel(self, text_color= "Grey",
                            textvariable=statusbar.status_strvar,
                            relief=tkinter.SUNKEN)
        env.VARS["status_var"].grid(row=0, column=0, padx=1, pady=1)
    
    def update(self, txt):
        statusbar.status_strvar.set(txt)

