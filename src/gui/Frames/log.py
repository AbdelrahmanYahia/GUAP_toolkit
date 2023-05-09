#!/usr/bin/python3
import customtkinter
from datetime import datetime
from .. import env
import tkinter

class logframe(customtkinter.CTkTextbox):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        # self.log_frame = customtkinter.CTkTextbox(parent,corner_radius=5)
        self.insert("insert", env.system_info + "\n")
        self.configure(state="disabled")  # configure textbox to be read-only

    def grid(self, sticky=("nswe"), **kwargs):
        super().grid(sticky=sticky, **kwargs)

    def update_log(self,txt):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        self.configure(state="normal")  # configure textbox to be read-only
        self.insert("end", current_time + "     " + str(txt) + "\n")
        self.configure(state="disabled")  # configure textbox to be read-only
        self.yview(tkinter.END)

