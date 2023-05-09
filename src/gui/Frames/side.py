#!/usr/bin/python3
import tkinter
import customtkinter
from .. import models
from .. import env
from ..Frames import options
from .. import main

class sideframe(customtkinter.CTkFrame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        #================ configure side frame ===============

        self.rowconfigure((0, 1, 2, 3), weight=1)
        env.VARS["MAIN_LABEL"] = customtkinter.CTkLabel(self,
                                        text="GUAP toolkit",
                                        text_font=("Roboto Medium", -16)) 
        env.VARS["MAIN_LABEL"].grid(row=0, column=0,padx=5,pady=5)

        env.VARS["button_16s"] = customtkinter.CTkButton(self,
                                                text="16s rRNA", fg_color=("gray75", "gray30"),  # <- custom tuple-color
                                                command=lambda: self.button_command("16s", 
                                                main.GUAP_GUI.set_options_frame(parent,"16s")))
        env.VARS["button_16s"].grid(row=1, column=0,padx=5,pady=5,sticky="ns")


        env.VARS["button_RNA"] = customtkinter.CTkButton(self,
                                                text="RNAseq", fg_color=("gray75", "gray30"), 
                                                command=lambda: self.button_command("RNA", 
                                               main.GUAP_GUI.set_options_frame(parent,"RNA")))
        env.VARS["button_RNA"].grid(row=2, column=0,padx=5,pady=5,sticky="ns")


        env.VARS["button_WES"] = customtkinter.CTkButton(self,
                                                text="WES", fg_color=("gray75", "gray30"), 
                                                command=lambda: self.button_command("WES", 
                                               main.GUAP_GUI.set_options_frame(parent,"WES")))
        env.VARS["button_WES"].grid(row=3, column=0,padx=5,pady=5,sticky="ns")


        env.VARS["button_WGS"] = customtkinter.CTkButton(self,
                                                text="WGS", fg_color=("gray75", "gray30"), 
                                                command=lambda: self.button_command("WGS", 
                                               main.GUAP_GUI.set_options_frame(self,"WGS")))
        env.VARS["button_WGS"].grid(row=4, column=0,padx=5,pady=5,sticky="ns")

        env.VARS["change_theme_switch"] = customtkinter.CTkSwitch(self,
                                                text="Dark Mode",
                                                command=parent.change_mode)

        env.VARS["change_theme_switch"].grid(row= 8, pady=10, padx=10, sticky=tkinter.S,rowspan=4)



    def grid(self, sticky=("nswe"), **kwargs):
        super().grid(sticky=sticky, **kwargs)

    def button_command(self, button, obj):
        if button in ["WES", "RNA", "16s", "WGS"]:
            env.frames_atr["statusbar_frame"].update(f"{button} started")
            
        else:
            print(f"ERROR: {button} is not defined!")
        return obj