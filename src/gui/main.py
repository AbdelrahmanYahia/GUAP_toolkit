#!/usr/bin/python3
import customtkinter
import tkinter as tk
from . import env
from .Frames import side, config, status, log, options
from .analysis import rRNA

# Main App class
class GUAP_GUI(customtkinter.CTk):
    WIDTH = 1200
    HEIGHT = 700


    def __init__(self):
        super().__init__()

        # main UI elements
        self.title("GUAP")
        self.geometry(f"{GUAP_GUI.WIDTH}x{GUAP_GUI.HEIGHT}")
        self.protocol("WM_DELETE_WINDOW", self.on_closing) # call .on_closing() when app gets closed
       
       # create side frame 
        env.frames_atr["side_frame"] =  side.sideframe(self, corner_radius=5)
        env.frames_atr["side_frame"].grid(row=0, column=0, columnspan=2, rowspan=5, sticky="nsew",padx=5, pady=5)

        # create config_frame
        env.frames_atr["config_frame"] = config.configframe(self,corner_radius=5)
        env.frames_atr["config_frame"].grid(row=0, column=2, columnspan=6, rowspan=5, sticky="nsew",padx=5, pady=5)

        # create options_frame
        env.frames_atr["options_frame"] = options.options_frame(self,corner_radius=5)
        env.frames_atr["options_frame"].grid(row=5, column=0, columnspan=8, rowspan=6, sticky="nsew",padx=5, pady=5)

        # create Statusbar frame
        env.frames_atr["statusbar_frame"] = status.statusbar(self,corner_radius=5)
        env.frames_atr["statusbar_frame"].grid(row=11, column=0, columnspan=20, rowspan=1, sticky= tk.W+tk.E,padx=5, pady=5)

        # create log_frame
        env.frames_atr["log_frame"] = log.logframe(self,corner_radius=5)
        env.frames_atr["log_frame"].grid(row=0,column=10, columnspan=8, rowspan=11, sticky="nswe",padx=5, pady=5)
        
        env.VARS["change_theme_switch"].toggle()
        env.VARS["set_snakemake"].toggle()
        env.VARS["bash_continue"].toggle()


        #================ configure status bar ===============

        # grid configure
        self.grid_columnconfigure(15,weight=1)
        self.grid_rowconfigure(8,weight=1)

    def on_closing(self, event=0):
        self.destroy()

    def start(self):
        self.mainloop()

    def change_mode(self):
        if env.VARS["change_theme_switch"].get() == 1:
            customtkinter.set_appearance_mode("dark")
        else:
            customtkinter.set_appearance_mode("light")

    def set_options_frame(self, txt):
        if txt not in ["16s", "WGS", "RNA", "WES"]:
            print("ERROR!!!")
            exit(1)

        elif txt == "16s":
            env.frames_atr["options_frame"] = rRNA._16s_analysis(self,corner_radius=5)
            env.frames_atr["options_frame"].grid(row=5, column=0, columnspan=8, rowspan=6, sticky="nsew",padx=5, pady=5)





