#!/usr/bin/python3
import customtkinter
import tkinter
from .. import env
from .. import models
import sys

class configframe(customtkinter.CTkFrame):

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        #============== configure options frame ===============

        # env.VARS["side_frame"].rowconfigure((0, 1, 2, 3,4,5,6,7,8), weight=1)
        # env.VARS["side_frame"].columnconfigure((0, 1, 2, 3,4,5,6,7,8), weight=1)

        env.INPUTS['threads'] = env.all_threads
        env.INPUTS['input'] = tkinter.StringVar()
        env.INPUTS['output'] = tkinter.StringVar()
        env.INPUTS['metadata'] = tkinter.StringVar()


        env.VARS["threads_str"] = tkinter.StringVar()
        env.VARS["threads_str"].set(f"Threads = {env.all_threads}")

        # redirct the analysis to log frame switch
        env.VARS["redswitch"] = customtkinter.CTkSwitch(self,
                                                text="stout/err redirect",
                                                command=self.redirect_out)
        env.VARS["redswitch"].grid(row= 0, column=2, pady=10, padx=10, sticky="s")


        # frame label
        env.VARS["label_2"] = customtkinter.CTkLabel(self,
                                        text="Basic configuration",
                                        text_font=("Roboto Medium", -16)) 
        env.VARS["label_2"].grid(row=0, column=0,padx=5,pady=5)
        # get input dir
        env.VARS["input_btn"] = customtkinter.CTkButton(self,
                                                text="In Dir",
                                                command= lambda:
                                                env.INPUTS['input'].set(tkinter.filedialog.askdirectory(initialdir= "~/", title='Please select a directory'))
                                                )

        env.VARS["input_btn"].grid(row=1, column=0,sticky="nswe", pady=5, padx=5)

        env.VARS["input_label"] = customtkinter.CTkLabel(self,
                                                  textvariable=env.INPUTS['input'])

        env.VARS["input_label"].grid(row=1, column=1, sticky="w", columnspan=2)


        # get output dir
        env.VARS["output_btn"] = customtkinter.CTkButton(self,
                                                text="Out Dir",
                                                command= lambda:
                                                env.INPUTS['output'].set(tkinter.filedialog.askdirectory(initialdir= "~/", title='Please select a directory'))
                                                )

        env.VARS["output_btn"].grid(row=2, column=0,sticky="nswe", pady=5, padx=5)
        env.VARS["output_label"] = customtkinter.CTkLabel(self,
                                                  textvariable=env.INPUTS['output'])
        env.VARS["output_label"].grid(row=2, column=1,sticky="w", columnspan=2, pady=5)


        # get metadata file 
        env.VARS["metadata_btn"] = customtkinter.CTkButton(self,
                                                text="metadata",
                                                command= lambda:
                                                env.INPUTS['metadata'].set(tkinter.filedialog.askopenfilename(initialdir= "~/", title='Please select a directory')))
        env.VARS["metadata_btn"].grid(row=3, column=0,sticky="nswe", pady=5, padx=5)
        env.VARS["metadata_btn"].configure(state="disabled", text="choose analysis to enable")

        env.VARS["metadata_LABEL"] = customtkinter.CTkLabel(self,
                                                  textvariable=env.INPUTS['metadata'])
        env.VARS["metadata_LABEL"].grid(row=3, column=1,sticky="w", columnspan=2, pady=5)


        # set number of threads
        env.VARS["N_threads"] = customtkinter.CTkSlider(self,from_=4, to=env.all_threads,number_of_steps=4-env.all_threads,
                                                command=self.get_N_threads)
        env.VARS["N_threads"].grid(row=4, column=1, columnspan=2, pady=5, padx=5, sticky="")
        env.VARS["N_threads"].set(env.all_threads)
        env.VARS["threads_label"] = customtkinter.CTkLabel(self,
                                                textvariable=env.VARS["threads_str"])              
        env.VARS["threads_label"].grid(row=4, column=0, columnspan=1, pady=5, padx=5)

        # set options 
        env.VARS["set_snakemake"] = customtkinter.CTkSwitch(self, 
                                                command= lambda:env.INPUTS['snakemake'] == True if env.VARS["set_snakemake"].get() == 1 else env.INPUTS['snakemake'] == False ,
                                                text="Use Snakemake",)
        env.VARS["set_snakemake"].grid(row= 5,column=0,columnspan=1, pady=5, padx=5, sticky="w")

        env.VARS["smk_dry_run"] = customtkinter.CTkSwitch(self,
                                                command= lambda:env.INPUTS['snakemake_dry_run'] == True if env.VARS["smk_dry_run"].get() == 1 else env.INPUTS['snakemake_dry_run'] == False ,
                                                text="Snakemake Dry Run")
        env.VARS["smk_dry_run"].grid(row= 5,column=1,columnspan=1, pady=5, padx=5, sticky="w")

        env.VARS["bash_continue"] = customtkinter.CTkSwitch(self,
                                                command= lambda:env.INPUTS['bash_continue'] == True if env.VARS["bash_continue"].get() == 1 else env.INPUTS['snakemake'] == False ,
                                                text="Continue")
        env.VARS["bash_continue"].grid(row= 5,column=2,columnspan=1, pady=5, padx=5, sticky="w")


    def grid(self, sticky=("nswe"), **kwargs):
        super().grid(sticky=sticky, **kwargs)


    def get_N_threads(self,value):
        env.VARS["threads_str"].set(f"Threads = {int(value)}")
        env.INPUTS['threads'] = value
        env.frames_atr["statusbar_frame"].update(txt=f"Number of Threads was set to {int(value)}")


    def redirect_out(self):
        if env.VARS["redswitch"].get() == 1:
            sys.stdout = models.TextRedirector(env.frames_atr["log_frame"], "stdout")
            sys.stderr = models.TextRedirector(env.frames_atr["log_frame"], "stderr")
            env.frames_atr["statusbar_frame"].update("stdout and stderr will be redirected to log screen")

        else:
            sys.stdout = env.old_stdout
            sys.stderr = env.old_stderr 
            env.frames_atr["statusbar_frame"].update("stdout and stderr will be redirected to terminal")
