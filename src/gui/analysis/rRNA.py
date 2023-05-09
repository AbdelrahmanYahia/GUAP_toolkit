from tkinter import filedialog
from tkinter import *
import customtkinter
from .. import env
from .. import models




Args_INPUTS = {
    'bashdownstream': False,
    'export_figs': False,
    'condition_name': "condition",
    'name': 'GUAP_Unnamed_run',
    'skip_QC': False,
    'skip_trimmomatic': False,
    'trim_min_length': 50,
    'remove_primers': False,
    'forward_primer' : "",
    'reverse_primer' : "",
    'min_length': 50,
    'trunc_f' :0,
    'trunc_r':0,
    'trim_l':0,
    'trim_r':0,
    'maxee_f':4,
    'maxee_r':5,
    'min_overlap':10,
    'chimera_method':"consensus",
    'deblur': False,
    'deblur_trim_length': 100,
    "use_QIIME2": False,
    "choose_classifier": "qiime",
    "train": False
}


class _16s_analysis(customtkinter.CTkFrame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        env.VARS["use_QIIME2"] = True
        env.VARS["redswitch"].select()
        env.frames_atr["config_frame"].redirect_out()
        env.VARS["metadata_btn"].configure(state="normal", text="Metadata") # "normal" (standard) or "disabled" 

        for key, value in Args_INPUTS.items():
            try:
                env.INPUTS[key]
            except:
                env.INPUTS[key] = value


        # create options_frame
        self.rowconfigure((0, 1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21), weight=1)
        self.columnconfigure((0,1,2,3), weight=1)
        env.frames_atr["log_frame"].update_log("""16s analysis
In advanced options you can set:
    1. Dada2 params
    2. deblur min len
    3. trimmomatic min len
    4. naive-bayes classifier params
    5. Remove primers
    6. downstream max-depth
    7. downstream sampling depth

for more info visit:
    URL://@url@info.wiki
""")

        # labels
        env.VARS["label_QC"] = customtkinter.CTkLabel(self, text="QC",
                                        text_font=("Roboto Medium", -16))
        env.VARS["label_QC"].grid(row=0, column=0,padx=5,pady=5,sticky="nsew")

        env.VARS["label_ASV"] = customtkinter.CTkLabel(self, text="ASV",
                                        text_font=("Roboto Medium", -16))
        env.VARS["label_ASV"].grid(row=3, column=0,padx=5,pady=5,sticky="nsew")

        env.VARS["label_class"] = customtkinter.CTkLabel(self, text="Classifier",
                                        text_font=("Roboto Medium", -16))
        env.VARS["label_class"].grid(row=0, column=1,padx=5,pady=5)

        env.VARS["label_figs"] = customtkinter.CTkLabel(self, text="Export Figs",
                                        text_font=("Roboto Medium", -16))
        env.VARS["label_figs"].grid(row=3, column=1,padx=5,pady=5)

        # QC options
        env.VARS["skip_QC"] = customtkinter.CTkCheckBox(self, 
                                command= lambda: env.INPUTS["skip_QC"] == True 
                                if env.VARS["skip_QC"].get() == 1 
                                else env.INPUTS["skip_QC"] == False,
                                 text="Skip QC")

        env.VARS["skip_QC"].grid(row= 1,column=0, pady=5, padx=10,sticky="nsew")
        
        env.VARS["skip_trimmomatic"] = customtkinter.CTkCheckBox(self, 
                                        command= lambda: env.INPUTS["skip_trimmomatic"] == True 
                                        if env.VARS["skip_trimmomatic"].get() == 1 
                                        else env.INPUTS["skip_trimmomatic"] == False, 
                                        text="Skip trimming")

        env.VARS["skip_trimmomatic"].grid(row= 2,column=0, pady=5, padx=10,sticky="nsew")
        env.VARS["skip_trimmomatic"].toggle()
        
        # classifier btn
        env.VARS["classifier_btn_label"] = customtkinter.CTkLabel(self,text="Choose Classifier")
        env.VARS["classifier_btn_label"].grid(row= 1,column=1, pady=5, padx=0,sticky="nsew")
        env.INPUTS["choose_classifier"] = StringVar()
        env.INPUTS["choose_classifier"].set("qiime")
        env.VARS["classifier_btn"] = customtkinter.CTkOptionMenu(self,values=["dada2", "QIIME2 naive-bayes"],
                        command=lambda choice: env.INPUTS["choose_classifier"].set("qiime") if choice == "QIIME2 naive-bayes" 
                        else env.INPUTS["choose_classifier"].set("dada"))

        env.VARS["classifier_btn"].grid(row= 1,column=2, pady=5, padx=0,sticky="nsew")
        env.INPUTS["classifier"] = StringVar()
        env.VARS["classifier_btn"].set("QIIME2 naive-bayes")  # set initial value
        env.VARS["classifier_file_btn"] = customtkinter.CTkButton(self,text="classifier file",
                                        command= lambda:env.INPUTS["classifier"].set(
                                            filedialog.askopenfilename(initialdir= "~/", 
                                        title='Please select a directory')))

        env.VARS["classifier_file_btn"].grid(row=2, column=1, pady=5, padx=5,sticky="nsew")
        env.VARS["classifier_file_LABEL"] = customtkinter.CTkLabel(self,textvariable=env.INPUTS["classifier"])
        env.VARS["classifier_file_LABEL"].grid(row=2, column=2,sticky="nsew", pady=5)

        # ASV options
        env.VARS["use_QIIME2"] = customtkinter.CTkCheckBox(self, 
                                command= lambda: env.INPUTS["use_QIIME2"] == True 
                                if env.VARS["use_QIIME2"].get() == 1 
                                else env.INPUTS["use_QIIME2"] == False, 
                                text="Use QIIME2")
        env.VARS["use_QIIME2"].grid(row= 4,column=0, pady=5, padx=10,sticky="nsew")
        
        env.VARS["use_deblur"] = customtkinter.CTkCheckBox(self, 
                                command= lambda: env.INPUTS["deblur"] == True 
                                if env.VARS["use_deblur"].get() == 1 
                                else  env.INPUTS["deblur"] == False, 
                                text="Use DEBLUR")

        env.VARS["use_deblur"].grid(row= 5,column=0, pady=5, padx=10,sticky="nsew")

        env.VARS["adv_options_btn"] = customtkinter.CTkButton(self, 
                            fg_color=["gray85", "gray15"],   # <- no fg_color
                            text="Advanced Options",command= lambda: _16s_analysis.create_adv_16s(self))

        env.VARS["adv_options_btn"].grid(row=21, column=0, padx=15, pady=15,sticky=W+S)

        # export figs 
        env.VARS["export_figs_btn"] = customtkinter.CTkCheckBox(self, 
                                    command= lambda: env.INPUTS["export_figs"] == True 
                                    if env.VARS["export_figs_btn"].get() == 1 
                                    else env.INPUTS["export_figs"] == False, 
                                    text="Export figs")

        env.VARS["export_figs_btn"].grid(row= 4,column=1, pady=5, padx=10,sticky="nsew")
        
        env.VARS["downstream_btn"] = customtkinter.CTkCheckBox(self,
                                     command= lambda:  env.INPUTS["downstream"] == True 
                                     if env.VARS["downstream_btn"].get() == 1 
                                     else env.INPUTS["downstream"] == False, 
                                     text="QIIME2 Downstream")

        env.VARS["downstream_btn"].grid(row= 5,column=1, pady=5, padx=10,sticky="nsew")
        # run btn
        env.VARS["Run_btn_16s"]= customtkinter.CTkButton(self, 
                                command=lambda: env.print_env_input(),
                                text="RUN")
        env.VARS["Run_btn_16s"].grid(row=21,column=3, sticky=E+S, padx=15, pady=15)

        # print(f"default values:")
        # for key, value in Args_INPUTS.items():
        #     try:
        #         print(f"{key} : {value.get()}")
        #     except:
        #         print(f"{key} : {value}")


    def create_adv_16s(self):

        advanced_options_16s = customtkinter.CTkToplevel(self)
        advanced_options_16s.geometry("750x400")
        advanced_options_16s.title("16s Advanced Options")
        advanced_options_16s.grid_columnconfigure((0, 1, 2, 3,4,5,6), weight=1)
        # advanced_options_16s.grid_rowconfigure((0, 1), weight=1)


        env.VARS["ASV_label"] = customtkinter.CTkLabel(master=advanced_options_16s,text="DADA2 options", text_font=("Roboto Medium", -16))
        env.VARS["ASV_label"].grid(row=1,column=0,padx=5, pady=5,sticky="nsew",columnspan=1)

        env.VARS["dada_tff_ent"] = models.label_entry(advanced_options_16s, varname="trunc_f",labelname="trunc len forward")
        env.VARS["dada_tff_ent"].grid(row=2, column=0)
        env.VARS["dada_tfr_ent"] = models.label_entry(advanced_options_16s, varname="trunc_r",labelname="trunc len reverse")
        env.VARS["dada_tfr_ent"].grid(row=3, column=0)
        env.VARS["dada_tmf_ent"] = models.label_entry(advanced_options_16s, varname="trim_l",labelname="trim len forward")
        env.VARS["dada_tmf_ent"].grid(row=4, column=0)
        env.VARS["dada_tmr_ent"] = models.label_entry(advanced_options_16s, varname="trim_r",labelname="trim len reverse")
        env.VARS["dada_tmr_ent"].grid(row=5, column=0)
        env.VARS["dada_eff_ent"] = models.label_entry(advanced_options_16s, varname="maxee_f",labelname="maxEE forward")
        env.VARS["dada_eff_ent"].grid(row=6, column=0)
        env.VARS["dada_efr_ent"] = models.label_entry(advanced_options_16s, varname="maxee_r",labelname="maxEE reverse")
        env.VARS["dada_efr_ent"].grid(row=7, column=0)

        env.VARS["chim_meth_label"] = customtkinter.CTkLabel(master=advanced_options_16s, text="Chimera Method")
        env.VARS["chim_meth_label"].grid(row=8, column=0, padx=5, pady=5,sticky="w",columnspan=1) 
        
        env.INPUTS["chimera_method"] = StringVar()
        env.INPUTS["chimera_method"].set("consensus")
        env.VARS["chim_meth_ent"] = customtkinter.CTkOptionMenu(master=advanced_options_16s,
                                    values=["consensus", "pooled"],
                                    command=lambda choice: env.INPUTS["chimera_method"].set(choice))

        env.VARS["chim_meth_ent"].grid(row=8, column=2, padx=5, pady=5,sticky="ew")
        env.VARS["chim_meth_ent"].set("consensus")


        env.VARS["primer_label"] = customtkinter.CTkLabel(master=advanced_options_16s,text="Primers sequence", text_font=("Roboto Medium", -16))
        env.VARS["primer_label"].grid(row=1,column=4,padx=5, pady=5,sticky="ew")

        env.VARS["rm_primers"] = customtkinter.CTkCheckBox(master=advanced_options_16s, 
                                                            command= lambda: 
                                                            env.INPUTS["remove_primers"] == True 
                                                            if env.VARS["rm_primers"].get() == 1 
                                                            else env.INPUTS["remove_primers"] == False, 
                                                            text="Remove Primers", text_font=("Roboto Medium", -16))

        env.VARS["rm_primers"].grid(row= 2,column=4, pady=5, padx=10,sticky="w")
        env.VARS["ffP_ent"] = models.label_entry(advanced_options_16s, varname="forward_primer",labelname="forward primer")
        env.VARS["ffP_ent"].grid(row=3, column=4)
        env.VARS["frP_ent"] = models.label_entry(advanced_options_16s, varname="reverse_primer",labelname="reverse primer")
        env.VARS["frP_ent"].grid( row=4, column=4)
        env.VARS["min_len_ent"] = models.label_entry(advanced_options_16s, varname="min_length",labelname="minimum sequence length")
        env.VARS["min_len_ent"].grid(row=5, column=4)


        env.VARS["train_set"] = customtkinter.CTkCheckBox(master=advanced_options_16s, 
                            command= lambda: env.INPUTS["train"] == True 
                            if env.VARS["train_set"].get() == 1 
                            else env.INPUTS["train"] == False, 
                            text="Train set (classifier)", 
                            text_font=("Roboto Medium", -16))

        env.VARS["train_set"].grid(row= 6,column=4, pady=5, padx=10,sticky="w")

        env.INPUTS["t_i_seq"] = StringVar()
        env.INPUTS["t_i_seq"].set("Please select file")
        env.VARS["train_seq_file_btn"] = customtkinter.CTkButton(master=advanced_options_16s,
                                        text="sequence file (qza)",
                                        command= lambda:env.INPUTS["t_i_seq"].set(
                                            filedialog.askopenfilename(
                                            initialdir= "~/",
                                            title='Please select sequence file (qza)')))
        env.VARS["train_seq_file_btn"].grid(row=7, column=4, pady=5, padx=5,sticky="w")
        env.VARS["train_seq_file_LABEL"] = customtkinter.CTkLabel(master=advanced_options_16s,
                                            textvariable=env.INPUTS["t_i_seq"],
                                            text_color="grey50")
        env.VARS["train_seq_file_LABEL"].grid(row=7, column=5,sticky="w", pady=5)

        env.INPUTS["t_i_taxa"] = StringVar()
        env.INPUTS["t_i_taxa"].set("Please select file")
        env.VARS["train_taxa_file_btn"] = customtkinter.CTkButton(master=advanced_options_16s,
                                            text="taxa file (qza)",
                                            command= lambda:env.INPUTS["t_i_taxa"].set(
                                            filedialog.askopenfilename(
                                            initialdir= "~/", title='Please select taxa file (qza)')))

        env.VARS["train_taxa_file_btn"].grid(row=8, column=4, pady=5, padx=5,sticky="w")
        env.VARS["train_taxa_file_LABEL"] = customtkinter.CTkLabel(master=advanced_options_16s,textvariable=env.INPUTS["t_i_taxa"],text_color="grey50")
        env.VARS["train_taxa_file_LABEL"].grid(row=8, column=5,sticky="w", pady=5)

        env.VARS["trainset_minlength"] = models.label_entry(advanced_options_16s, varname="trainset_minlength",labelname="minimum sequence length")
        env.VARS["trainset_minlength"].grid(row=9,column=4)
        env.VARS["trainset_maxlength"] = models.label_entry(advanced_options_16s, varname="trainset_maxlength",labelname="maximum sequence length")
        env.VARS["trainset_maxlength"].grid(row=10,column=4)

        env.VARS['done'] = customtkinter.CTkButton(advanced_options_16s, text= "Done",
                    command=lambda: update_adv()).grid(row=10, column=0, padx=5, pady=5,sticky="nsew",columnspan=6, rowspan= 1)
        
        
        def update_adv():
            elements = [
                "trainset_minlength",
                "trainset_maxlength", "trunc_f","trunc_r",
                "trim_l","trim_r","maxee_f",
                "maxee_r","forward_primer","reverse_primer","min_length"
                ]
            for element in elements:
                key = element
                widget = env.INPUTS[key]
                try:
                    env.INPUTS[key] = widget.get()
                    print(f"{key}: {widget.get()}")
                except:
                    print(f"{key} wasn't updates, default: {env.INPUTS[key]}")
            
            advanced_options_16s.destroy()
        
        
    def grid(self, sticky=("nswe"), **kwargs):
        super().grid(sticky=sticky, **kwargs)

         # line 1
        # self.inputs['Date'] = LabelInput(
        #     recordinfo, "Date",
        #     input_class=DateEntry,
        #     input_var=tk.StringVar()
        # )
