#!/usr/bin/python3
from datetime import datetime
from tkinter import *
import tkinter as tk
import customtkinter
from . import env

# class to redirect out
class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        super().__init__()
        self.widget = widget
        self.tag = tag
    def write(self, strg):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        self.widget.configure(state="normal")
        self.widget.insert("end", current_time +"        " + strg + "\n")
        self.widget.configure(state="disabled")
        self.widget.yview(END)

class label_entry(customtkinter.CTkFrame):
    def __init__(self, parent, varname, labelname, defualtvalue="", **kwargs):
        super().__init__(parent, corner_radius=5, **kwargs)
        # creating env var
        if varname in env.INPUTS:
            print(f"{varname} found!")
        else:
            print(f"{varname} NOT FOUND!!!")
        try:
            value = env.INPUTS[varname]
        except:
            value = defualtvalue

        env.INPUTS[varname] = StringVar()
        env.INPUTS[varname].set(value)
        print(f"varname: {varname}\n{env.INPUTS[varname]} defualts: {value}")
        # create label and entry

        self.label = customtkinter.CTkLabel(self, corner_radius=5, text= labelname)
        self.label.grid(row=0, column=0, pady=2, padx=2,sticky="nsew")
        self.entry = customtkinter.CTkEntry(self,  placeholder_text=value, textvariable=env.INPUTS[varname])
        self.entry.grid(row=0, column=1, padx=2, pady=2,sticky="ne")
        print(f"get() value for {varname} : {self.entry.get()}")

    def grid(self, row, column, rowspan= 1, columnspan = 3, sticky=("ns"), padx=5, pady=5, **kwargs):
        super().grid(sticky=sticky, row=row, column=column, rowspan= rowspan, columnspan = columnspan, padx=padx, pady=pady, **kwargs)

    def get(self,**kwargs):
        self.entry.get(**kwargs)

def optionmenu_callback(var,choice):
    env.INPUTS[var].set(choice)
    
# class Widget_input(customtkinter.CTkFrame):
#     """A widget containing a label and input together."""

#     def __init__(self, parent, label='', input_class=None,
#                  input_var=None, input_args=None, label_args=None, **kwargs):
#         super().__init__(parent, **kwargs)
#         input_args = input_args or {}
#         label_args = label_args or {}

#         self.variable = input_var

#         if input_class in (customtkinter.CTkCheckBox, customtkinter.CTkButton, customtkinter.CTkRadioButton):
#             input_args["text"] = label
#             input_args["variable"] = self.variable
#         else:
#             self.label = customtkinter.CTkLabel(self, text=label, **label_args)
#             self.label.grid(row=0, column=0)
#             input_args["textvariable"] = self.variable

#         self.input = input_class(self, **input_args)
#         self.input.grid(row=0, column=1)
#         self.columnconfigure(0, weight=1)
#         self.error = getattr(self.input, 'error', tk.StringVar())
#         self.error_label = customtkinter.CTkLabel(self, textvariable=self.error)
#         self.error_label.grid(row=2, column=0)

#     def grid(self, sticky=(tk.E + tk.W), **kwargs):
#         super().grid(sticky=sticky, **kwargs)

#     def get(self):
#         try:
#             if self.variable:
#                 return self.variable.get()
#             elif type(self.input) == tk.Text:
#                 return self.input.get('1.0', tk.END)
#             else:
#                 return self.input.get()
#         except (TypeError, tk.TclError):
#             # happens when numeric fields are empty.
#             return ''

#     def set(self, value, *args, **kwargs):
#         if type(self.variable) == tk.BooleanVar:
#                 self.variable.set(bool(value))
#         elif self.variable:
#                 self.variable.set(value, *args, **kwargs)
#         elif type(self.input) in (customtkinter.CTkCheckbutton, customtkinter.CTkRadiobutton):
#             if value:
#                 self.input.select()
#             else:
#                 self.input.deselect()
#         elif type(self.input) == customtkinter.CTkText:
#             self.input.delete('1.0', tk.END)
#             self.input.insert('1.0', value)
#         else:
#             self.input.delete(0, tk.END)
#             self.input.insert(0, value)
