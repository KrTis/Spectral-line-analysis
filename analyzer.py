# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 23:04:48 2016

@author: Krešimir
"""

try:
    from Tkinter import *
except ImportError:
    from tkinter import *
import ttk
import tkMessageBox
import tkFileDialog
import os
import shutil
from uncertainties import ufloat
from math import *
from operator import itemgetter, attrgetter, methodcaller
from scipy.optimize import leastsq
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from collections import OrderedDict


class App:
    size=25
    '''
    The main program
    '''
    def __init__(self,master):
        ### Dictionary containing data on all frames in the main window
        self.saving_var={'Files':[StringVar(),
                             OrderedDict([ ('Spectral lines',StringVar()), ('Spectrum',StringVar()), ('Removed',StringVar()), ('Cut',StringVar())])],
                    'Source':[StringVar(),
                              OrderedDict([('Source name',StringVar()), ('LSR',StringVar()), ('Min Frequency',StringVar()), ('Max Frequency',StringVar()), ( 'Column density of H2',StringVar()), ('Column density error',StringVar()) ])],
                    'Statistics':[StringVar(),
                                  OrderedDict([('Pearson limit',StringVar()),('Cut Temperature',StringVar()),(' Max Temperature',StringVar()),('Error cut',StringVar())])]}
        
        #### Menu
        self.menubar=Menu(master)
        self.file_menu=Menu(self.menubar,tearoff=0)
        self.file_menu.add_command(label="Hello!",command=self.open_file)
        self.file_menu.add_command(label="Save properties", command=lambda:self.save_window(self.saving_var))
        self.file_menu.add_command(label="Open properties", command=lambda:self.open_window(self.saving_var))
        self.file_menu.add_command(label="Rotational diagrams", command=self.analyze_input)
        self.menubar.add_cascade(label="File", menu=self.file_menu)
        master.config(menu=self.menubar)
        
        #### frames
        self.Frame=Frame()
        self.Frame.grid(row=0,column=0)
        self.frames={key:LabelFrame(self.Frame,text=key,width=290,height=180) for key in self.saving_var.keys()}
        j_count=0
        for j in self.saving_var.keys():
            self.frames[j].grid(row=0,column=j_count)
            self.frames[j].grid_propagate(False)
            i=0
            j_count+=1
            for key in self.saving_var[j][1].keys():
                self.input_entry(self.frames[j],key,self.saving_var[j][1][key],i,App.size)
                i=i+1
        # Progress bar at the bottom of the screen
        self.progress = ttk.Progressbar(master, orient="horizontal",length=870, style="TProgressbar", mode="determinate")
        self.progress["value"] = 0
        self.progress.grid(row=1,column=0)
    def open_file(self):
        # empty function
        print 'x'
        
    def input_entry(self,frame_name,entry_title,entry_variable,position,width,**kwargs):
        # function for creating lists of entry widgets
        Label(frame_name,text=entry_title).grid(row=position,column=0)
        Entry(frame_name,textvariable=entry_variable,width=width,**kwargs).grid(row=position,column=1)

    def save_window(self,saving_var):
        # saving window
        topl=Toplevel(master)
        topl.wm_title("Save")
        tframe=Frame(topl)
        tframe.grid(row=0,column=0)
        i=0
        for j in saving_var.keys():
            self.input_entry(tframe,j,saving_var[j][0],i,App.size)
            ttk.Button(tframe,text="Open",command=lambda tframe=tframe,saving_var=saving_var,j=j:self.file_name(tframe,saving_var,j),width=10).grid(row=i,column=2)
            i=i+1
        ttk.Button(topl,text='Save',command=lambda:self.save_entries(topl,saving_var),width=10).grid(row=1,column=0)
        
    def save_entries(self,wind,save_var):
        # saving function
        for i in save_var.keys():
            if len(save_var[i][0].get())>0:
                f=open(save_var[i][0].get(),'w')
                for j in save_var[i][1].keys():
                    f.write(str(j+':'+save_var[i][1][j].get()+'\n'))
                f.close()
        wind.destroy()
                
    def open_window(self,saving_var):
        # opening window
        topl=Toplevel(master)
        topl.wm_title("Open")
        tframe=Frame(topl)
        tframe.grid(row=0,column=0)
        i=0
        for j in saving_var.keys():
            self.input_entry(tframe,j,saving_var[j][0],i,App.size)
            ttk.Button(tframe,text="Open",command=lambda tframe=tframe,saving_var=saving_var,j=j:self.file_name(tframe,saving_var,j),width=10).grid(row=i,column=2)
            i=i+1
        i=0
        ttk.Button(topl,text='Open',command=lambda:self.open_entries(topl,saving_var),width=10).grid(row=1,column=0)
        
    def open_entries(self,wind,save_var):
        # opening function
        self.progress["value"] = 0
        for i in save_var.keys():
            if len(save_var[i][0].get())>0:
                f=open(save_var[i][0].get(),'r')
                for line in f:
                    line1=line.split(':')
                    for j in sorted(save_var[i][1].keys()):
                        if line1[0]==j:
                            save_var[i][1][j].set(line1[1])
                f.close()
        wind.destroy()
        
    def file_name(self,par,saving_var,k):
            # a file dialog for opening files
            tf=tkFileDialog.askopenfilename(parent=par,title='Choose a file')
            saving_var[k][0].set(tf)
            
    def analyze_input(self):
        # Reading the spectral lines input file
        self.progress["maximum"]=400
        self.progress["value"] = 0
        file_candidate=self.saving_var['Files'][1]['Spectral lines'].get()
        if len(file_candidate)<1:
            tkMessageBox.showerror(message="No Spectral line file chosen!")
            return
        self.input_name=self.saving_var['Files'][1]['Spectral lines'].get()
        self.input_file=open(self.input_name,'r')
        self.input_dict={0:'Area',1:'AreaE',2:'LamPos',3:'LamPosE',4:'Wid',5:'WidE',6:'Tpeak',7:'Noise',8:'LamTrans',9:'LamTransE',12:'Aij',13:'Eu',14:'gu',17:'F0'}
        xaux=[]
        k=0
        for line in self.input_file:
            line=line.strip()
            line=line.split(':')
            for j in self.input_dict.keys():
                line[j]=float(line[j])
            xaux.append(line)
        self.progress["value"]+=100
        self.molecules=list(set([y[16] for y in xaux]))
        self.cnames=list(set([y[15] for y in xaux]))
        self.species={}
        for mol in self.molecules:
            z=list(set([y[19] for y in xaux if y[16]==mol]))
            self.species[mol]=z
        self.cnames_dict={}
        self.progress["value"] +=100
        for mol in self.cnames:
            z=list(set([y[16] for y in xaux if y[15]==mol]))
            self.cnames_dict[mol]=z
        self.x={}
        self.progress["value"]+=100
        for mol in self.molecules:
            aux1={}
            for spec in self.species[mol]:
                aux={}
                for key,vals in self.input_dict.items():
                    aux[vals]=[y[key] for y in xaux]
                aux1[spec]=aux
            self.x[mol]=aux1
        self.progress["value"] +=100
        self.progress.after(800, self.terminate_progress)
        
    def terminate_progress(self):
        # resetting progressbar to 0
        self.progress["value"]=0
master = Tk()
App(master)

# display the menu
master.wm_title("Line analysis")
master.mainloop()
'''
        plt.figure()
        plt.plot([self.x['CH3OCHOv=0'][s]['Eu'] for s in self.species['CH3OCHOv=0']],[self.x['CH3OCHOv=0'][s]['Tpeak'] for s in self.species['CH3OCHOv=0']],'k.')
        plt.savefig('proba.png')
'''