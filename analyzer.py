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
#from shutil import rmtree
from uncertainties import ufloat
from math import *
from operator import itemgetter, attrgetter, methodcaller
from matplotlib.colors import LogNorm
import numpy as np
from matplotlib.pyplot import figure, plot,show, contour, hist2d,hist,colorbar,savefig,errorbar,close,ioff,ion
from collections import OrderedDict
from scipy import constants
from json import dump as jdump
from json import load as jload
from astroML.plotting.mcmc import convert_to_stdev
from sklearn import linear_model
from sklearn.gaussian_process import GaussianProcess
from scipy.optimize import curve_fit
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
class VerticalScrolledFrame(Frame):
    """A pure Tkinter scrollable frame that actually works!
    * Use the 'interior' attribute to place widgets inside the scrollable frame
    * Construct and pack/place/grid normally
    * This frame only allows vertical scrolling

    """
    def __init__(self, parent, *args, **kw):
        Frame.__init__(self, parent, *args, **kw)            

        # create a canvas object and a vertical scrollbar for scrolling it
        vscrollbar = Scrollbar(self, orient=VERTICAL)
        vscrollbar.pack(fill=Y, side=RIGHT, expand=FALSE)
        canvas = Canvas(self, bd=0, highlightthickness=0,
                        yscrollcommand=vscrollbar.set)
        canvas.config(width=666, height=666)
        canvas.pack(side=LEFT, fill=BOTH, expand=TRUE)
        vscrollbar.config(command=canvas.yview)

        # reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # create a frame inside the canvas which will be scrolled with it
        self.interior = interior = Frame(canvas)
        interior_id = canvas.create_window(0, 0, window=interior,
                                           anchor=NW)

        # track changes to the canvas and frame width and sync them,
        # also updating the scrollbar
        def _configure_interior(event):
            # update the scrollbars to match the size of the inner frame
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the canvas's width to fit the inner frame
                canvas.config(width=interior.winfo_reqwidth())
        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the inner frame's width to fill the canvas
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())
        canvas.bind('<Configure>', _configure_canvas)
        
class plotting(object):
    @classmethod
    def simple_plot(self,plo1,mol,species,xx,A,B,opts,opts1,opts2,opts3,scale,args):
        ioff()
        if mol!='All':
            x,y=comp.combine_species(mol, species[mol],xx,A,B)
        else:
            x=[]
            y=[]
            for m in species.keys():
                x1,y1=comp.combine_species(m, species[m],xx,A,B)
                x+=x1
                y+=y1
                
        with open('rotdiag.json','r') as outfile:
           rot=jload(outfile,encoding='utf-8')

        fig=figure()
        ax=fig.add_subplot(111)
        if opts3[0]=='Y':
            ax.plot(x,y,args)
        if opts[0]=='Y':
            counts, xedges, yedges, im=ax.hist2d(x,y,norm=LogNorm(),bins=scale)
            colorbar(im)
        if opts1[0]=='Y':
            H, xbins, ybins = np.histogram2d(x,y,bins=scale)
            Nsigma = convert_to_stdev(np.log(H))
            ax.contour(0.5 * (xbins[1:] + xbins[:-1]),0.5 * (ybins[1:] + ybins[:-1]),Nsigma.T,levels=[0.6827,0.6827,0.9545, 0.9545], colors=['.25','.25','0.5','0.5'],linewidths=2)
        
        if A=='Eu' and B=='Nu/gu':
            y=[comp.rot(x1,rot[mol]['Rot'][0],rot[mol]['Rot'][1]) for x1 in x]
            plot(x,y)
        elif opts2=='Yes':
            x=np.array(x).reshape((len(x),1))
            y=np.array(y).reshape((len(y),1))
            regr = linear_model.LinearRegression()
            regr.fit(x,y)        
            plot(x, regr.predict(x), color='blue',linewidth=3)
        ax.set_xlabel(App.tex_output_dict[A])
        ax.set_ylabel(App.tex_output_dict[B])
        pframe=Frame(plo1)
        pframe.grid(row=2,column=0)
        canvas = FigureCanvasTkAgg(fig, master=pframe)
        #canvas.show()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2TkAgg(canvas, pframe)
        toolbar.update()
        canvas._tkcanvas.pack()
        close(fig)
        ion()
        def on_key_event(event):
            key_press_handler(event, canvas, toolbar)
        
        canvas.mpl_connect('key_press_event', on_key_event)

        close(fig)
        #ion()
    @classmethod
    def saving_plot(self,pclass,Eu, Nu, NuE,ymol,pmol,ax,label,ypos,**kwargs):
            ax.set_xlim(0,450.)
            lims=ax.get_ylim()
            if lims[0]>.9*min(Nu) or lims[0]<15.:
                 ax.set_ylim(bottom=.9*min(Nu))
            if lims[1]<1.1*max(Nu):
                ax.set_ylim(top=1.1*max(Nu))
            ax.set_title("$ "+pmol+"$")
            ax.set_xlabel(App.tex_output_dict['Eu'])
            ax.set_ylabel(App.tex_output_dict['Nu/gu'])
            
            
            if ymol['Rot'][0]>0 and ymol['Rot'][0]<500. :
                if pclass.mol in pclass.z_file.keys():
                    ymol['Z'],ymol['ZStdev']=comp.get_z(ymol['Rot'][0],ymol['RotSigmas'][0],pclass.z_file[pclass.mol]['Fit'],pclass.z_file[pclass.mol]['Stdev'])
                    ymol['N'], ymol['Nerr']=comp.column_density(ymol['Z'],ymol['ZStdev'],ymol['Rot'][1],ymol['RotSigmas'][1])
            
                    if  ymol['N']!=float("inf") and ymol['Nerr']!=float("inf"):
                        text_plot=label+" ${:.2uL}$ K".format(ufloat(ymol['Rot'][0],ymol['RotSigmas'][0]))
    
                        ax.errorbar(Eu,Nu,yerr=NuE,**kwargs)
                        x_line=np.linspace(0, 450, 100)
                        y_line=[comp.rot(x1,ymol['Rot'][0],ymol['Rot'][1]) for x1 in x_line]
                        ax.plot(x_line,y_line,label=label)
                        
                        text_plot+="\n${:.2uL}$".format(ufloat(ymol['N'],ymol['Nerr']))+r"$\,\mathrm{cm^{-2}}$"
    
                        ax.text(.65,ypos,text_plot,transform = ax.transAxes)
                        return True
                    else:
                        ax.errorbar(Eu,Nu,yerr=NuE,fmt='k.',ecolor='black')
                else:
                    ax.errorbar(Eu,Nu,yerr=NuE,fmt='k.',ecolor='black')
            else:
                ax.errorbar(Eu,Nu,yerr=NuE,fmt='k.',ecolor='black')
            

class comp:
    @classmethod
    def len1(self,x):
        try:
            len(x)
            return len(x)
        except TypeError:
            return 1
    @classmethod
    def cutting_list(self,limit,data):
        length=len(data)
        remove=[]
        for i in range(length):
            if data[i]>limit:
                remove.append(i)
        return remove
    @classmethod
    def finding_in_list(self,data,what):
        length=len(data)
        remove=[]
        for i in range(length):
            if data[i]==what:
                remove.append(i)
        return remove
    @classmethod
    def delete_list(self,input_list,remove_list):
        return [input_list[i] for i in range(len(input_list)) if i not in remove_list]
    @classmethod
    def delete_dict_list(self, input_dict,remove_list):
        return {key:self.delete_list(input_dict[key],remove_list) for key in input_dict.keys()}
    @classmethod
    def doppler_velocity(self,measured,errmeas,database,errdb):
        return [[0.001*constants.c*(database[i]-float(measured[i]))/database[i] for i in range(len(measured))],
                 [sqrt((errdb[i]/database[i])**2+(errmeas[i]/float(measured[i]))**2)*0.001*constants.c*(database[i]-float(measured[i]))/database[i] for i in range(len(measured))]]
    @classmethod
    def save_to_json(self, saving_var, progress):
        # Reading the spectral lines input file
        progress["value"]=0
        progress["maximum"]=500
        file_candidate=saving_var['Files'][1]['Spectral lines'].get()
        if len(file_candidate)<1:
            tkMessageBox.showerror(message="No Spectral line file chosen!")
            return
        fc_testing=file_candidate.split('.')
        if fc_testing[1]=='json':
            return
        input_name=saving_var['Files'][1]['Spectral lines'].get()
        input_file=open(input_name,'r')
        progress["value"] +=100
        
        xaux=[]
        for line in input_file:
            line=line.strip()
            line=line.split(':')
            for j in App.input_dict.keys():
                if j!=2:                
                    line[j]=float(line[j])
            xaux.append(line)
        
        progress["value"] +=100
        
        molecules=list(set([y[16] for y in xaux]))
        cnames=list(set([y[15] for y in xaux]))
        species={}
        for mol in molecules:
            species[mol]=list(set([y[19] for y in xaux if y[16]==mol]))   
        cnames_dict={}
        progress["value"] +=100
        
        for mol in cnames:
            cnames_dict[mol]=list(set([y[16] for y in xaux if y[15]==mol]))
        progress["value"]+=100
        x={}
        
        for mol in molecules:
            aux1={}
            for spec in species[mol]:
                aux={}
                for key,vals in App.input_dict.items():
                    aux[vals]=[y[key]  for y in xaux if (y[16]==mol and y[19]==spec) ]
                aux1[spec]=aux
                
            x[mol]=aux1
        for mol in molecules:
            for spec in species[mol]:
                x[mol][spec]['NuPos']=[float(y) for y in x[mol][spec]['NuPosS']]
                x[mol][spec]['SpecPlot']=comp.quantum_numbers(spec,mol,cnames_dict)
        fc=file_candidate.split('.')
        fc=fc[0]+'.json'
        with open(fc, 'w') as outfile:
		jdump(x, outfile, encoding='utf-8')
        saving_var['Files'][1]['Spectral lines'].set(fc)
        
        with open('cnames.json','w') as outfile:
            jdump(cnames_dict,outfile,encoding='utf-8')
        progress["value"] +=100
        return x, molecules, species, cnames_dict
    @classmethod
    def evodd(self,n):
	if(n%2==0):
		return 'e'
	else:
		return 'o'
    
    @classmethod
    def W(self,lin,i,mol,spec):
        if(lin['gu'][i]<1.):
            yy=log(lin['Area'][i]*lin['NuPos'][i])-lin['Aij'][i]*log(10.)+6.36762145-log(0.75)  
            y,upstate=comp.gcorr(mol,spec,yy)		
        else:
		yy=log(lin['Area'][i]*lin['NuPos'][i]/lin['gu'][i])-lin['Aij'][i]*log(10.)+6.36762145-log(0.75) 
        ye=np.sqrt((lin['AreaE'][i]/lin['Area'][i])**2+(2*lin['NuPosE'][i]/lin['NuPos'][i])**2)
        return yy,ye
    @classmethod
    def Nu_gu(self,lin,mol,spec):
        y=[]
        ye=[]
        for i in range(len(lin['Area'])):
            y1,ye1=comp.W(lin,i,mol,spec)
            y.append(y1)
            ye.append(ye1)
        return y,ye
    @classmethod
    def gcorr(self,chemn,qns,yy):
	upstate=1
    	if(chemn.strip()=='CH3OCH3'or chemn.strip()=='(CH3)2COv=0'):
		qns1=qns.replace('(-','(*').replace(',-',',*').replace('+-','+*').replace('-+','*+').replace('--','**')
		qns2=qns1.replace('(',',').replace(')',',').replace('-','').replace(',v=','').replace(',F=','').split(',')
		qns3=[m.replace('*','-') for m in qns2]
		if(len(qns3)>6):	
			nu=int(qns3[0])
			kau=int(qns3[1])
			kcu=int(qns3[2])
			nd=int(qns3[3])
			kad=int(qns3[4])
			kcd=int(qns3[5])
			extra1=qns3[6]
			if((comp.evodd(kau)=='e' and comp.evodd(kcu)=='e')or(comp.evodd(kau)=='o' and comp.evodd(kcu)=='o')):
				if(extra1=='AA'):
					korekcijaK=6
				if(extra1=='AE'):
					korekcijaK=2
				if(extra1=='EA'):
					korekcijaK=4
				if(extra1=='EE'):
					korekcijaK=16
				
			if((comp.evodd(kau)=='e' and comp.evodd(kcu)=='o')or(comp.evodd(kau)=='o' and comp.evodd(kcu)=='e')):
				if(extra1=='AA'):
					korekcijaK=10
				if(extra1=='AE'):
					korekcijaK=6
				if(extra1=='EA'):
					korekcijaK=4
				if(extra1=='EE'):
					korekcijaK=16
			upstate=korekcijaK*(2*nu+1)*.5
    	if(chemn.strip() in ['13CH3OHvt=0','CH3OHvt=0','CH3OCHOv=0','CH3OCHOv=1']):
		qns1=qns.replace('(-','(*').replace(',-',',*').replace('+-','+*').replace('-+','*+').replace('--','**')
		qns2=qns1.replace('(',',').replace(')',',').replace('-','').replace(',v=','').replace(',F=','').split(',')
		qns3=[m.replace('*','-') for m in qns2]
		Ju=int(qns3[0])
		upstate=(2*Ju+1)
	if(chemn.strip() in ['CH3CNv=0','CH3CCHv=0']):
		qns1=qns.replace('(-','(*').replace(',-',',*').replace('+-','+*').replace('-+','*+').replace('--','**')
		qns2=qns1.replace('(',',').replace(')',',').replace('-','').replace(',v=','').replace(',F=','').split(',')
		qns3=[m.replace('*','-') for m in qns2]
		if(len(qns3)>3):
			ku=int(qns3[0])
			upstate=(2*ku+1)
			if(ku>0):
				upstate=2*upstate
	return yy-log(upstate) , upstate        
    @classmethod
    def read_json(self,file_name):
        with open(file_name,'r') as infile:
            x=jload(infile,encoding='utf-8')
        molecules=x.keys()
        species={}
        for mol in molecules:
            species[mol]=x[mol].keys()
        k=0
        for mol in molecules:
            for spec in species[mol]:
                    column_names=x[mol][spec].keys()
        return x, molecules, species,column_names
    @classmethod
    def combine_species(self,mol, species,x,A,B):
        xx=[]
        xy=[]
        for spec in species:
            xx+=x[mol][spec][A]
            xy+=x[mol][spec][B]
        return xx, xy
    @classmethod
    def z_fit(self,T, a, b):
            ### A fit for partition functions
            return a * np.exp(b*T)*(T**1.5)
    @classmethod
    def get_z(self,T,M_T,pars,M_pars):
        ### A function that returns the partition function at temperature T
        if comp.len1(pars)==1:
            return pars*T**1.5, M_pars*np.sqrt(T)*M_T
        if comp.len1(pars)==2:
            fz=comp.z_fit(T,pars[0],pars[1])
            return fz, fz*np.sqrt((M_pars[0]/pars[0])**2+(T*M_pars[1])**2+((1.5/T+pars[1])*M_T)**2)
    @classmethod
    def fit_population(self,Eu,Nu,NuE):
        Eu=np.asarray(Eu)
        Nu=np.asarray(Nu)
        NuE=np.asarray(NuE)
        params,covar=curve_fit(comp.rot, Eu,Nu,[100, 20.],sigma=NuE)
        sigmas=np.sqrt(np.diag(covar))
        return params.tolist(), covar.tolist(), sigmas.tolist()
    @classmethod
    def fit_rotational(self, species, x):
        Eu=[]
        Nu=[]
        NuE=[]
        for spec in species:
                Eu+=x[spec]['Eu']
                Nu+=x[spec]['Nu/gu']
                NuE+=x[spec]['Nu/guE']
        a, b, c=comp.fit_population(Eu, Nu, NuE)
        return Eu, Nu, NuE, a, b,c
    @classmethod        
    def rot(self,Eu,a,b):
            return -Eu/a+b
    @classmethod
    def column_density(self,Z, Zerr, NZ, NZerr):
        N=Z*np.exp(NZ)
        return N, N*np.sqrt((Zerr/Z)**2+NZerr**2)

    @classmethod
    def quantum_numbers(self, qns,ch,cnames_dict):
        empty=[0,0,'','','']
        extra=''
        for mol in cnames_dict.keys():
            print ch, cnames_dict[mol]
            if ch in cnames_dict[mol]:
                qns=qns.replace('(-','(*').replace(',-',',*').replace('+-','+*').replace('-+','*+').replace('--','**')
                if comp.len1(qns)>0 and not 'N/A' in qns and mol in App.molecular_specifier.keys():
                    
                    if App.molecular_specifier[mol] in [303,1404,14041]:
                        exiting_string='$K_a^u-K_c^u$'
                        
                        qns=qns.replace('(',',').replace(')',',').replace('-','').replace(',v=','').replace(',F=','').split(',')
                        qns=[m.replace('*','-') for m in qns] 
                        if comp.len1(qns)<7:
                            return empty
                        if App.molecular_specifier[mol] in [1404,14041]:
                            extra=qns[6]                            
                        return [100*int(qns[0])+int(qns[3]),100*int(qns[1])+int(qns[4]),str(comp.evodd(int(qns[1]))+comp.evodd(int(qns[2]))+'-'+comp.evodd(int(qns[4]))+comp.evodd(int(qns[5]))),extra,exiting_string]
                        
                    if App.molecular_specifier[mol] in [1202,101]:
                        if App.molecular_specifier[mol]==1202:
                            exiting_string='$J-J$'
                        else:
                            exiting_string='$N-N$'
                            
                        qns=qns.split('-')
                        qns=[m.replace('*','-') for m in qns]
                        if comp.len1(qns[0])>3 or comp.len1(qns[1])>3:
                            return empty
                        return [100*int(qns[0])+int(qns[1]),100*int(qns[0])+int(qns[1]),str(qns[0]+'-'+qns[1]),extra,exiting_string]
                        
                    if App.molecular_specifier[mol] in [202,203]: 
                       exiting_string='$K-K$'
                       
                       qns=qns.replace('(',',').replace(')',',').replace(',F=',"").replace(',-',',').replace(',,',',').split(',')
                       qns=[m.replace('*','-') for m in qns]
                       return [100*int(qns[0])+int(qns[2]),100*int(qns[1])+int(qns[3]),str(qns[1]+'-'+qns[3]),extra,exiting_string]
                       
                    if App.molecular_specifier[mol] in [1021,102]:
                       if App.molecular_specifier[mol]==102:
                           exiting_string='$J-J$'
                       else:
                           exiting_string='$F-F$'
                           
                       qns=qns.replace('-',',').replace('N=','').replace('J=','').replace(',F=',',').replace(',,',',').split(',')
                       qns=[m.replace('*','-') for m in qns] 
                       return [100*int(qns[0])+int(qns[1]),100*int(qns[2])+int(qns[3]),str(qns[2]+'-'+qns[3]),extra,exiting_string]
                    
                    if App.molecular_specifier[mol] in [234]:
                        exiting_string='$F-F$'
                        
                        qns=qns.replace('N=','').replace('J=','').replace('F=','').replace('/2-',',').replace('/2,',',').replace('/2',',').replace('-',',').split(',')
                        qns=[m.replace('*','-') for m in qns]
                        return [100*int(qns[0])+int(qns[1]),100*int(qns[2])+int(qns[3]),str(qns[2]+'-'+qns[3]),extra,exiting_string]
                        
        return empty
        
class App:
    '''
    The main program
    '''
    size=25
    input_dict={0:'Area',1:'AreaE',2:'NuPosS',3:'NuPosE',4:'Wid',5:'WidE',6:'Tpeak',7:'Noise',8:'NuTrans',9:'NuTransE',12:'Aij',13:'Eu',14:'gu',17:'F0'}
    tex_output_dict={'Area':'$A$','AreaE':'$M_A$','NuPos':'$\\nu\,[\mathrm{MHz}]$','NuPosE':'$M_{\\nu}\,[\mathrm{MHz}]$','Wid':'$W$','WidE':'$M_W$','Tpeak':'$T_{\mathrm{peak}}\,[\mathrm{K}]$',
    'Noise':'$\sigma_{\mathrm{RMS}}$','NuTrans':'$\\nu_{ul}\,[\mathrm{MHz}]$','NuTransE':'$M_{\\nu_{ul}}\,[\mathrm{MHz}]$','Aij':'$A_{ul}$','Eu':'$E_u\,[\mathrm{K}]$','gu':'$g_u$',
    'F0':'$\\nu_0\,\mathrm{MHz}$','Nu/gu':'$\ln(N_u/g_u)$','Nu/guE':'$M_{\ln(N_u/g_u)}$','Relative error A':'$M_A/A$','Relative error W':'$M_W/W$','Doppler V':'$v\,[\mathrm{km/s}]$','Doppler VE':'$M_v\,[\mathrm{km/s}]$','Sigma':'$\sigma$'}
    changeable_list=[ u'Noise', u'Relative error W', u'Area', u'Tpeak',  u'Wid', u'Eu', u'Doppler V', u'Sigma', u'Relative error A']    
    files_list=['Spectral lines','Cleaned Spectral lines','Spectrum','Removed','Special Eu','Plotting names','Partition functions']
    source_list=['Source name','LSR','Min Frequency','Max Frequency','Column density of H2','Column density error'] 
    statistics_list=['Pearson limit','Max Eu','Error cut']
    molecular_specifier={"2-Propynal": 303,
            "Acetaldehyde": 1404,
            "Acetone": 1404,
            "Butatrienylidene":  303,
            "Carbon Monosulfide":  1202,
            "Carbonyl Sulfide":1202,
            "Cyanamide": 14041,
            "Cyanoallene": 303,
            "Cyanoformaldehyde": 303,
            "Cyclopropenylidene": 303,
            "Dimethyl ether": 1404,
            "Ethyl Cyanide": 303,
            "Ethylene Oxide": 303,
            "Formamide":303,
            "gauche-Ethanol":1404,
            "Glycolaldehyde":303,
            "Isocyanic Acid":303,
            "Ketene": 303,
            "Ketenimine":303,
            "Methanimine":303,
            "Methanol": 1404,
            "Methylamine":0,
            "Methyl Acetylene": 202,
            "Methyl diacetylene":202,
            "Methyl Cyanide": 202,
            "Methyl Formate": 1404,
            "Potassium cyanide, potassium isocyanide":303,
            "Propenal":303,
            "Protonate 2-proynenitrile": 1021,
            "Silicon Cyanide": 234,
            "Silicon Monoxide": 1202,
            "Sulfur Dioxide": 303,
            "Thioformaldehyde": 303,
            "Thioformylium": 101,
            "Thioxoethenylidene":102,
            "Vinyl Cyanide":303,
            "Water": 303}


    def __init__(self,master):
        ### Dictionary containing data on all frames in the main window
        self.saving_var={'Files':[StringVar(),
                             OrderedDict([ (x,StringVar()) for x in App.files_list])],
                    'Source':[StringVar(),
                              OrderedDict([(x,StringVar()) for x in App.source_list])],
                    'Statistics':[StringVar(),
                                  OrderedDict([(x,StringVar()) for x in App.statistics_list])]}
        
        #### Menu
        self.menubar=Menu(master)
        self.file_menu=Menu(self.menubar,tearoff=0)
        self.analysis_menu=Menu(self.menubar,tearoff=0)
        self.file_menu.add_command(label="Hello!",command=self.open_file)
        self.file_menu.add_command(label="Save properties", command=lambda:self.save_window(self.saving_var))
        self.file_menu.add_command(label="Open properties", command=lambda:self.open_window(self.saving_var))
        self.analysis_menu.add_command(label="Save to JSON", command=lambda:comp.save_to_json(self.saving_var,self.progress))
        self.analysis_menu.add_command(label="Save Partition functions to JSON", command=self.reading_partition_functions)
        self.analysis_menu.add_command(label="Histograms", command=self.hist_cut)
        self.analysis_menu.add_command(label="Apply criteria", command=self.clean)
        self.analysis_menu.add_command(label="Simple Plot", command=self.plotting_window)
        self.menubar.add_cascade(label="File", menu=self.file_menu)
        self.menubar.add_cascade(label="Analysis", menu=self.analysis_menu)
        master.config(menu=self.menubar)
        
        #### frames
        self.Frame=Frame()
        self.Frame.grid(row=0,column=0)
        self.frames={key:ttk.LabelFrame(self.Frame,text=key,width=290,height=180) for key in self.saving_var.keys()}
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
                            save_var[i][1][j].set(line1[1].strip())
                f.close()
        wind.destroy()
        
    def file_name(self,par,saving_var,k):
            # a file dialog for opening files
            tf=tkFileDialog.askopenfilename(parent=par,title='Choose a file')
            saving_var[k][0].set(tf)
            
    def check_cut_files(self,text,else_set):
        try:
                f=float(text)
        except ValueError:
                f=else_set
        return f
            
    

	
    def clean(self):
        ioff()
        ### cleans spectra by criteria from the Statistics column
        self.progress["value"]=0
        self.progress["maximum"]=1000
        
        if len(self.saving_var['Files'][1]['Removed'].get())>0:
            removed_file=open(self.saving_var['Files'][1]['Removed'].get(),'r')
            removed_list=[line.strip() for line in removed_file]
        else:
            removed_list=[]
            
        self.ErrorCut=self.check_cut_files(self.saving_var['Statistics'][1]['Error cut'].get(),100.)
        self.MaxEu=self.check_cut_files(self.saving_var['Statistics'][1]['Max Eu'].get(),1000.)
        self.special_Eu={}
        file_candidate=self.saving_var['Files'][1]['Spectral lines'].get()
        self.progress["value"] +=100
        
        if len(file_candidate)<1:
            tkMessageBox.showerror(message="No Spectral line file chosen!")
            return
        if ".json" not in file_candidate:
            x,molecules,species,cnames_dict=comp.save_to_json(self.saving_var,self.progress)
            progress["value"]=100
        else:
            with open(file_candidate, 'r') as outfile:
                x=jload( outfile, encoding='utf-8')
            with open('cnames.json','r') as outfile:
                cnames_dict=jload(outfile, encoding='utf-8')
        self.progress["value"] +=100
        
        partfun=self.saving_var['Files'][1]['Partition functions'].get()
        if len(partfun)>0:
            partfun=partfun.split('.')
            if partfun[1]!='json':
                 self.reading_partition_functions
            with open(self.saving_var['Files'][1]['Partition functions'].get()) as outfile:
                    self.z_file=jload(outfile,encoding='utf-8')
        else:
            tkMessageBox.showerror(message="Partition functions not given!")
            return
        self.progress["value"] +=100
        
        molecules=x.keys()
        molecules=[mol for mol in molecules if mol not in removed_list]
        x={mol:x[mol] for mol in x.keys() if mol in molecules}
        species={mol:x[mol].keys() for mol in molecules}
        self.progress["value"] +=100
        
        for mol in molecules:
            for spec in species[mol]:
                x[mol][spec]['Relative error A']=[x[mol][spec]['AreaE'][i]/x[mol][spec]['Area'][i] for i in range(len(x[mol][spec]['Area']))]
                x[mol][spec]['Relative error W']=[x[mol][spec]['WidE'][i]/x[mol][spec]['Wid'][i] for i in range(len(x[mol][spec]['Wid']))]
                x[mol][spec]['Sigma']=[x[mol][spec]['Tpeak'][i]/x[mol][spec]['Noise'][i] for i in range(len(x[mol][spec]['Wid']))]
                x[mol][spec]['Doppler V'],x[mol][spec]['Doppler VE']=comp.doppler_velocity(x[mol][spec]['NuPos'],x[mol][spec]['NuPosE'],x[mol][spec]['NuTrans'],x[mol][spec]['NuTransE'])
                x[mol][spec]['Nu/gu'], x[mol][spec]['Nu/guE']=comp.Nu_gu(x[mol][spec],mol,spec)      
                if mol not in self.special_Eu.keys():
                    remove=comp.cutting_list(self.MaxEu,x[mol][spec]['Eu'])
                else:
                    remove=comp.cutting_list(self.special_Eu[mol],x[mol][spec]['Eu'])
                remove+=comp.cutting_list(self.ErrorCut,x[mol][spec]['Relative error A'])
                remove+=comp.cutting_list(self.ErrorCut,x[mol][spec]['Relative error W'])
                remove+=comp.cutting_list(-3.,[(-1.)*z for z in x[mol][spec]['Sigma']])
                x[mol][spec]=comp.delete_dict_list(x[mol][spec],remove)
        self.progress["value"] +=100
        
        for mol in molecules:
            for i in range(len(species[mol])):
                if species[mol][i] not in x[mol].keys():
                    del species[mol][i]
        remove=[]
        self.progress["value"] +=100
        
        nu_list={}
        for mol in molecules:
            for spec in species[mol]:
                nu_list[mol]=list(set([y  for spec in species[mol] for y in x[mol][spec]['NuPosS'] ]))
        
        for mol in molecules:
            for mol1 in molecules:
                for nu_cand in nu_list[mol]:
                    if (nu_cand in nu_list[mol1]) and mol1!=mol:
                        remove.append(nu_cand)
        remove=list(set(remove))
        self.progress["value"] +=100
        
        for mol in molecules:
            for spec in species[mol]:
                x[mol][spec]=comp.delete_dict_list(x[mol][spec],remove)
        
        for mol in molecules:
                if len([x[mol][spec]['Area'][i] for spec in species[mol] for i in range(len(x[mol][spec]['Area'])) ])<3:
                    del x[mol]
        molecules=[mol for mol in molecules if mol in x.keys()]
        species={key:val for key,val in species.items() if key in molecules}
        self.progress["value"] +=100
        
       
        plotting_mol={}
        for mol in molecules:
            plotting_mol[mol]=mol.strip().replace(' ','\,').replace('13C',r'^{13}C').replace('H2','H_{2}').replace('H3','H_{3}').replace('H4','H_{4}').replace('C2','C_{2}').replace(')2',')_{2}')
        plot_names_fc=self.saving_var['Files'][1]['Plotting names'].get()
        if len(plot_names_fc)==0:
            plot_names_fc='plotting_names.json'
            self.saving_var['Files'][1]['Plotting names'].set('plotting_names.json')
        with open(plot_names_fc, 'w') as outfile:
		jdump( plot_names_fc, outfile, encoding='utf-8')
        self.progress["value"] +=100
        
        cleaned_file_candidate=self.saving_var['Files'][1]['Cleaned Spectral lines'].get()        
        if len(cleaned_file_candidate)==0:
            cleaned_file_candidate=file_candidate.split('.')
            cleaned_file_candidate=cleaned_file_candidate[0]+'-clean.json'
            self.saving_var['Files'][1]['Cleaned Spectral lines'].set(cleaned_file_candidate)
        with open(cleaned_file_candidate, 'w') as outfile:
		jdump( x, outfile, encoding='utf-8')
        y={mol:{} for mol in molecules}
        for self.mol in molecules:
            Eu, Nu, NuE, y[self.mol]['Rot'],y[self.mol]['RotCovar'],y[self.mol]['RotSigmas']=comp.fit_rotational( species[self.mol],x[self.mol])
            fig=figure()
            ax=fig.add_subplot(111)
            plotting.saving_plot(self,Eu, Nu, NuE,y[self.mol],plotting_mol[self.mol],ax,'',0.9,fmt='k.',ecolor='black')
            fig.savefig(str('mols/'+self.mol+'.png'))
            close(fig)
        with open('rotdiag.json', 'w') as outfile:
		jdump( y, outfile, encoding='utf-8')
        y_list={}
        y={}
        for mol in molecules:
            y_list[mol]={}
            for spec in species[mol]:
                y_list[mol][x[mol][spec]['SpecPlot'][2]]=[]
        for mol in molecules:
            for spec in species[mol]:
                y_list[mol][x[mol][spec]['SpecPlot'][2]].append(spec)
        for self.mol in molecules:
          y={}
          fig=figure()
          ax=fig.add_subplot(111)
          i=0
          for spec in y_list[self.mol].keys():
            Eu=[]
            Nu=[]
            NuE=[]
            for spec1 in y_list[self.mol][spec]:
                Eu+=x[self.mol][spec1]['Eu']
                Nu+=x[self.mol][spec1]['Nu/gu']
                NuE+=x[self.mol][spec1]['Nu/guE']
            try:
                y['Rot'],y['RotCovar'],y['RotSigmas']=comp.fit_population(Eu, Nu, NuE)
                Eu, Nu, NuE, y['Rot'],y['RotCovar'],y['RotSigmas']=comp.fit_rotational( y_list[self.mol][spec],x[self.mol])
                
                incr=plotting.saving_plot(self,Eu, Nu, NuE,y,plotting_mol[self.mol],ax,spec,0.9-0.1*i,fmt='.')
                if incr:                
                    i+=1
                
            except (RuntimeError,TypeError):
                pass
          if i>0:
              fig.savefig(str('mols/'+self.mol+'_color.png'))
          close(fig)
        
        self.progress["value"] +=100
        ion()
       
        self.progress["value"] +=100
        for mol in molecules:
            for spec in species[mol]:
                print mol, len(x[mol][spec]['Area'])
    def plotting_window(self):

        x,molecules,species,column_names=comp.read_json(self.saving_var['Files'][1]['Cleaned Spectral lines'].get())
        plo1=Toplevel(master)
        plo1.wm_title('Simple plotting')
        plo=Frame(plo1)
        plo.grid(row=0,column=0)

        val_x=StringVar()
        val_y=StringVar()
        set_mol=StringVar()
        typep=StringVar()
        typep1=StringVar()
        typep2=StringVar()
        typep3=StringVar()
        mol0=['All']
        mol0+=sorted(molecules)
        column_names=sorted(column_names)
        opt_mol=ttk.OptionMenu(plo, set_mol, 'Select Molecule', *mol0 )
        opt_x=ttk.OptionMenu(plo, val_x, 'X', *column_names)
        opt_y=ttk.OptionMenu(plo, val_y, 'Y', *column_names)
        opt_mol.grid(row=0,column=0)
        opt_x.grid(row=0,column=1)
        opt_y.grid(row=0,column=2)     
        
        self.scale=ttk.Scale(plo, from_=5, to=200, orient=HORIZONTAL,command=self.update_scale,value=5)
        
        opt_type=ttk.OptionMenu(plo, typep, 'Histogram', *['Yes', 'No'] )
        opt_type1=ttk.OptionMenu(plo, typep1, 'Contours', *['Yes', 'No'] )
        opt_type2=ttk.OptionMenu(plo, typep2, 'Fit', *['Yes', 'No'] )
        opt_type3=ttk.OptionMenu(plo, typep3, 'Scatter', *['Yes', 'No'] )
        self.label_scale=Label(plo)
        
        opt_type.grid(row=0,column=4)
        opt_type1.grid(row=0,column=5)
        opt_type2.grid(row=0,column=6)
        opt_type3.grid(row=0,column=3)
        self.scale.grid(row=0,column=7)
        self.label_scale.grid(row=0,column=8)
        ttk.Button(plo1,text='Plot',command=lambda:plotting.simple_plot(plo1,set_mol.get(),species,x,val_x.get(),val_y.get() ,typep.get(),typep1.get(),typep2.get(),typep3.get(),int(self.scale.get()),'k.'),width=10).grid(row=1,column=0)
        
    def there_is_a_scale(self,evt):
        if evt=='Yes':
            self.label_scale.grid(row=0,column=7)
            self.scale.grid(row=1,column=0)
        else:
            self.label_scale.grid_forget()
            self.scale.grid_forget()
            
    def update_scale(self,evt):
        self.label_scale.configure(text="Bins: {}".format(int(float(evt))))
        
    def reading_partition_functions(self):
        ### A function for fitting partition functions
        self.progress["value"]=0
        self.progress["maximum"]=100
        
        z_name=self.saving_var['Files'][1]['Partition functions'].get()
        if len(z_name)==0:
            tkMessageBox.showerror(message="No Partition function data chosen!")
            return
        zz=z_name.split('.')
        if zz[1]=='json':
            tkMessageBox.showerror(message="Partition functions have already been read!")
            return
        self.z_file=open(z_name,'r')
        z_name=zz[0]+'.json'
        self.progress["value"]+=100
        
        x=np.asarray([300,225,150,75,37.5,18.75,9.375])
        z_actions={'Fitting':9,'Rotational':5}
        l=1
        partition_functions={}
        self.progress["value"]+=0
        self.progress["maximum"]=20
        for line in self.z_file:
            line=line.split(':')
            if len(line) not in z_actions.values():
                tkMessageBox.showerror(message="Incorrect Partition function in line {}!".format(l))
                return
            l+=1
            
            y=np.asarray([float(line[mm]) for mm in range(1,len(line)-1)])
            mol_weight=float(line[len(line)-1])
            if len(line)==z_actions['Fitting']:
                popt, pcov = curve_fit(comp.z_fit, x, y,[1.0, 0.0001])
                partition_functions[line[0]]={'Fit':popt.tolist(),'Covar':pcov.tolist(),'Stdev':np.sqrt(np.diag(pcov)).tolist()}
            if len(line)==z_actions['Rotational']:
                con=(5.331035*10**6)/np.sqrt(y[0]*y[1]*y[2])
                partition_functions[line[0]]={'Fit':con,'Stdev':1.5*con}
            self.progress["value"]+=1
            
        #save partition_functions
        with open(z_name, 'w') as outfile:
                jdump(partition_functions , outfile, encoding='utf-8')
        self.saving_var['Files'][1]['Partition functions'].set(z_name)
        
    def hist_cut(self):
        ioff()

        plotting_w=Toplevel(master)
        plotting_w.minsize(width=666, height=666)
        pframe=VerticalScrolledFrame(plotting_w)
        pframe.grid(row=0,column=0)
        file_candidate=self.saving_var['Files'][1]['Spectral lines'].get()        
        if len(file_candidate)<1:
            tkMessageBox.showerror(message="No Spectral line file chosen!")
            return
        if ".json" not in file_candidate:
            x,molecules,species,cnames_dict=comp.save_to_json(self.saving_var,self.progress)
        else:
            with open(file_candidate, 'r') as outfile:
                x=jload( outfile, encoding='utf-8')
        y={}
        for mol in x.keys():
            for spec in x[mol].keys():
                for key in x[mol][spec].keys():
                    if key not in ['SpecPlot','NuPosS']:
                        if key not in y.keys():
                            y[key]=[]
                        y[key]+=x[mol][spec][key]
        i=0

        for key in App.changeable_list:
            if key in y.keys():
                fig=figure()
                ax=fig.add_subplot(111)            
                ax.set_xlabel(App.tex_output_dict[key])
                hist(y[key], 50, normed=1)
                canvas = FigureCanvasTkAgg(fig, master=pframe.interior)
                canvas.get_tk_widget().grid(row=i,column=0)
                close(fig)
                i+=1


    
                

master = Tk()
App(master)

# display the menu
master.wm_title("Line analysis")
master.mainloop()
