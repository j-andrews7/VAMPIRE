#!/usr/bin/env python

#For gui
import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askdirectory
from tkinter import messagebox

#For command line
import os

class Tfsites_checker_options():
    def __init__(self):
        self.r = None
        self.i = None
        self.m = None
        self.o = None
        
        #Optional arguments
        self.ci = None
        #not implemented self.co = None
        self.bp = None
        
        # filter arguments
        self.fm = True
        self.fc = False
        self.fn = False
        
    def missing_args(self):
        return (self.r == None or self.i == None or self.m == None or
            self.o == None)
        

def choose_o_dir():
    dirname = askdirectory(title="Select output VCF location")
    checker_op.o = dirname
    checker_o.set("Output directory: "+final_fn(dirname))
    
def choose_i_dir():
    dirname = askdirectory(title="Select file containing VCFs to be processed")
    checker_op.i = dirname
    checker_i.set("Input VCF directory: "+final_fn(dirname))
    
def choose_m():
    filename = askopenfilename(title="Select motif file (JASPAR file format required)",
        filetypes=[('Text file', '.txt'),('All file types', '.*')])
    checker_op.m = filename
    checker_m.set("Motif list file: "+final_fn(filename))
    
def choose_r():
    filename = askopenfilename(title="Select reference file",
        filetypes=[('Fasta file', '.fa'),('All file types', '.*')])
    checker_op.r = filename
    checker_r.set("Reference sequence file: "+final_fn(filename))

def choose_ci(strVar):
    filename = askopenfilename(title="Select ChIP input file",
        filetypes=[('Merged TF ChIP peaks', '.bed'),('All file types', '.*')])
    checker_op.ci = filename
    strVar.set("Transcription Factor ChIP peaks (Merged .bed file:\n"+
        final_fn(filename))   

def choose_bp(strVar):
    filename = askopenfilename(
        title="Select baseline nucleotide frequency file",
        filetypes=[('Text file', '.txt'),('All file types', '.*')])
    checker_op.bp = filename
    strVar.set("Baseline nucleotide probabilities file:\n"+
        final_fn(filename))
        
def final_fn(filename):
    if filename == None:
        return "None"
    output = filename.replace('\\','/').split('/')[-1]
    if output == "":
        return "None"
    return output


    
def optional_window():
    master.withdraw()
    
    root = Toplevel()
    
    #Makes sure the previous window opens again if the options window is closed
    root.protocol("WM_DELETE_WINDOW", lambda: close_options(root))
    
    root.wm_title("NEED A NAME Optional Arguments")
    
    #Variables to display current files selected
    checker_ci = StringVar()
    checker_ci.set("Transcription Factor ChIP peaks (Merged .bed file:\n"+
        final_fn(checker_op.ci))
        
    checker_bp = StringVar()
    checker_bp.set("Baseline nucleotide probabilities file:\n"+
        final_fn(checker_op.bp))
    
    #Variables for check boxes with filter options
    checker_fm = BooleanVar()
    checker_fc = BooleanVar()
    checker_fn = BooleanVar()
    
    header = Label(root,text="Tfsites_checker.py optional arguments:")
    header.grid(row=0,column=0,sticky=W,pady=10)
    
    bci = Button(root, text='Choose ChIP input file', 
        command= lambda: choose_ci(checker_ci), width=25)
    bci.grid(row=2,column=0,sticky=W,padx=5,pady=2)
    
    Label(root, textvariable=checker_ci).grid(row=2,column=1,sticky=W,padx=5)
    
    bbp = Button(root, text='Choose baseline probabilities file', 
        command= lambda: choose_bp(checker_bp), width=25)
    bbp.grid(row=3,column=0,sticky=W,padx=5,pady=2)
    
    Label(root, textvariable=checker_bp).grid(row=3,column=1,sticky=W,padx=5)    
    
    bfc = Checkbutton(root, text="Filter VCF output by ChIP peak\n"+
        "(Only output variants that match a ChIP peak)", variable=checker_fc,
        onvalue=True, offvalue=False)
    if checker_op.fc:
        bfc.select()
    bfc.grid(row=4, column=1, sticky=W, padx=5,pady=2)
    
    bfn = Checkbutton(root, text="Filter VCF output by no ChIP peak\n"+
        "(Only output variants that do not match a ChIP peak)", variable=checker_fn,
        onvalue=True, offvalue=False)
    if checker_op.fn:
        bfn.select()
    bfn.grid(row=5, column=1, sticky=W, padx=5,pady=2)
    
    bfm = Checkbutton(root, text="Filter VCF output by motif match\n"+
        "(Only output variants that match a motif)", variable=checker_fm,
        onvalue=True, offvalue=False)
    #The button should be checked if it is present in
    # the Tfsites_checker_options object
    if checker_op.fm:
        bfm.select()
    else:
        bfm.deselect()
        bfc.config(state=DISABLED)
        bfn.config(state=DISABLED)
    bfm.grid(row=4, column=0, sticky=W, rowspan=2, padx=5,pady=2)

    bfm.config(command = lambda: select_fm(checker_fm, bfm, bfc, bfn))
    bfc.config(command = lambda: select_fc(checker_fc, bfc, bfn))
    bfn.config(command = lambda: select_fn(checker_fn, bfc, bfn))
    
    b_close = Button(root,text=' Back ',
        command= lambda: close_options( root), width=10)
    b_close.grid(row=7,column=0,pady=5)
    
    root.focus_force()
    

def select_fm(boolvar, bfm, bfc, bfn):
    #If the button is selected
    if boolvar.get():
        checker_op.fm = True
        bfc.config(state=ACTIVE)
        bfn.config(state=ACTIVE)
    else:
        checker_op.fm = False
        checker_op.fc = False
        checker_op.fn = False
        bfc.deselect()
        bfn.deselect()
        bfc.config(state=DISABLED)
        bfn.config(state=DISABLED)

def select_fc(boolvar, bfc, bfn):
    #If the button is selected (on)
    if boolvar.get():
        checker_op.fc = True
        checker_op.fn = False
        bfn.deselect()
    else:
        checker_op.fc = False
        checker_op.fn = True
    
def select_fn(boolvar, bfc, bfn):
    #If the button is selected (on)
    if boolvar.get():
        checker_op.fc = False
        checker_op.fn = True
        bfc.deselect()
    else:
        checker_op.fc = True
        checker_op.fn = False
    
def close_options(window):
    window.destroy()
    master.deiconify()  
    #debug print(checker_op.fm)
    #debug print(checker_op.fc)
    #debug print(checker_op.fn)
    sys.stdout.flush()
    
    
    
    
def run_tfsites_checker():
    #Check that all required arguments are present
    if checker_op.missing_args():
        messagebox.showwarning("Missing required arguments",
            "Please enter all required arguments")
        return
    
    master.quit()
    master.destroy()
    
    cmd_str = "python tfsites_checker.py -r "+checker_op.r+" -m "+checker_op.m
    
    if checker_op.fm:
        cmd_str += " -fm"
    if checker_op.fc:
        cmd_str += " -fc"
    elif checker_op.fn:
        cmd_str += " -fn"
    
    if checker_op.bp != None:
        cmd_str += " -bp "+checker_op.bp
    
    if checker_op.ci != None:
        cmd_str += " -ci "+checker_op.ci
    
    #Easier to use variables
    vcf_i_dir_path = checker_op.i
    vcf_o_dir_path = checker_op.o
    
    #Set of (string) command lines
    cmd_lines = []
    
    for filename in os.listdir(vcf_i_dir_path):
        vcf_i_path = os.path.join(vcf_i_dir_path,filename)
        if os.path.isfile(vcf_i_path) and filename.endswith(".vcf"):
            vcf_o_path = os.path.join(vcf_o_dir_path,filename)
            cmd_add = " -i "+os.path.normpath(vcf_i_path)
            cmd_add += " -o "+os.path.normpath(vcf_o_path)
            cmd_lines.append( cmd_str + cmd_add )

    print("VCFs to read: "+str(len(cmd_lines)))
    sys.stdout.flush()
    
    for cmd_line in cmd_lines:
        #print(cmd_line)
        #sys.stdout.flush()
        os.system(cmd_line)
        #break
        
            
        
    
    
from tkinter import *
master = Tk()
master.wm_title("NEED A NAME")

#Options class to be used when calling tfsites_checker
checker_op = Tfsites_checker_options()

#String variables to be displayed in the window
checker_r = StringVar()
checker_r.set("Reference sequence file: None")
checker_m = StringVar()
checker_m.set("Motif list file: None")
checker_i = StringVar()
checker_i.set("Input VCF directory: None")
checker_o = StringVar()
checker_o.set("Output directory: None")
   
Label(master,text="Tfsites_checker.py arguments:").grid(row=0,sticky=W,pady=10)
#var1 = IntVar()
#Checkbutton(master, text="male", variable=var1).grid(row=1, sticky=W)
#var2 = IntVar()
#Checkbutton(master, text="female", variable=var2).grid(row=2, sticky=W)


br = Button(master, text='Choose reference file', command=choose_r, width=25)
br.grid(row=2,sticky=W,padx=5,pady=2)
Label(master, textvariable=checker_r).grid(row=2,column=1,sticky=W,padx=5)

bm = Button(master, text='Choose motif list', command=choose_m, width=25)
bm.grid(row=3,sticky=W,padx=5,pady=2)
Label(master, textvariable=checker_m).grid(row=3,column=1,sticky=W,padx=5)

bi = Button(master, text='Choose input VCF directory', command=choose_i_dir)
bi.config(width=25)
bi.grid(row=4,sticky=W,padx=5,pady=2)
Label(master, textvariable=checker_i).grid(row=4,column=1,sticky=W,padx=5)

bo = Button(master, text='Choose output directory', command=choose_o_dir)
bo.config(width=25)
bo.grid(row=5,sticky=W,padx=5,pady=2)
Label(master, textvariable=checker_o).grid(row=5,column=1,sticky=W,padx=5)

b_options = Button(master,text='Optional arguments',command=optional_window)
b_options.config(width=20)
b_options.grid(row=6,column=0,pady=10,padx=5)

b_quit = Button(master,text='Exit',command=master.quit, width=10)
b_quit.grid(row=7,column=1,sticky=E, pady=5, padx=5)

b_run = Button(master,text=' Run ',command=run_tfsites_checker, width=10)
b_run.grid(row=7,column=0,pady=5)

mainloop()