#
# Copyright 2013 Thomas Trogdon
#                Krithika Manohar
# This software is distributed under the terms of the GNU General Public License
#

import pdb

import sys
import string
from libhill import *
from Tkinter import *
from Panel import *
from SessionPanel import *
from AnimationPanel import *

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

class EntryPanel(Panel):

    def __init__(self, *args, **kwargs):
           
        def _enterSession(currentrow):

            Label(master=sessionDialog,text='System of Equations:')

            # Vector Size
            V = int(vbox.get())
            # Number of Fourier Modes
            N = int(nbox.get())
            # Period Length
            L = eval(lbox.get())
            # Number of Periods
            P = 2; 
            # Number of mu values
            Z = int(zbox.get())
           
            vbox.configure(state='disabled')
            nbox.configure(state='disabled')
            lbox.configure(state='disabled')
            zbox.configure(state='disabled')

            # Boolean for Eigenvectors (force for efcn plotting)
            vecflag = True
            # Boolean for using existing Fourier Coefficients
            coefbool = False 

            def _enterSystem(s):        

                def _enterCoefficients(order):

                    coeffboxes = []

                    # loop through each coefficient box within one equation
                    for i in range(order,-1,-1):
                        existing = coeffvector[hashmap[s]].split(':')
                        derivstr = ''
                        box = Entry(master=eqnDialog)
                        coeffboxes.append(box)
                        
                        if (len(existing) == order+1):
                            box.insert(0,existing[i])
                        else:
                            box.insert(0,1)

                        for j in range(0,i):
                            derivstr = derivstr + "'"

                        if (i == 0):
                            derivstr = 'y' 
                        else:
                            derivstr = 'y'+derivstr+'+'

                        Label(master=eqnDialog, text=derivstr).grid(row=order-i+1,column=2,sticky=W)
                        box.grid(row=order-i+1,column=1)
            
                    def _saveCoefficients():
                        order = len(coeffboxes)-1

                        str_array = []
                        for i in reversed(range(0,len(coeffboxes))):
                            str_array.append(coeffboxes[i].get())

                        coeffvector[hashmap[s]] += ':'.join(str_array)
                        eqnDialog.destroy()

                    Button(master=eqnDialog,text='Continue',command=_saveCoefficients).grid(row=order+2,column=0)
                    Button(master=eqnDialog,text='Change Order',command=lambda:_changeOrder()).grid(row=order+2,column=1)

                    def _changeOrder():
                        _enterSystem(s)
                        coeffvector[hashmap[s]]=''
                        eqnDialog.destroy()


                eqnDialog = Panel()
                eqnDialog.wm_title(s+": Enter Coefficients")
                if (len(coeffvector) and len(coeffvector[hashmap[s]])):
                    s_coeff = coeffvector[hashmap[s]].split(':')
                    _enterCoefficients(len(s_coeff)-1)
                else:
                    Label(master=eqnDialog, text='Maximum Differential Order:').grid(row=0)
                    orderstr =  StringVar()
                    orderbox = Entry(master=eqnDialog, textvariable=orderstr)
                    orderbox.insert(0, "0")
                    orderbox.grid(row=0,column=1)
                    b = Button(master=eqnDialog, text='Setup Equation',command=lambda:_enterCoefficients(int(orderbox.get())))
                    b.grid(row=0,column=2)
                    b.configure(command=lambda:[_enterCoefficients(int(orderbox.get())), b.destroy(), orderbox.configure(state='disabled')])

      
            coeffvector = []
            hashmap = dict()
            counter = 0

            Label(master=sessionDialog, text='Enter Equations:').grid(row=currentrow,column=0)
            Label(master=sessionDialog, text=u' = \u03bb y').grid(row=currentrow+2,column=V+1)

            for i in range(0,V):
                currentrow = currentrow+1
                for j in range(0,V):
                    coeffvector.append('') 
                    txt = 'L['+str(i+1)+','+str(j+1)+']'
                    hashmap[txt] = counter
                
                    b = Button(master=sessionDialog, text=txt)
                    b.grid(row=1+currentrow,column=j)
                    b.configure(command=lambda widget=b.cget('text'): _enterSystem(widget))
                    counter += 1


            currentrow = currentrow+1

            def _plotSession(Z):
                coeffstr = ';'.join(coeffvector)

                param = paramEntry.get()
                if len(param) > 0:
                    start = startEntry.get()
                    end = endEntry.get()
                    steps = stepsEntry.get()
                    param = param+':'+start+':'+end+':'+steps
                
                debug = False
                sysHashmap = { "VectorSize" : V, 
                               "FourierModes" : N, 
                               "Period" : L, 
                               "NumPeriods" : P, 
                               "MuVals" : Z , 
                               "Coefficients" : coeffstr ,
                               "Parameters" : param}
                print coeffstr
                sys.stdout.flush()
                if debug:
                    sysHashmap = { "VectorSize" : 2, 
                                   "FourierModes" : 10, 
                                   "Period" : pi, 
                                   "NumPeriods" : 2, 
                                   "MuVals" : 50 , 
                                   "Coefficients" : "sqrt(1+v0)*sin(2*x);-0.5+(1+v0)*(1-cos(2*x)):0:-0.5;0.5-(1+cos(2*x)):0:0.5;-sqrt(1+v0)*sin(2*x)",
                                   "Parameters" : "v0:-1:1:10"}

                    # sysHashmap = { "VectorSize" : 1, 
                    #                "FourierModes" : 10, 
                    #                "Period" : pi, 
                    #                "NumPeriods" : 2, 
                    #                "MuVals" : 50 , 
                    #                "Orders" : "4", 
                    #                "Coefficients" : "(1-a^2)*sin(x)*cos(x):a:-sin(x):0:-1",
                    #                "Parameters" : "a:0:1:10" }

                sessionDialog.destroy()
#                 if len(param) > 0:
#                     anipanel = AnimationPanel(sysHashmap)
#                     anipanel.withdraw()
                SessionPanel(sysHashmap)

            currentrow = currentrow + 1
            Label(master=sessionDialog, text='Parameter name:').grid(row=currentrow)
            paramEntry = Entry(master=sessionDialog)
            paramEntry.grid(row=currentrow,column=1)

            currentrow = currentrow + 1
            Label(master=sessionDialog, text='Start value:').grid(row=currentrow)
            startEntry = Entry(master=sessionDialog)
            startEntry.grid(row=currentrow,column=1)

            currentrow = currentrow + 1
            Label(master=sessionDialog, text='End value:').grid(row=currentrow)
            endEntry = Entry(master=sessionDialog)
            endEntry.grid(row=currentrow,column=1)

            currentrow = currentrow + 1
            Label(master=sessionDialog, text='Number of values:').grid(row=currentrow)
            stepsEntry = Entry(master=sessionDialog)
            stepsEntry.grid(row=currentrow,column=1)
                        
            Button(master=sessionDialog, text='Plot Spectrum', 
                   command=lambda:_plotSession(Z)).grid(row=V+1+currentrow)

        # System parameter window
        sessionDialog = Panel()
        sessionDialog.wm_title("System parameters")

        Label(master=sessionDialog, text='Vector size:').grid(row=0)
        vstr = StringVar()
        vbox = Entry(master=sessionDialog, textvariable=vstr)
        vbox.insert(0, "1")
        vbox.grid(row=0,column=1)
        
        Label(master=sessionDialog, text='Number of Fourier Modes:').grid(row=1,sticky=W)
        nbox = Entry(master=sessionDialog)
        nbox.insert(0,"10")
        nbox.grid(row=1,column=1)
        
        Label(master=sessionDialog, text='Period Length:').grid(row=2,sticky=W)
        lbox = Entry(master=sessionDialog)
        lbox.insert(0,"2*pi")
        lbox.grid(row=2,column=1)
        
    # Label(master=sessionDialog, text='Number of Periods').grid(row=3,sticky=W)
    # pbox = Entry(master=sessionDialog)
    # pbox.insert(0,"2")
    # pbox.grid(row=3,column=1)

        Label(master=sessionDialog, text='Number of mu values:').grid(row=4,sticky=W)
        zbox = Entry(master=sessionDialog)
        zbox.insert(0,"50")
        zbox.grid(row=4,column=1)
        row = 5
    # vecvar = IntVar()
    # vecbutton = Checkbutton(master=sessionDialog, text='Eigenvectors', variable=vecvar).grid(row=5,sticky=W)
    # coefvar = IntVar()
    # coefbutton = Checkbutton(master=sessionDialog, text='Use existing Fourier coefficients',variable=coefvar).grid(row=6,sticky=W)
    
        bContinue = Button(master=sessionDialog, text='Continue')
        bContinue.grid(row=5,column=1)
        bContinue.configure(command=lambda:[_enterSession(row+1), bContinue.destroy()])

        Tk.mainloop()

if __name__ == '__main__':
    EntryPanel()
