#
# Copyright 2013 Thomas Trogdon
#                Krithika Manohar
# This software is distributed under the terms of the GNU General Public License
#

import sys
import math
import string
from libhill import *
from Tkinter import *
import pickle
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pylab
import matplotlib.pyplot as plt
from PointBrowser import PointBrowser

from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler

from matplotlib.figure import Figure

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

def _continue():
   
    # TO DO: exception handling for invalid entries
    
    # Vector Size
    V = int(vbox.get())
    # Number of Fourier Modes
    N = int(nbox.get())
    # Period Length
    L = eval(lbox.get())
    # Number of Periods
    P = int(pbox.get())
    # Number of mu values
    Z = int(zbox.get())
    # Boolean for Eigenvectors
    vecflag = bool(vecvar.get()) 
    # Boolean for using existing Fourier Coefficients
    coefbool = bool(coefvar.get())
    w_param.destroy()

    def _entercoeffmat(s):        

        def _entercoeffs():
            order = int(orderbox.get())

            w_coeff = Tk.Tk()
            w_coeff.wm_title("Enter coefficients")
            coeffboxes = []

            for i in range(order,-1,-1):
                derivstr = ''
                box = Entry(master=w_coeff)
                coeffboxes.append(box)
                box.insert(0, str(i))
                
                for j in range(0,i):
                    derivstr = derivstr + "'"

                if (i == 0):
                    derivstr = u'y = \u03bby' # unicode greek letter lambda
                else:
                    derivstr = 'y'+derivstr+'+'

                Label(master=w_coeff, text=derivstr).grid(row=order-i,column=1,sticky=W)
                box.grid(row=order-i)
            
            def _saveCoeffs():
                str_array = []
                for i in reversed(range(0,len(coeffboxes))):
                    str_array.append(coeffboxes[i].get())

                coeffvector[hashmap[s]] += ':'.join(str_array)
                w_coeff.destroy()
                w_dorder.destroy()

            Tk.Button(master=w_coeff,text='OK',command=_saveCoeffs).grid(row=order+1)
            
        w_dorder = Tk.Tk()
        w_dorder.wm_title("Coefficient equation")
        Label(master=w_dorder, text='Differential order:').grid(row=0)
        orderstr =  StringVar()
        orderbox = Entry(master=w_dorder, textvariable=orderstr)
        orderbox.insert(0, "2")
        orderbox.grid(row=0,column=1)
        Tk.Button(master=w_dorder, text='OK',command=_entercoeffs).grid(row=1,column=2)

    def _plot(Z):                                

        coeffstr = ';'.join(coeffvector)
        print(coeffstr)

        #U = reformat('(1-.26^2)*sin(x)*cos(x):.26:-sin(x):0:-1',V)
        # String in quotes containing coefficients, then formatted
        U = reformat(coeffstr,V)
            
        # Initial variable for Fourier Coefficients
        C = None

        # Look to see if non-local coefficients are present, and format it
        try:
            R = reformat(sys.argv[9],V)
            print 'Non-Local'
        except:
            R = None
            print 'Local'

        if mod(Z,2) == 0:
            Z += -1

        MU = linspace(-pi/(L*P),pi/(L*P),Z)


        # If Fourier Coefficients are already present in a file.
        if coefbool:
            pyDat = open("pyDat.dat")
            C = pickle.load(pyDat)
            pyDat.close()
            # Compute Spectrum Depending on local/non-local
            # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
            if len(C) == 2:
                A = hill_comp(MU,L,N,P,U,R,vecflag,C[0],C[1])
            else:
                A = hill_comp(MU,L,N,P,U,R,vecflag,C[0])
        else:
            # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
            A = hill_comp(MU,L,N,P,U,R,evecs=vecflag)


        # Open all storage files
        EigenValuesReal = open("EigenValuesReal.dat","w")
        EigenValuesIm = open("EigenValuesIm.dat","w")
        EigenVectorsReal = open("EigenVectorsReal.dat","w")
        EigenVectorsIm = open("EigenVectorsIm.dat","w")
        mu = open("mu.dat","w")

        spectrumx = [];
        spectrumy = [];

        print 'Writing Data'
        for i in range(len(MU)):
            mu.write(str(MU[i])+'\n') #should this come after 'for i...'
            for j in range(len(A[0][i])):
                EigenValuesReal.write(str(real(A[0][i][j]))+'\n')
                EigenValuesIm.write(str(imag(A[0][i][j]))+'\n')        
                spectrumx.append(real(A[0][i][j]));
                spectrumy.append(imag(A[0][i][j]));
                        
                #Efficient way to turn vector into delimited string
                EigenVectorsReal.write(string.join([`num` for num in real(A[2][i][:,j])],'\t')+'\n')
                EigenVectorsIm.write(string.join([`num` for num in imag(A[2][i][:,j])],'\t')+'\n')


        # If Fourier coefficients were not given
        if C is None:
            print 'Writing Coefficients'
            pyDat = open("pyDat.dat","w")
            if R is not None:
                pickle.dump([A[1],A[3]],pyDat)
            else:
                pickle.dump([A[1]],pyDat)
        
        # All files automagically closed.
        def _quit():
            w_spectrum.quit()     # stops mainloop
            w_spectrum.destroy()  # this is necessary on Windows to prevent
            # Fatal Python Error: PyEval_RestoreThread: NULL tstate
                    
        w_spectrum = Tk.Tk()
        w_spectrum.wm_title("Spectrum")
                
        # TO DO: parametrize plot dimensions
        plot = Figure(figsize=(5,4), dpi=100)
        plot.add_subplot(111).plot(spectrumx,spectrumy,'b.')


        #####################################        
        # plot = plt.figure()
        # ax = plot.add_subplot(211)
        # ax.set_title('click on point to plot time series')
        # line, = ax.plot(spectrumx, spectrumy, 'o', picker=5) #5 points tolerance
        # ax2 = plot.add_subplot(212)

        # browser = PointBrowser(ax,ax2,spectrumx,spectrumy,line,A[0],plot)

        # plot.canvas.mpl_connect('pick_event', browser.onpick)
        # plot.canvas.mpl_connect('key_press_event', browser.onpress)
        
        # plt.show()                
        #####################################


        canvas = FigureCanvasTkAgg(plot, master=w_spectrum)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
                
        toolbar = NavigationToolbar2TkAgg( canvas, w_spectrum )
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
                
        Tk.Button(master=w_spectrum, text='Quit', command=_quit).pack(side=Tk.BOTTOM)   
        w_matrix.destroy()
      
    w_matrix = Tk.Tk()
    w_matrix.wm_title("Coefficient Matrix")
    coeffvector = []
    hashmap = dict()
    counter = 0

    for i in range(0,V):
        for j in range(0,V):
            coeffvector.append('') 
            txt = 'L['+str(i+1)+','+str(j+1)+']'
            hashmap[txt] = counter

            b = Tk.Button(master=w_matrix, text=txt)
            b.grid(row=i,column=j)
            b.configure(command=lambda widget=b.cget('text'): _entercoeffmat(widget))
            counter += 1

    Tk.Button(master=w_matrix, text='Plot Spectrum', command=lambda:_plot(Z)).grid(row=V+1)

# System parameter window
w_param = Tk.Tk()
w_param.wm_title("System parameters")

Label(master=w_param, text='Vector size:').grid(row=0)
vstr = StringVar()
vbox = Entry(master=w_param, textvariable=vstr)
vbox.insert(0, "1")
vbox.grid(row=0,column=1)

Label(master=w_param, text='Number of Fourier Modes').grid(row=1,sticky=W)
nbox = Entry(master=w_param)
nbox.insert(0,"10")
nbox.grid(row=1,column=1)

Label(master=w_param, text='Period Length').grid(row=2,sticky=W)
lbox = Entry(master=w_param)
lbox.insert(0,"2*pi")
lbox.grid(row=2,column=1)

Label(master=w_param, text='Number of Periods').grid(row=3,sticky=W)
pbox = Entry(master=w_param)
pbox.insert(0,"2")
pbox.grid(row=3,column=1)

Label(master=w_param, text='Number of mu values').grid(row=4,sticky=W)
zbox = Entry(master=w_param)
zbox.insert(0,"50")
zbox.grid(row=4,column=1)

vecvar = IntVar()
vecbutton = Checkbutton(master=w_param, text='Eigenvectors', variable=vecvar).grid(row=5,sticky=W)
coefvar = IntVar()
coefbutton = Checkbutton(master=w_param, text='Use existing Fourier coefficients',variable=coefvar).grid(row=6,sticky=W)
    
Tk.Button(master=w_param, text='Continue', command=_continue).grid(row=8,column=1)

Tk.mainloop()

