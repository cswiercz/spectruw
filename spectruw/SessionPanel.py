import math
import string
from libhill import *
from Tkinter import *
import matplotlib, sys
import matplotlib.pyplot as plt
import numpy
import cmath
from Panel import *
from numpy import arange, sin, pi
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

class SessionPanel(Panel):
    gPolar = False
    columnindex = 0
    rowindex = 0

    def __init__(self,sysHashmap, *args, **kwargs):

        Panel.__init__(self, *args, **kwargs)
        self.title('SpectruW')

        vecSize = sysHashmap["VectorSize"]
        coeffStr = sysHashmap["Coefficients"]
        numFModes = sysHashmap["FourierModes"]
        period = sysHashmap["Period"]
        muSize = sysHashmap["MuVals"]
        params = sysHashmap["Parameters"]
        
        if len(params):
            param_str = params.split(':')
            name = param_str[0]
            start = float(param_str[1])
            coeffStr = insert_parameter(coeffStr,name,start)

        P = 2
        if mod(muSize,2) == 0:
            muSize += -1
            
        MU = linspace(-pi/(period*P),pi/(period*P),muSize)

        result = self.computeSpectrum(vecSize,numFModes,period,muSize,reformat(coeffStr,vecSize),MU)
        equations = self.createEquations(vecSize,coeffStr)
        plotptrs = self.createPlot(result,MU,numFModes,P,period)
        self.createInterface(equations,vecSize,MU,numFModes,P,period,plotptrs[0],plotptrs[1],plotptrs[2], result)
        
        Tk.mainloop()


    def computeSpectrum(root,vecSize,numFModes,period,muSize,coeffStr,MU):

        # To do: Add support for non-local coefficients
        R = None

        # Boolean for Eigenvectors (force for efcn plotting)
        vecflag = True
        # Boolean for using existing Fourier Coefficients
        coefbool = False 

        # Number of Periods
        P = 2;

        # Initial variable for Fourier Coefficients
        C = None
        
        # If Fourier Coefficients are already present in a file.
        if coefbool:
            pyDat = open("pyDat.dat")
            C = pickle.load(pyDat)
            pyDat.close()
            # Compute Spectrum Depending on local/non-local
            # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
            if len(C) == 2:
                A = hill_comp(MU,period,numFModes,P,coeffStr,R,vecflag,C[0],C[1])
            else:
                A = hill_comp(MU,period,numFModes,P,coeffStr,R,vecflag,C[0])
        else:
            # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
            A = hill_comp(MU,period,numFModes,P,coeffStr,R,evecs=vecflag)

            return A
        
        
    def createEquations(root,vect,coeffs):
        
        eqns = reformat(coeffs,vect)
        strs = []
        for eqn in eqns:
            for fcn in eqn:
                c = ''
                order = 0
                for coef in fcn:
                    temp = ''
                    for k in range(0,order): temp = temp+ "'"
                    if len(c) > 0: 
                        plus=' + '
                    else:
                        plus=''
                    c = str(coef)+'y' + temp + plus + c
                    order = order + 1
                strs.append(c+' = 0')

        return strs


    def createPlot(root,specMat,MU,N,P,L):
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
            for j in range(len(specMat[0][i])):
                EigenValuesReal.write(str(real(specMat[0][i][j]))+'\n')
                EigenValuesIm.write(str(imag(specMat[0][i][j]))+'\n')        
                spectrumx.append(real(specMat[0][i][j]));
                spectrumy.append(imag(specMat[0][i][j]));
                        
                #Efficient way to turn vector into delimited string
                EigenVectorsReal.write(string.join([`num` for num in real(specMat[2][i][:,j])],'\t')+'\n')
                EigenVectorsIm.write(string.join([`num` for num in imag(specMat[2][i][:,j])],'\t')+'\n')

        plot = Figure(figsize=(5,4), dpi=100, facecolor='whitesmoke',edgecolor=None)

        spectrumplot = plot.add_subplot(211)
        spectrumplot.set_title('$\mathrm{Spectrum}$')
        line = spectrumplot.plot(spectrumx,spectrumy, 'bo', picker=5) #5 points tolerance
        spectrumplot.set_xlabel('$\mathrm{Re }\lambda$')
        spectrumplot.set_ylabel('$\mathrm{Im }\lambda$')

        efcnplot = plot.add_subplot(212)
        efcnplot.set_xlabel('$x$')
        efcnplot.set_ylabel('$y(x)$')
                
        lastind = 0
        xs = spectrumx
        ys = spectrumy
        selected,  = spectrumplot.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.75,
                                       color='lightgoldenrodyellow', visible=False)

        def onclick(event):
            n = len(event.ind)
            if not n: return True
        
            # the click locations
            x = event.mouseevent.xdata
            y = event.mouseevent.ydata
            distances = numpy.hypot(x-xs[event.ind[0]], y-ys[event.ind[0]])
            indmin = distances.argmin()
            dataind = event.ind[indmin]
            lastind = dataind

            root.columnindex = dataind/(2*N+1)
            root.rowindex=dataind%(2*N+1)

            selected.set_visible(True)
            selected.set_data(xs[dataind], ys[dataind])

            root.plotEigenfunction(efcnplot,plot,specMat,MU,N,P,L)
            
            

        plotWindow = FigureCanvasTkAgg(plot,master=root)    
        plotWindow.show() 

        toolbar = NavigationToolbar2TkAgg(plotWindow,root)
        toolbar.pack(side=TOP)

        plotWindow.get_tk_widget().pack(side=LEFT,fill=BOTH,expand=True)
        plotWindow.mpl_connect('pick_event', onclick)

        return [spectrumplot, plot, efcnplot]


    def createInterface(root,equations,vecSize,MU,N,P,L,spectrumplot,plot,efcnplot,A):

        def _select_mu(widget):
            m = float(widget.cget('text'))
            intvar = int(root.getvar(widget.cget('variable')))
            itemindex = numpy.where(abs(MU - m) < numpy.finfo(numpy.single).eps)
            if (intvar==1):
                spectrumplot.plot(real(A[0][itemindex[0][0]]),imag(A[0][itemindex[0][0]]),'o',color='skyblue',picker=5)
            else:
                spectrumplot.plot(real(A[0][itemindex[0][0]]),imag(A[0][itemindex[0][0]]),'bo',picker=5)
            plot.canvas.draw()
        
        def _changePlot(widget):
            
            root.gPolar = int(root.getvar(widget.cget('variable')))
            root.plotEigenfunction(efcnplot,plot,A,MU,N,P,L)
                            
        Label(master=root,text='Equations:').pack(side=TOP,anchor='w')

        for eqn in equations:
            Label(master=root,text=eqn).pack(side=TOP)

        # Button(master=root, text='New').pack(side=TOP,anchor='w')
        # Button(master=root, text='Close').pack(side=TOP,anchor='w')

        Label(master=root,text='Eigenfunction Plot:').pack(side=TOP,anchor='w')

        modVar=IntVar()            
        modbutton = Checkbutton(master=root, text='Phase-plane plot',variable=modVar)
        modbutton.pack(side=TOP)
        modbutton.configure(command=lambda cbutton=modbutton:_changePlot(cbutton))
        
        Label(master=root, text='Mu Values:').pack(side=TOP,anchor='w')

        canvas = Canvas(root,width=150, borderwidth=0)
        frame = Frame(canvas, background='white')
        vsb = Scrollbar(root, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side=LEFT, fill=Y)
        canvas.pack(side=LEFT, fill=BOTH, expand=True)
        canvas.create_window((4,2), width=100,window=frame, anchor="nw", 
                             tags="frame")
        
        frame.bind("<Configure>", lambda event:canvas.configure(scrollregion=canvas.bbox("all")))

        for muval in MU:
            mu_selected = IntVar()
            cb=Checkbutton(master=frame, text=str(muval), variable=mu_selected)
            cb.pack(side=TOP,fill=BOTH,anchor='w')
            cb.configure(command=lambda widget=cb:_select_mu(widget))

    def plotEigenfunction(root,efcnplot,plot,specMat,MU,N,P,L):
        columnindex = root.columnindex
        rowindex = root.rowindex
        x = arange(0,2*P*L,.01)
        y = numpy.zeros(len(x),dtype=numpy.complex128)
        
        for k in range(0,len(x)):
            m = 0
            
            for j in range(-N,N+1):
                f_exp = x[k]*(MU[columnindex]+pi*j/L)                
                y[k] += specMat[2][columnindex][rowindex][m] * cmath.exp(1j*f_exp)
                m += 1

        efcnplot.cla()
        if root.gPolar == False:
            efcnplot.plot(x,numpy.real(y),label='real')
            efcnplot.plot(x,numpy.imag(y),'r',label='imag')
            efcnplot.legend(fontsize=10)

            efcnplot.set_xlim(0, P*L)
            efcnplot.set_xlabel('$x$')
            efcnplot.set_ylabel('$y(x)$')

        else:
            efcnplot.plot(numpy.real(y),numpy.imag(y))
            efcnplot.set_xlabel('$\mathrm{Re}$ $y(x)$')
            efcnplot.set_ylabel('$\mathrm{Im}$ $y(x)$')
        
        currentEval = specMat[0][columnindex][rowindex]
        efcnplot.text(0.05, 0.9, '$\mu=%1.3f$\n$\lambda=%1.3f + %1.3fi$' %(MU[columnindex], currentEval.real, currentEval.imag), transform=efcnplot.transAxes, va='top')
        
        plot.canvas.draw()

#    def saveSessionState():
