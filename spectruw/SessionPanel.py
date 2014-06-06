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
from fractions import gcd
import pdb
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

cMidground = 'MediumBlue'
cForeground = 'Aqua'
cReal = 'MediumBlue'
cImag = 'Crimson' #IndianRed

class SessionPanel(Panel):
    gPolar = False
    columnindex = 0
    rowindex = 0
    MU = None
    N = None
    L = None
    P = None
    U = None
    V = None
    cid = None   
    plot_handles = None
    colors = None
    def __init__(self,sysHashmap, *args, **kwargs):

        Panel.__init__(self, *args, **kwargs)
        self.title('SpectruW')

        vecSize = sysHashmap["VectorSize"]
        coeffStr = sysHashmap["Coefficients"]
        numFModes = sysHashmap["FourierModes"]
        period = sysHashmap["Period"]
        muSize = sysHashmap["MuVals"]
        params = sysHashmap["Parameters"]
        
        name = ''
        paramValues = None
        if len(params)>0:
            param_str = params.split(':')

            name = param_str[0]
            start = float(param_str[1])
            end = float(param_str[2])
            numsteps = float(param_str[3])
            paramValues = linspace(start,end,numsteps)
            coeffeqns = insert_parameter(coeffStr,name,start)

        self.P = 2
        if mod(muSize,2) == 0:
            muSize += -1

        self.V = vecSize
        self.L = period
        self.N = numFModes
        self.MU = linspace(-pi/(period*self.P),pi/(period*self.P),muSize)
        self.U = reformat(coeffStr,vecSize)

        result = self.computeSpectrum(vecSize,coeffStr,name,paramValues)

        equations = self.createEquations(vecSize,coeffStr)
        plotptrs = self.createPlot(result,name,paramValues)
        self.createInterface(equations,vecSize,plotptrs[0],plotptrs[1],plotptrs[2], plotptrs[3],result,name,paramValues)
        Tk.mainloop()


    def computeSpectrum(self,vecSize,coeffStr,param_str,paramValues):

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
        
        def calculate(*fcoefs):
            B = []
            if len(param_str)==0:
                coeffs = reformat(coeffStr,vecSize)
                return [hill_comp(self.MU,self.L,self.N,self.P,coeffs,R,vecflag)]
            for param in paramValues:
                coeffparam = insert_parameter(coeffStr,param_str,param)
                print coeffparam

                coeffs = reformat(coeffparam,vecSize)
                B.append(hill_comp(self.MU,self.L,self.N,self.P,coeffs,R,vecflag))
            return B
                
        # If Fourier Coefficients are already present in a file.
        if coefbool:
            pyDat = open("pyDat.dat")
            C = pickle.load(pyDat)
            pyDat.close()
            # Compute Spectrum Depending on local/non-local
            # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
            if len(C) == 2:
                return calculate(C[0],C[1])
            else:
                return calculate(C[0])
        else:
            # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
            return calculate()        
        
    def createEquations(self,vect,coeffs):
        
        eqns = raw_reformat(coeffs,vect)
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


    def createPlot(self,specMat,param_str,paramValues):
        
        s = len(self.U)
        spectra = []
            
        plot = Figure(figsize=(5,4), dpi=100, facecolor='whitesmoke',edgecolor=None)

        spectrumplot = plot.add_subplot(211)
        spectrumplot.set_title('$\mathrm{Spectrum}$')
        for A in specMat:
            spectrum = []
            for i in range(len(self.MU)):
                for j in range(len(A[0][i])):
                    spectrum.append(A[0][i][j])
            spectra.append(spectrum)
        
        spectrumplot.set_xlabel('$\mathrm{Re }\lambda$')
        spectrumplot.set_ylabel('$\mathrm{Im }\lambda$')

        efcnplot = plot.add_subplot(212)
        efcnplot.set_xlabel('$x$')
        efcnplot.set_ylabel('$y(x)$')

        lastind = 0
            
        xs = real(spectra[0])
        ys = imag(spectra[0])

        
        def efcnchange(ax):
            print 'called'
            #Get the range for the new area
            xstart,ystart,xdelta,ydelta = ax.viewLim.bounds
            xend = xstart + xdelta
            yend = ystart + ydelta

            x = arange(xstart,xend,.01)
            self.plotEigenfunction(efcnplot,plot,specMat[0],x)
            
        # Connect for changing the view limits
        # retain a reference to avoid garbage collection
        cid = efcnplot.callbacks.connect('xlim_changed', efcnchange)
        #efcnplot.callbacks.connect('ylim_changed', efcnchange)

        plotWindow = FigureCanvasTkAgg(plot,master=self)    
        self.plotSpectrum(specMat,0,spectrumplot,efcnplot, plotWindow)
        
        plotWindow.show() 

        toolbar = NavigationToolbar2TkAgg(plotWindow,self)
        toolbar.pack(side=TOP,padx=10,pady=10)

        plotWindow.get_tk_widget().pack(side=LEFT,fill=BOTH,expand=True,padx=10,pady=10)

        return [spectrumplot, plot, efcnplot, plotWindow]
    
    def plotSpectrum(self,specMat,ind, spectrumplot, efcnplot, fig):
        
        fig.mpl_disconnect(self.cid)

        spectrumplot.cla()
        spectrum = []
        A = specMat[ind]

        self.colors = [cMidground]*(len(self.MU)*len(A[0][0]))

        for i in range(len(self.MU)):
            for j in range(len(A[0][i])):
                spectrum.append(A[0][i][j])
        self.plot_handles = spectrumplot.scatter(real(A[0]),imag(A[0]),c=self.colors,edgecolor='black',linewidth=0.15,alpha=0.7,picker=True)

        self.plot_handles.set_facecolors(self.colors)


        xs = real(spectrum)
        ys = imag(spectrum)

        #line = spectrumplot.plot(xs,ys, 'bo', picker=5)

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

            s = self.V
            self.columnindex = dataind/(s*(2*self.N+1))
            self.rowindex=dataind%(s*(2*self.N+1))

            selected.set_visible(True)
            selected.set_data(xs[dataind], ys[dataind])

            x = arange(0,self.P*self.L,.01)
            self.plotEigenfunction(efcnplot,fig,specMat[ind],x)
            
        fig.draw()
        self.cid = fig.mpl_connect('pick_event', onclick)


    def createInterface(self,equations,vecSize,spectrumplot,plot,efcnplot,plotWindow,A,param_str,paramValues):

        ind = 0

        def _select_mu(widget):
            if len(param_str)>0:
                p = float(slider.get())

                for i in range(0,len(paramValues)):
                    if paramValues[i] >= p:
                        break
            else:
                i = 0

            m = float(widget.cget('text'))
            intvar = int(self.getvar(widget.cget('variable')))
            itemindex = numpy.where(abs(self.MU - m) < numpy.finfo(numpy.single).eps)
            B = A[i]

            M = itemindex[0][0]
            num = len(B[0][0])

            if (intvar==1):
                self.colors[M*num:M*num+num] = [cForeground]*num                
            else:
                self.colors[M*num:M*num+num] = [cMidground]*num
            self.plot_handles.set_facecolors(self.colors)
            plotWindow.draw()

        def _changePlot(widget):
            
            self.gPolar = int(self.getvar(widget.cget('variable')))
            x = None

            self.plotEigenfunction(efcnplot,plotWindow,A[0],x)

        def _changeSpectrum(widget):

            p = float(slider.get())

            for i in range(0,len(paramValues)):
                if paramValues[i] >= p:
                    break
            
            #ind = numpy.where(paramValues==p)[0][0] #paramValues.index(p)
            ind = i
            self.plotSpectrum(A,ind,spectrumplot,efcnplot,plotWindow)

        Label(master=self,text='Equations:',
              font="Verdana 10 bold").pack(side=TOP,anchor=W)

        for eqn in equations:
            Label(master=self,text=eqn).pack(side=TOP)


        Label(master=self,text='Eigenfunction Plot:',
              font="Verdana 10 bold").pack(side=TOP,anchor=W)

        modVar=IntVar()            
        modbutton = Checkbutton(master=self, text='Polar plot',variable=modVar)
        modbutton.pack(side=TOP)
        modbutton.configure(command=lambda cbutton=modbutton:_changePlot(cbutton))
        
        Label(master=self, text='Mu Values:',
              font="Verdana 10 bold").pack(side=TOP,anchor='w')

        canvas = Canvas(self,width=150, borderwidth=0)
        frame = Frame(canvas)#, background='white')
        vsb = Scrollbar(self, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side=LEFT, fill=Y)
        canvas.pack(side=LEFT, fill=Y, expand=True)
        canvas.create_window((4,2), width=100,window=frame, anchor="nw", 
                             tags="frame")

        if len(param_str):
            res = (paramValues[-1]-paramValues[0])/len(paramValues)
            paramvar = DoubleVar()
            slider = Scale(self, from_=paramValues[0], to =paramValues[-1],resolution=res,variable=paramvar,command=_changeSpectrum)
            slider.pack(side=LEFT, fill=Y)
        
        frame.bind("<Configure>", lambda event:canvas.configure(scrollregion=canvas.bbox("all")))
        i = 0
        for muval in self.MU:
            mu_selected = IntVar()
            cb=Checkbutton(master=frame, text=str(muval), variable=mu_selected)
            #cb.pack(side=TOP,fill=BOTH,anchor='w')
            cb.grid(row=i,column=0, sticky=W,columnspan=2)
            i = i+1
            cb.configure(command=lambda widget=cb:_select_mu(widget))

    def plotEigenfunction(self,efcnplot,plot,specMat,x):
        columnindex = self.columnindex
        rowindex = self.rowindex
        
        N = len(self.MU)
        n = self.columnindex

        if (n == N/2): #corresponds to mu=0
            s = 1
        else:
            N = N-1
            g = gcd(2*N,2*n-N)
            s = 2*N/g


        x = arange(0,abs(self.L*self.P*s),.01)            
     
        
        y = numpy.zeros(len(x),dtype=numpy.complex128)
        
        for k in range(0,len(x)):
            m = 0
            
            for j in range(-self.N,self.N+1):
                f_exp = x[k]*(self.MU[columnindex]+pi*j/self.L)   

                y[k] += specMat[2][columnindex][rowindex][m] * cmath.exp(1j*f_exp)
                m += 1

        efcnplot.cla()
        if self.gPolar == False:
            efcnplot.plot(x,numpy.real(y),cReal,label='real')
            efcnplot.plot(x,numpy.imag(y),cImag,label='imag')
            efcnplot.legend(fontsize=10)

            efcnplot.set_xlabel('$x$')
            efcnplot.set_ylabel('$y(x)$')

        else:            
            efcnplot.plot(numpy.real(y),numpy.imag(y),cReal)
            efcnplot.set_xlabel('$\mathrm{Re}$ $y(x)$')
            efcnplot.set_ylabel('$\mathrm{Im}$ $y(x)$')
        
        currentEval = specMat[0][columnindex][rowindex]
        efcnplot.text(0.05, 0.9, '$\mu=%1.3f$\n$\lambda=%1.3f + %1.3fi$' %(self.MU[columnindex], currentEval.real, currentEval.imag), transform=efcnplot.transAxes, va='top')
        plot.draw()

