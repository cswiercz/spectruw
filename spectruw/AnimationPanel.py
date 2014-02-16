import sys
import math
import string
from libhill import *
from Tkinter import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.widgets import Slider
import pylab
import numpy
from Panel import *

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

class AnimationPanel(Panel):

    def __init__(self,sysHashmap, *args, **kwargs):

        Panel.__init__(self, *args, **kwargs)
        self.title('Parameter Animation')

        vecSize = sysHashmap["VectorSize"]
        coeffStr = sysHashmap["Coefficients"]
        numFModes = sysHashmap["FourierModes"]
        period = sysHashmap["Period"]
        muSize = sysHashmap["MuVals"]
        params = sysHashmap["Parameters"]
        
        if len(params) > 0:
            params=params.split(':')

            param_str = params[0]
            start = float(params[1])
            end = float(params[2])
            numsteps = float(params[3])
            paramValues = linspace(start,end,numsteps)

        
        P = 2
        if mod(muSize,2) == 0:
            muSize += -1
            
        MU = linspace(-pi/(period*P),pi/(period*P),muSize)

        result = self.computeSpectrum(vecSize,numFModes,period,muSize,
                                      coeffStr,MU,param_str,paramValues)
        plotptrs = self.createPlot(result,MU,numFModes,P,period,paramValues,param_str)

    def computeSpectrum(root,vecSize,numFModes,period,muSize,coeffStr,MU,param_str,paramValues):

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
            for param in paramValues:
                coeffparam = insert_parameter(coeffStr,param_str,param)
                coeffs = reformat(coeffparam,vecSize)
                B.append(hill_comp(MU,period,numFModes,P,coeffs,R,vecflag))

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
        
    def createPlot(root,specMat,MU,N,P,L,paramValues,param):
        
        plot = plt.figure('Parameter Animation',facecolor='whitesmoke')

        cmap = plt.cm.jet
        numcolors = len(specMat)

        spectrumplot = plot.add_subplot(211)
        spectrumplot.set_title('Animation Window')
                
        spectra = []
        k=0
        for A in specMat:
            spectrum = []
            
            for i in range(len(MU)):
                for j in range(len(A[0][i])):
                    spectrum.append(A[0][i][j])
            spectra.append(spectrum)
            spectrumplot.plot(real(spectrum),imag(spectrum),linestyle='None',
                              color=cmap(1.*k/numcolors),marker='.')
            k = k+1

        colorax = plot.add_axes([0.92, 0.1, 0.03, 0.8])
        cb = matplotlib.colorbar.ColorbarBase(colorax, cmap=cmap, spacing='proportional', 
                                              ticks=paramValues, boundaries=paramValues, 
                                              format='$%.2f$')
        colorax.set_xlabel('$' + param + '$')
        
        sliderplot = plot.add_subplot(212)
        sliderplot.plot(real(spectra[0]),imag(spectra[0]),'b.')
        
        sliderplot.set_ylabel('$\mathrm{Im }\lambda$')


        axparam = plot.add_axes([0.2, 0, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        sparam = Slider(axparam, 'param', paramValues[0], paramValues[-1], 
                        valinit=paramValues[-1])
        def _update(widget):
            sliderplot.cla()
            val = sparam.val
            step = paramValues[1]-paramValues[0]
            index = int((val-paramValues[0])/step)
            sliderplot.plot(real(spectra[index]),imag(spectra[index]),'b.')
            sliderplot.set_xlabel('$\mathrm{Re }\lambda$')
          
        sparam.on_changed(_update)

        spectrumplot.set_xlabel('$\mathrm{Re }\lambda$')
        spectrumplot.set_ylabel('$\mathrm{Im }\lambda$')
        sliderplot.set_xlabel('$\mathrm{Re }\lambda$')
        sliderplot.set_ylabel('$\mathrm{Im }\lambda$')
                    
        plot.show()
        
        return [spectrumplot, plot]


