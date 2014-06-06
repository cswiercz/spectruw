#
# Copyright 2013 Krithika Manohar
# 
# This software is distributed under the terms of the GNU General Public License
#

import numpy as np
from numpy import arange, sin, cos, pi

ax2 = None
xs = None
ys = None
line = None
fig = None
A = None
P = None
L = None
N = None
MU = None

class PointBrowser:
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """
    def __init__(self, ax, axis2, xs0, ys0, l0, A0,P0,L0,N0,MU0,fig0):
        global ax2, xs, ys, line,fig,A,P,L,N,MU
        ax2 = axis2
        xs = xs0
        ys = ys0
        line = l0
        fig = fig0
        A = A0
        P = P0
        L = L0
        N = N0
        MU = MU0

        self.lastind = 0

        self.text = ax.text(0.05, 0.95, 'selected: none',
                            transform=ax.transAxes, va='top')
        self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.75,
                                  color='yellow', visible=False)

    def onpress(self, event):
        if self.lastind is None: return
        if event.key not in ('n', 'p'): return
        if event.key=='n': inc = 1
        else:  inc = -1


        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
        self.update()

    def onpick(self, event):

        print 'hello'
        if event.artist!=line: return True

        n = len(event.ind)
        if not n: return True
        
        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        
        distances = np.hypot(x-xs[event.ind[0]], y-ys[event.ind[0]])
        indmin = distances.argmin()
        dataind = event.ind[indmin]

        self.lastind = dataind
        self.update()

    def update(self):
        global ax2
        if self.lastind is None: return

        dataind = self.lastind

        ax2.cla()
        div = 2*N+1
        x = arange(0,P*L,.01)
        y = np.zeros(len(x),dtype=np.complex128)
        columnindex = dataind/(2*N+1)
        rowindex = dataind%(2*N+1)
        
        for k in range(0,len(x)):
            m = 0

            for j in range(-N,N+1):
                f_exp = x[k]*(MU[columnindex]+pi*j/L)

                y[k] += A[2][columnindex][rowindex][m] * (cos(f_exp)+complex(0,1)*sin(f_exp))                
                m += 1


        ax2.plot(x,y)

        ax2.text(0.05, 0.9, 'mu=%1.3f\nevalue=%1.3f + %1.3fi' %(MU[columnindex], A[0][columnindex][rowindex].real, A[0][columnindex][rowindex].imag), transform=ax2.transAxes, va='top')
        ax2.set_xlim(0, P*L)

        self.selected.set_visible(True)
        self.selected.set_data(xs[dataind], ys[dataind])

        self.text.set_text('selected: %d'%dataind)
        fig.canvas.draw()



