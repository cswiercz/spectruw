#
# Copyright 2013 Thomas Trogdon
#                Krithika Manohar
# This software is distributed under the terms of the GNU General Public License
#

import sys
from Tkinter import *
from EntryPanel import *
from SessionPanel import *

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

if __name__ == '__main__':

    mainPanel = Panel()
    mainPanel.wm_title("SpectruW")

    Button(master=mainPanel, text='New', command=lambda:EntryPanel()).pack()
    Button(master=mainPanel, text='Import').pack()
    Tk.mainloop()
