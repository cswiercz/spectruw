#
# Copyright 2013 Thomas Trogdon
#                Krithika Manohar
# This software is distributed under the terms of the GNU General Public License
#


from Tkinter import *
from EntryPanel import *
from SessionPanel import *

import sys
import csv
import pdb
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

if __name__ == '__main__':

    mainPanel = Panel()
    mainPanel.wm_title("SpectruW")

    def _import():
        filename = filebox.get()
        reader = csv.DictReader(open(filename,'rb'),delimiter=',')
        dict =  {row['key']:row['value'] for row in reader}
        for k,v in dict.iteritems():
            try:
                dict[k] = float(v) if '.' in v else int(v)
            except ValueError:
                dict[k] = v
        SessionPanel(dict)
                
    Button(master=mainPanel, text='New', command=lambda:EntryPanel()).grid(row=0)
    Button(master=mainPanel, text='Import',command=_import).grid(row=1)
    filebox = Entry(master=mainPanel)
    filebox.grid(row=1,column=1)
    Tk.mainloop()
