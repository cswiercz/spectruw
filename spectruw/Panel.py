import Tkinter as Tk

class Panel(Tk.Tk):     
    def __init__(self):
        Tk.Tk.__init__(self)
        #self.configure(background='whitesmoke')
        self.tk_setPalette(background='whitesmoke')
#               activeBackground='gray', activeforeground='white')

    def close_window(self):
        # children destroyed first
        self.destroy()

# TO DO: make active window, show window, hide window
