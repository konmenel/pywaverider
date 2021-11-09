import tkinter as tk
import windowConfig as wcfg


def createWindow():
    root = tk.Tk()
    # root.config(bg=wcfg.DARK_GRAY)
    root.geometry(wcfg.SIZE)
    root.title('Waverider Generation')
    root.configure(bg=wcfg.WHITE)


    frame_top = tk.Frame(root)
    frame_top.pack(fill='x')

    plotBaseButton = tk.Button(frame_top, text='Plot Base View', command=drawToCanvas)
    plotBaseButton.pack(side='left', fill='x', expand=True)

    plotTopButton = tk.Button(frame_top, text='Plot Top View', command=drawToCanvas)
    plotTopButton.pack(side='left', fill='x', expand=True)

    plot3DButton = tk.Button(frame_top, text='Plot 3D View', command=drawToCanvas)
    plot3DButton.pack(side='left', fill='x', expand=True)

    wcfg.window = root
    return


def drawToCanvas():
    try: 
        wcfg.plotCanvas.get_tk_widget().pack_forget()
    except AttributeError: 
        pass

    x = np.linspace(0, 5*np.pi, 100)
    y = np.sin(x)
    fig = plt.Figure(figsize=(5, 4), dpi=100)
    fig.add_subplot(111).plot(x, y)
    canvas = FigureCanvasTkAgg(fig, master=wcfg.window)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    if wcfg.plotToolbar is None:
        toolbar = NavigationToolbar2Tk(canvas, wcfg.window)
        toolbar.update()
        wcfg.plotToolbar = toolbar

    wcfg.plotCanvas = canvas
    


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    createWindow()
    window = wcfg.window

    

    
    window.mainloop()
