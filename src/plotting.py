"""
Functionality for plotting of individual snapshots of curves.

"""
import numpy as np
import matplotlib.pyplot as plt

class LineImage:
    def __init__(self, xmin=-2.5, xmax=2.5, ymin=-2.5, ymax=2.5, 
                 figsize=(3,3), linewidth=4, show_axes=False,
                 fg_color='black', bg_color=None):

        self.dimensions = [xmin, xmax, ymin, ymax]
        self.figsize = figsize
        # self.bg_color = # TODO: Implement this
        self.fg_color = fg_color
        self.linewidth = linewidth
        self.show_axes = show_axes

        self.fig = plt.figure(figsize=self.figsize)
        self.ax  = self.fig.add_subplot(111)

    def plot_curve(self, c):
        self.ax.plot(c.points[:, 0], c.points[:, 1], color=self.fg_color,
                     linewidth=self.linewidth, solid_capstyle='round')

    def finalize(self):

        if not self.show_axes:
            self.ax.set_frame_on(False)
            self.ax.axes.get_yaxis().set_visible(False)
            self.ax.axes.get_xaxis().set_visible(False)

        self.ax.axis('equal')    
        self.ax.axis(self.dimensions)

    def save(self, f): # TODO kwargs
        self.finalize()
        plt.savefig(f)

    def show(self):  # TODO kwargs
        self.finalize()
        plt.show()

    def close(self):
        plt.close()
