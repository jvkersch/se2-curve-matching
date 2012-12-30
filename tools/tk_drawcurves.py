from Tkinter import *
import tkFileDialog

import scipy.io as sio

# TODO store curves in normalized form, i.e. in RH cartesian frame with
# coords between 0 and 1

def tk_coords(points):
        """Returns coordinates in a flat list, understood by Tk."""
        return [c for pt in points for c in pt]

class Curves:
    """Store points on a number of curves in 2D Cartesian coordinate system."""

    def __init__(self):
        self.curves = [[]]
        self.n = 0

    def num_points(self):
        """Returns number of points."""
        return self.n

    def new_curve(self):
        """Adds new curve."""
        self.curves.append([])

    def del_curve(self, n=-1):
        """Delete selected curve."""
        try: 
            c = self.curves[n]
            self.n -= len(c)
            del self.curves[n]
        except IndexError:
            pass

        if len(self.curves) == 0:
            self.curves.append([])

    def add_point(self, coords):
        """Adds point."""
        x, y = coords
        self.curves[-1].append((x, y))
        self.n += 1

    def clear(self):
        """Removes all points."""
        self.curves = [[]]
        self.n = 0

    def save(self, filename):
        """Saves curves to matlab file."""
        d = {}
        for n, curve in enumerate(self.curves):
            if len(curve) > 0:
                d['c%d' % n] = curve
        sio.savemat(filename, d, oned_as='row')



class StatusBar(Frame):
    """
    Status bar class, taken from 

    http://www.pythonware.com/library/tkinter/introduction/x996-status-bars.htm

    """
    def __init__(self, master):
        Frame.__init__(self, master)
        self.label = Label(self, bd=1, relief=SUNKEN, anchor=W)
        self.label.pack(fill=X)

    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()

    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()


class App: 
    WIDTH  = 300
    HEIGHT = 300

    def __init__(self, master, curves):

        self.curves = curves

        frame = Frame(master)
        frame.pack()

        self.canvas = Canvas(frame, width=self.WIDTH, height=self.HEIGHT, 
                             background='#FF9')
        self.canvas.configure(cursor="crosshair")
        self.canvas.bind("<B1-Motion>", self.on_mouse_down)
        self.canvas.pack(side=TOP)

        self.new_curve = Button(frame, text="New Curve", command=self.new_curve)
        self.new_curve.pack(side=LEFT)

        self.del_curve = Button(frame, text="Delete Last Curve", 
                                command=self.del_curve)
        self.del_curve.pack(side=LEFT)

        self.clear = Button(frame, text="Clear", command=self.clear_drawing)
        self.clear.pack(side=LEFT)

        self.save = Button(frame, text="Save", command=self.save_drawing)
        self.save.pack(side=LEFT)

        self.status = StatusBar(master)
        self.status.pack(side=BOTTOM, fill=X)

    def new_curve(self):
        self.curves.new_curve()
        self.redraw()

    def del_curve(self):
        self.curves.del_curve()
        self.redraw()

    def clear_drawing(self):
        """Discard current drawing."""
        self.curves.clear()
        self.redraw()

    def redraw(self):
        """Redraw current curves and set status."""
        self.canvas.delete(ALL)

        # Draw
        for curve in self.curves.curves:
            if len(curve) > 1:
                pts = tk_coords(curve)
                self.canvas.create_line(pts)

        # Update status 
        self.status.set("Points: %d", self.curves.num_points())

    def save_drawing(self):
        
        filename = tkFileDialog.asksaveasfilename(defaultextension=".mat")
        self.curves.save(filename)

    def on_mouse_down(self, event):
        """Record location of cursor if inside canvas."""

        x, y = event.x, event.y
        if x < 0 or x > self.WIDTH or y < 0 or y > self.HEIGHT:
            return 

        self.curves.add_point((event.x, event.y))
        self.redraw()
        
root = Tk()
app = App(root, Curves())
root.mainloop()
