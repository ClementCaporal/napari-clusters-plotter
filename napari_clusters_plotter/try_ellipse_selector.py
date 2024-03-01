
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import EllipseSelector, RectangleSelector
import matplotlib.patches as patches
from matplotlib.backend_bases import MouseButton

from nap_plot_tools import make_cat10_mod_cmap
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.path import Path


class AspectRatioEllipse(patches.Ellipse):
    """Class to handle fixed ellipse creation and properties with a crosshair."""
    
    def __init__(self, ax, xy, width, crosshair_color, height2width_ratio = 0.5, **kwargs):
        # Calculate height as a fixed proportion of the width
        self.height2width_ratio = height2width_ratio
        height = self.height2width_ratio * width
        # Initialize the parent class
        super().__init__(xy, width, height, **kwargs)
        self.ax = ax
        # Add a crosshair at the center
        self.center_crosshair = ax.scatter(*xy, color=crosshair_color, marker='+')
        ax.add_patch(self)
        # Redraw the canvas to reflect the changes
        ax.figure.canvas.draw_idle()

    @property
    def center(self):
        return self.get_center()
    
    @property
    def width(self):
        return self.get_width()

    def contains_event(self, event):
        return self.contains(event)[0]
    
    def remove(self):
        self.remove()
        self.center_crosshair.remove()
        self.ax.figure.canvas.draw_idle()
    
    def set_center(self, center):
        super().set_center(center)
        self.center_crosshair.set_offsets(center)
        self.ax.figure.canvas.draw_idle()
    
    def set_size(self, width):
        self.set_width(width)
        self.set_height(self.height2width_ratio*width)
        self.ax.figure.canvas.draw_idle()

    def set_edge_style(self, linestyle, linewidth):
        self.set_linestyle(linestyle)
        self.set_linewidth(linewidth)
        self.ax.figure.canvas.draw_idle()

    def set_edge_color(self, color):
        self.set_edgecolor(color)
        self.center_crosshair.set_color(color)
        self.ax.figure.canvas.draw_idle()
    
class CustomEllipseSelector:
    """Class to handle multi ellipse selections

    Returns
    -------
    _type_
        _description_
    """    
    SCROLL_STEP = 0.05
    MIN_RADIUS = 0.05
    def __init__(self, ax, full_data, parent=None):
        self.plotter = parent
        self.colormap = make_cat10_mod_cmap()
        fig.canvas.mpl_connect('button_press_event', self.on_press)
        fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        fig.canvas.mpl_connect('button_release_event', self.on_release)
        fig.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.axes = ax
        self.ellipses = []  # List to store all ellipses
        self.active_ellipse = None
        self.pressevent = None
        
        self.full_data = full_data


    def get_active_ellipse(self, event):
        """Get the ellipse that contains the event, if any

        Parameters
        ----------
        event : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """        
        for ellipse in self.ellipses:
            if ellipse.contains_event(event):
                return ellipse
        return None

    def on_press(self, event, color='yellow'):
        """Create a ellipse with left click or move an existing ellipse with right click

        Parameters
        ----------
        event : _type_
            _description_
        """        
        if event.inaxes != self.axes:
            return
        if event.button == MouseButton.RIGHT:
            self.clear_all_ellipse_selections()
        elif event.button == MouseButton.LEFT:
            self.clear_all_ellipse_selections()
            self.active_ellipse = self.get_active_ellipse(event)
            if self.active_ellipse is not None:
                self.pressevent = event
                # Store the initial press coordinates and the ellipse's initial center
                self.initial_press_xdata = event.xdata
                self.initial_press_ydata = event.ydata
                self.initial_ellipse_center = self.active_ellipse.center
                # Add dashed contour to active ellipse
                self.active_ellipse.set_edge_style('dashed', 3)
                self.axes.figure.canvas.draw_idle()
            else:
                new_ellipse = AspectRatioEllipse(ax=self.axes, 
                                    xy=(event.xdata, event.ydata), 
                                    width=0.25, 
                                    crosshair_color=self.colormap(self.plotter.color_idx),
                                    alpha=0.5,
                                    fill=False,
                                    ec=self.colormap(self.plotter.color_idx)
                                    )
                self.ellipses.append(new_ellipse)
                self.active_ellipse = self.get_active_ellipse(event)
                self.active_ellipse.set_edge_style('dashed', 3)
                self.update_hue_inside_ellipse()
            

    def clear_all_ellipse_selections(self):
        """Clear all ellipse selections"""   
        for ellipse in self.ellipses:
            ellipse.set_edge_style('solid', 2)
        self.axes.figure.canvas.draw_idle()

    def on_scroll(self, event):
        """Increase or decrease ellipse radius_x with scroll wheel

        Parameters
        ----------
        event : _type_
            _description_
        """
        if event.inaxes != self.axes:
            return
        if self.active_ellipse is None:
            return
        if self.active_ellipse.contains_event(event):
            increment = self.SCROLL_STEP if event.button == 'up' else -self.SCROLL_STEP
            self.active_ellipse.set_size(self.active_ellipse.width + increment)
            self.update_hue_inside_ellipse()
            if self.active_ellipse.width < self.MIN_RADIUS:
                self.active_ellipse.remove()
                self.ellipses.remove(self.active_ellipse)
                self.active_ellipse = None

    def on_move(self, event):
        """Move ellipse with right click

        Parameters
        ----------
        event : _type_
            _description_
        """
        if not hasattr(self, 'pressevent'):
            return
        # If there is no press event (ellipse was not clicked) or the mouse is not in the press event axes, return
        if self.pressevent is None or event.inaxes != self.pressevent.inaxes:
            return
        # If there is no ellipse drawn or the mouse is not in the ellipse, return
        if self.active_ellipse is None:
            return
        dx = event.xdata - self.initial_press_xdata
        dy = event.ydata - self.initial_press_ydata
        new_center = (self.initial_ellipse_center[0] + dx, self.initial_ellipse_center[1] + dy)
        self.active_ellipse.set_center(new_center)
        self.update_hue_inside_ellipse()

    def update_hue_inside_ellipse(self):
        path = self.active_ellipse.get_path()
        affine_transf = self.active_ellipse.get_patch_transform()
        # Apply Affine 2D to transform path from axes coordinates to data coordinates
        path = affine_transf.transform_path(path)
        # Optionally plot path
        # patch = patches.PathPatch(path, edgecolor='red', lw=2, facecolor='none')
        # self.axes.add_patch(patch)

        self.ind_mask = path.contains_points(self.full_data)
        # Get coordinates of points inside the ellipse and plot scatter with different color
        # Find scatter plot in axes
        scatter = self.axes.collections[0]
        scatter_array = scatter.get_array()
        # Set color of points inside the ellipse
        new_scatter_array = self.ind_mask.astype(int)*self.plotter.color_idx
        # Set new positions that became encompassed by ellipse to the current color index
        scatter_array[new_scatter_array == self.plotter.color_idx] = self.plotter.color_idx
        # Then, set to 0 positions that were previously inside the ellipse and now are not
        mask = (scatter_array == self.plotter.color_idx) & (new_scatter_array == 0)
        scatter_array[mask] = 0
        
        # Update scatter colors
        scatter.set_array(scatter_array)
        self.axes.figure.canvas.draw_idle()



    def on_release(self, event):
        """Reset press event

        Parameters
        ----------
        event : _type_
            _description_
        """        
        if self.pressevent is not None:
            self.pressevent = None

class ColorNumber:
    def __init__(self):
        self.color_idx = 0
    def new_color_idx(self, value):
        self.color_idx = value
        return self.color_idx
            
fig = plt.figure(layout='constrained')
axs = fig.subplots(1)

N = 1000 
x = np.linspace(0, 1, N)
y = np.random.gamma(2, size=N) + x


from matplotlib.colors import ListedColormap, Normalize
my_cmap = make_cat10_mod_cmap(first_color_transparent=False)

normalizer = Normalize(vmin=0, vmax=my_cmap.N - 1)

sct = axs.scatter(x, y, s=10, cmap=my_cmap, c=np.zeros(N), norm=normalizer)

# Class to temporarily mimic the plotter object containing the color spinbox
color_number = ColorNumber()

axs.axis([0, 1, 0, 0.5])
handler = CustomEllipseSelector(axs, full_data=np.array([x, y]).T, parent=color_number)
color_number.new_color_idx(2)
plt.show()