import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import EllipseSelector, RectangleSelector
import matplotlib.patches as patches
from matplotlib.backend_bases import MouseButton
from matplotlib.collections import PathCollection

from nap_plot_tools import make_cat10_mod_cmap
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.path import Path

class AspectRatioEllipse(patches.Ellipse):
    """Class to handle fixed ellipse creation and properties with a crosshair."""
    
    def __init__(self, ax, xy, width, color, color_idx, height2width_ratio = 0.5, **kwargs):
        # Store color index
        self.color_idx = color_idx
        # Calculate height as a fixed proportion of the width
        self.height2width_ratio = height2width_ratio
        height = self.height2width_ratio * width
        # Initialize the parent class
        super().__init__(xy, width, height, **kwargs)
        self.ax = ax
        # Creates an array to store the mask of the points inside the ellipse (if scatter plot is present in the axes)
        self._inside_mask_array = None
        # Add a crosshair at the center
        self.center_crosshair = ax.scatter(*xy, color=color, marker='+', label='crosshair')
        ax.add_patch(self)
        # Redraw the canvas to reflect the changes
        ax.figure.canvas.draw_idle()

    @property
    def center(self):
        return self.get_center()
    
    @property
    def width(self):
        return self.get_width()
    
    @property
    def inside_mask_array(self):
        return self._inside_mask_array
    
    @inside_mask_array.setter
    def inside_mask_array(self, values):
        self._inside_mask_array = values

    def contains_event(self, event):
        return self.contains(event)[0]
    
    def remove(self):
        # TODO: fix crosshair removal error
        self.center_crosshair.remove()
        self.remove()
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

    # def set_edge_color(self, color):
    #     super().set_edgecolor(color)
    #     self.center_crosshair.set_color(color)
    #     self.ax.figure.canvas.draw_idle()

    # def find_scatter_plot(self):
    #     for collection in self.ax.collections:
    #         if isinstance(collection, plt.collections.PathCollection) and collection.get_label() != 'crosshair':
    #             return collection
    #     return None
    
    # def create_array_colors(self):
    #     self.array_colors = None
    #     # Find scatter plot in axes (if any) and store it
    #     self.scatter = self.find_scatter_plot()
    #     if self.scatter is not None:
    #         self.array_colors = np.zeros(len(self.scatter.get_array()))
    
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

    def on_press(self, event):
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
                color = self.colormap(self.plotter.color_idx)
                new_ellipse = AspectRatioEllipse(ax=self.axes, 
                                    xy=(event.xdata, event.ydata), 
                                    width=0.25,
                                    color_idx=self.plotter.color_idx,
                                    color=color,
                                    alpha=0.5,
                                    fill=False,
                                    ec=color
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
        """Update the color of the points inside the ellipse"""
        path = self.active_ellipse.get_path()
        affine_transf = self.active_ellipse.get_patch_transform()
        # Apply Affine 2D to transform path from axes coordinates to data coordinates
        path = affine_transf.transform_path(path)
        # Check if points are inside the ellipse
        self.ind_mask = path.contains_points(self.full_data)
        # Find scatter plot in axes (if any) and store it
        scatter = self.find_scatter_plot()
        if scatter is None:
            return
        # Get current scatter color indices array
        scatter_array = scatter.get_array()
        # If active ellipse has no inside_mask_array yet, create it
        if self.active_ellipse.inside_mask_array is None:
            self.active_ellipse.inside_mask_array = np.zeros(len(scatter_array)).astype(bool)
        
        # Store mask of points inside the ellipse
        self.active_ellipse.inside_mask_array = self.ind_mask
        new_scatter_array = self.active_ellipse.inside_mask_array.astype(int)*self.active_ellipse.color_idx
        
        # Set new positions that became encompassed by ellipse to the current color index
        scatter_array[new_scatter_array == self.active_ellipse.color_idx] = self.active_ellipse.color_idx
        # Then, set to 0 positions that were previously inside the ellipse and now are not
        mask = (scatter_array == self.active_ellipse.color_idx) & (new_scatter_array == 0)
        # TODO: Check if these positions are inside other ellipses and preserve their color
        scatter_array[mask] = 0
        
        # Update scatter colors
        scatter.set_array(scatter_array)
        self.axes.figure.canvas.draw_idle()

    def find_scatter_plot(self):
        for collection in self.axes.collections:
            if isinstance(collection, PathCollection) and collection.get_label() != 'crosshair':
                return collection
        return None

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
color_number.new_color_idx(1)
plt.show()

a=1