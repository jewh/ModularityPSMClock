# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# This function takes a .csv file as input (corresponding to the saved output of our simulations), and renders the output into a 3D animation with cells coloured according to their clock phase.

import matplotlib
matplotlib.use('Agg') # Activate 'agg' backend for off-screen plotting.
import ovito 
import seaborn
import numpy
import os
from tracks_to_xyz import convert_file
import matplotlib.pyplot as plt
import PySide2.QtGui
import time
import PySide2.QtCore

# set global font
# 
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = "Arial"

def find_nearest(array, value):
    """
    Finds the index of the entry in an array with the closest value to input 'value'.

    Input: array (array-like) - n x 1 array of ints or floats.
    Output: index (int). 
    """
    array = numpy.asarray(array)
    index = (numpy.abs(array - value)).argmin()
    return index

def colour_scheme(rgb_array):
    """
    Computes the new colour for a colour value.
    Uses a defined colourmap ('c_map') to compare values against.

    Note that values in rgb_array must be in [0,1].

    Input: rgb_array (tuple) - 3-tuple of floats with values 0<=value<=1.
    Output: 3-tuple of floats with values 0<=value<=1. 
    """   
    red = float(rgb_array[0])
    c_map = seaborn.color_palette("flare", as_cmap=True)
    # c_map = seaborn.light_palette("seagreen", as_cmap=True)
    # c_map = seaborn.cubehelix_palette(start=-.5, rot=-.75, dark=0, light=0.95, as_cmap=False)
    return c_map(red)[0:3]   

def compute_colour(frame, data, colour_scheme=colour_scheme):
    """
    Compute the colour for particles in the pipeline.
    Uses the syntax for defining an Ovito modifier.

    Input: 
    frame (?) - needed syntactically. 
    data (ovito data object) - modified by the present function.
    colour_scheme (function) - tuple-valued function calculating 
    the colour according to a given scheme.

    Output: None
    """
    raw_colour_values = [list(color_object) for color_object in data.particles['Color']]
    processed_colour_values = [colour_scheme(color_object) for color_object in raw_colour_values]
    # assign the particle colors
    data.particles_.create_property('Color', data=processed_colour_values)
    
def somite(t):
    return 6.0 + (t / 24.7) - 0.5001 * numpy.exp(0.0049 * (t - 300.0))

def cell_diameter(t):
    return 11.0
    if t < 253.6:
        return 9.2
    elif t >= 253.6:
        return 9.2 - 0.2 * (somite(t) - 16)

def set_cell_size(frame, data, radius=5.5):
    """Ovito modifier to set the particle size. Sets all particles to the same size."""
    t = data.attributes['Time']
    radius = cell_diameter(t) / 2.0
    # new_radii = [radius for xyz in data.particles['Position']]
    new_radii = []
    for xyz in data.particles['Position']:
        if xyz[0] < 25.0 + xa(t):
            # tl = (25.0 + xa(t) - xyz[0]) / 1.67
            # new_radii.append(cell_diameter(t - tl))
            new_radii.append(radius)
        else:
            new_radii.append(radius)
    data.particles_.create_property('Radius', data=new_radii)

def reflect_x_axis(frame, data):
    """Ovito modifier that reflects the cell positions in the x-axis."""
    positions = [(-x,y,z) for (x,y,z) in data.particles_['Position']]
    data.particles_.create_property('Position', data=positions)

def reflect_y_axis(frame, data):
    """Ovito modifier that reflects the cell positions in the y-axis."""
    positions = [(x,-y,z) for (x,y,z) in data.particles_['Position']]
    data.particles_.create_property('Position', data=positions)

def reflect_z_axis(frame, data):
    """Ovito modifier that reflects the cell positions in the z-axis."""
    positions = [(x,y,-z) for (x,y,z) in data.particles_['Position']]
    data.particles_.create_property('Position', data=positions)

def set_origin(frame, data):
    """Ovito modifier that sets origin of data to (0, 0, 0)."""
    positions = [(x-25,y-25,z-25) for (x,y,z) in data.particles_['Position']]
    data.particles_.create_property('Position', data=positions)

def set_transparency(frame, data, transparency=0.0):
    """
    Ovito modifier to set the particle transparency. Sets all particles to the same transparency.
    """
    trans = [transparency for xyz in data.particles['Position']]
    data.particles_.create_property('Transparency', data=trans)

def find_between(string, start, end):
    return (string.split(start))[1].split(end)[0]

def get_frames(directory):
    """
    Returns a list of files which are 'frames' for the animation, in chronological order.
    """
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    files.sort(key=lambda f: float(find_between(f, "snapshot_t", ".xyz")))
    files = [os.path.join(directory, f) for f in files]
    return files

def plot_colorbar(args):
    """
    Plot a colour bar on the figure.
    """

    # define size of plot
    dpi = 160
    # needs to be 0.4 in both axes to include colorbar label
    plot_width, plot_height = 0.4 * args.size[0] / dpi, 0.4 * args.size[1] / dpi

    # create figure
    fig, ax = plt.subplots(figsize=(plot_width, plot_height), dpi=dpi)
    fig.set_facecolor((0.99, 0.99, 0.99)) # so that we don't get a yellow halo
    # hide axes
    fig.gca().set_visible(False)
    # make colorbar
    c_map = seaborn.color_palette("flare", as_cmap=True)
    # c_map = seaborn.light_palette("seagreen", as_cmap=True)
    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=c_map))
    # set properties
    cb.set_label(label=r'$\frac{1 + \sin \theta_{i}}{2}$', weight='bold', size=20)
    tick_font_size = 10
    cb.ax.tick_params(labelsize=tick_font_size)

    # render figure to buffer
    buf = fig.canvas.print_to_buffer()
    plt.close(fig)

    # Create a QImage from the memory buffer
    res_x, res_y = buf[1]
    img = PySide2.QtGui.QImage(buf[0], res_x, res_y, PySide2.QtGui.QImage.Format_RGBA8888)

    # Paint QImage onto viewport canvas
    args.painter.drawImage(0.5 * args.size[0], 0.5 * args.size[1] - 0.5 * plot_height * dpi, img)
    
def somite(t):
    return 6. + (t / 24.7) - 0.5001 * numpy.exp(0.0049 * (t - 300.))

def xa(t):
    # return 0
    if t < 253.6:
        return 0.0
    else:
        return 14.1 * (somite(t) - somite(253.6))

def plot_xa(args):
    """
    Plot a line showing xa on the figure.
    """
    bar_length = 225  # Simulation units (e.g. Angstroms)
    # bar_length = 75 # for side-on rendering
    bar_color = PySide2.QtGui.QColor(100, 100, 100)

    if args.is_perspective:
        raise Exception("This overlay only works with non-perspective viewports.")

    # Compute length of bar in screen space

    screen_length = args.project_size((0,0,0), bar_length)

    data = args.scene.pipelines[0].compute(args.frame)
    t = data.attributes['Time']

    # Define geometry of bar in screen space
    width = 0.002 * args.painter.window().width()
    xmargin, ymargin = args.project_point([xa(t), -25, -25])
    ymargin -= screen_length # silence for side-on
    rect = PySide2.QtCore.QRectF(xmargin, ymargin, width, screen_length)

    # Render bar rectangle
    args.painter.fillRect(rect, bar_color)

    # plot the label
    font = args.painter.font()
    font.setPixelSize(40)
    args.painter.setFont(font)
    args.painter.setPen(PySide2.QtGui.QPen(bar_color))
    args.painter.drawText(xmargin - 45, ymargin, "x = xₐ")

def plot_x0(args):
    """
    Plot a black line showing xa = 0 on the figure.
    """
    bar_length = 225  # Simulation units (e.g. Angstroms)
    bar_color = PySide2.QtGui.QColor(0, 0, 0)

    if args.is_perspective:
        raise Exception("This overlay only works with non-perspective viewports.")

    # Compute length of bar in screen space

    screen_length = args.project_size((0,0,0), bar_length)

    data = args.scene.pipelines[0].compute(args.frame)
    t = data.attributes['Time']

    # Define geometry of bar in screen space
    width = 0.002 * args.painter.window().width()
    xmargin, ymargin = args.project_point([0, -25, 0])
    ymargin -= screen_length - 10
    rect = PySide2.QtCore.QRectF(xmargin, ymargin, width, screen_length)

    # Render bar rectangle
    args.painter.fillRect(rect, bar_color)

def plot_xd(args):
    """
    Plot a line showing xd on the figure.
    """
    bar_length = 225  # Simulation units (e.g. Angstroms)
    bar_color = PySide2.QtGui.QColor(100, 100, 100)

    if args.is_perspective:
        raise Exception("This overlay only works with non-perspective viewports.")

    # Compute length of bar in screen space

    screen_length = args.project_size((0,0,0), bar_length)

    data = args.scene.pipelines[0].compute(args.frame)
    t = data.attributes['Time']

    # Define geometry of bar in screen space
    width = 0.002 * args.painter.window().width()
    xmargin, ymargin = args.project_point([100.0, -25, 0])
    ymargin -= screen_length - 10
    rect = PySide2.QtCore.QRectF(xmargin, ymargin, width, screen_length)

    # Render bar rectangle
    args.painter.fillRect(rect, bar_color)

    # plot the label
    font = args.painter.font()
    font.setPixelSize(40)
    args.painter.setFont(font)
    pen = PySide2.QtGui.QPen(PySide2.QtCore.Qt.DashLine)
    pen.setColor(bar_color)
    args.painter.setPen(pen)
    args.painter.drawText(xmargin - 45, ymargin, "x = xd")
    
def plot_title(args):
    """
    Plot a title for the figure.
    """

    # define size of plot
    dpi = 160
    # needs to be 0.4 in both axes to include colorbar label
    plot_width, plot_height = 0.5 * args.size[0] / dpi, 0.3 * args.size[1] / dpi

    # create figure
    fig, ax = plt.subplots(figsize=(plot_width, plot_height), dpi=dpi)
    fig.set_facecolor((0.99, 0.99, 0.99)) # so that we don't get a yellow halo
    # hide axes
    fig.gca().set_visible(False)

    fig.suptitle(f"Compaction-extension", fontsize=20, fontweight='bold')

    # render figure to buffer
    buf = fig.canvas.print_to_buffer()
    plt.close(fig)

    # Create a QImage from the memory buffer
    res_x, res_y = buf[1]
    img = PySide2.QtGui.QImage(buf[0], res_x, res_y, PySide2.QtGui.QImage.Format_RGBA8888)

    # Paint QImage onto viewport canvas
    args.painter.drawImage(0.5 * args.size[0] - 0.5 * plot_width * dpi, 0.5 * args.size[1] - plot_height * dpi, img)

def assign_time(frame, data):
    """
    Infer from the .xyz file what the time of the frame is.
    """
    _comment = data.attributes['Comment']
    _to_parse = find_between(_comment, '= ', '.')
    data.attributes['Time'] = float(_to_parse)

# This function is called by OVITO on every viewport update.
def render_scalebar(args):
    bar_length = 100   # Simulation units (e.g. Angstroms)
    bar_color = PySide2.QtGui.QColor(0,0,0)
    label_text = f"{bar_length} μm"
    label_color = PySide2.QtGui.QColor(0, 0, 0)
    if args.is_perspective:
        raise Exception("This overlay only works with non-perspective viewports.")

    # Compute length of bar in screen space
    screen_length = args.project_size((0,0,0), bar_length)

    # Define geometry of bar in screen space
    height = 0.01 * args.painter.window().height()
    xmargin = args.painter.window().width() - screen_length - 0.15 * args.painter.window().width()
    ymargin = 0.82 * args.painter.window().height()
    rect = PySide2.QtCore.QRectF(xmargin, ymargin, screen_length, height)

    # Render bar rectangle
    args.painter.fillRect(rect, bar_color)

    # Render text lab
    textrect = PySide2.QtCore.QRectF(xmargin, ymargin - height * 3 , screen_length, 3 * height)
    font = args.painter.font()
    font.setPixelSize(30)
    args.painter.setFont(font)
    args.painter.setPen(PySide2.QtGui.QPen(label_color))
    args.painter.drawText(textrect, PySide2.QtCore.Qt.AlignCenter, label_text)



# Generate a directory of .xyz files
data_file = r"data/20231216_frequency_tracking_dynamic_tissue_simulation_constantadvection__constantdensity_constant_dc_shrinking_randomnewcells_seed9300.csv"

t1 = time.time()
# convert_file(data_file, 'Time', 'color', 1)
# convert_file(data_file, 'Time_Mins', 'pher1', 1)
convert_file(data_file, 'Time_Mins', 'Phase_sim', 1)
print(f"\nWriting .xyz files took {time.time() - t1}s\n")

# want to import the files into ovito 
t1 = time.time()
try:
    pipeline = ovito.io.import_file("snapshots/snapshot_t*.*.xyz",
    columns = ["Position.X", "Position.Y", "Position.Z", "Color"])
except RuntimeError: # yields this when no files with the above pattern are found.
    pipeline = ovito.io.import_file("snapshots/snapshot_t*.xyz",
    columns = ["Position.X", "Position.Y", "Position.Z", "Color"])
print(f"Importing files into pipeline took {time.time() - t1}s\n")


t1 = time.time()
# specify the time
pipeline.modifiers.append(assign_time)
# Specify the colour and size of the particles. 
pipeline.modifiers.append(compute_colour)
pipeline.modifiers.append(set_cell_size)
# set origin to zero (may not be necessary)
pipeline.modifiers.append(set_origin)
# Invert the x,y axes (necessary for some tracks)
# pipeline.modifiers.append(reflect_x_axis)
# pipeline.modifiers.append(reflect_y_axis)
# Set the transparency of the particles in the simulation
# pipeline.modifiers.append(set_transparency)

# remove the black frame around the image
cell_vis = pipeline.source.data.cell.vis
cell_vis.line_width = 0

# Now load the data into the renderer
pipeline.add_to_scene()

# zoom_factor = 10
# vp = ovito.vis.Viewport(camera_pos = (-153,-15,36.5), camera_dir = (0, 0, -1)) # M4 right hand side of PSM
vp = ovito.vis.Viewport()
vp.type = ovito.vis.Viewport.Type.Ortho
vp.camera_pos = (-300, 85, 500) # for long figures
# vp.camera_pos = (250, 85, 0)
# vp.camera_pos = (550, 335, 275)
# vp.camera_dir = (1, 1, -1)
vp.camera_dir = (0, 0, -1)
# append to underlays so that the figure's white background doesn't plot over the animation
# vp.underlays.append(ovito.vis.PythonViewportOverlay(function = plot_xa))
vp.underlays.append(ovito.vis.PythonViewportOverlay(function = plot_x0))
# vp.underlays.append(ovito.vis.PythonViewportOverlay(function = plot_xd))
# vp.underlays.append(ovito.vis.PythonViewportOverlay(function = plot_colorbar))
# vp.underlays.append(ovito.vis.PythonViewportOverlay(function = plot_title))
vp.overlays.append(ovito.vis.PythonViewportOverlay(function = render_scalebar))

# plot the time of the frame
text_overlay = ovito.vis.TextLabelOverlay(
    text = 't = [Time] mins',
    # alignment = PySide2.QtCore.Qt.AlignHCenter ^ PySide2.QtCore.Qt.AlignBottom,
    offset_x = 0.45,
    offset_y = -0.8,
    font_size = 0.01,
    text_color = (0,0,0))

# Attach the overlay to viewport:
vp.overlays.append(text_overlay)


# set the background colour of the animation here, using RGB format
background_colour = [0.99,0.99,0.99] # default grey is 0.239, 0.239, 0.239 (higher for lighter)
# Set the resolution here - put it to 1080p for default
resolution = [3000, 520] # this is the standard 1920x1080 pixels (x,y)
# resolution = [1920, 1080]
# Now render the animation and save it as you wish - can be .mp4, .mov, .avi, .gif 
file_name = "animations/20231216_frequency_tracking_dynamic_tissue_simulation_constantadvection__constantdensity_constant_dc_shrinking_randomnewcells_seed9300.mp4"
# Now render the animation - can adjust the fps here if you wish
# TODO choose a best renderer - ovito gives us two options
print(f"Setting animation parameters took {time.time() - t1}s\n")
t1 = time.time()
vp.render_anim(file_name, size=resolution, fps=20, background=background_colour)

# Or, if we wish to only render stills, use this
# for frame in [t + 1 for t in range(444, 470)]:
for frame in [600]:

    filename = f"figures/20231216_frequency_tracking_dynamic_tissue_simulation_constantadvection__constantdensity_constant_dc_shrinking_randomnewcells_seed9300_t{frame}.png"

    vp.render_image(size=resolution, filename=filename, background=background_colour, frame=frame)

print(f"Rendering animation took {time.time() - t1}s\n")
# Works!
