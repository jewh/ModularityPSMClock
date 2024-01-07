
# Modularity of the segmentation clock and morphogenesis



## Visualisation - Animation rendering using Ovito 

Converts cell position data into animations, created using Ovito, a piece of software used for visualising and simulating particle dynamics in chemistry [ovito.org](https://www.ovito.org/). Here's an example of an animation made using Ovito:

![Animation made using Ovito](frontpageanimation.gif)

Ovito can be run using Python scripts. We're using this functionality here because it allows easy replication of figures. 


Ovito accepts file formats more typical of computational chemistry (.cif, .xyz, etc.). Cell tracking data (and associated gene expression) are typically stored in .csv or .xlsx files, so we need a way to convert the raw data to a file format readable by Ovito. Here I use .xyz, as it's the most straightforward. 

An .xyz file (at least as read by ovito) assumes the first three columns are the cartesian coordinates x, y, z. Ovito can also read a fourth column, which my code reads out as the colour of the particles in Ovito. We thus use this column to store gene expression values for visualising simulations on the cell track.

[tracks_to_xyz.py](tracks_to_xyz.py) will convert a .csv file with cartesian coordinates for cell nuclei over time into a folder of .xyz files, each representing a snapshot of embryo at a specific time. This folder is then read by Ovito in [to_ovito.py](to_ovito.py) to an animation. You can tweak the aesthetic parameters for Ovito (frame rate, colours, particle size, renderer, etc.) in this script, as well as the file format for the animation. Ovito allows .mp4, .avi, .mov, and .gif formats. Alternatively specifying a .jpeg or .png file will create a series of image files, one per frame. 

## Open Issues:

This is far from a smooth and clear workflow, so for now I will focus on allowing easy adjustment of the aesthetic parameters, and clarifying the python scripts so that users may easily animate their data with minimal coding. 
