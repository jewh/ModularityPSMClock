# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Converts .csv cell position and phase data to .xyz format so as to plot as an ovito animation.

import numpy as np
import pandas as pd
import os 
from defunct.kuramoto import translated_phase

def to_xyz(array, time):
    """
    Input: Array with four columns. Assumes the first three are cartesian coordinates. 
    Fourth is assumed to be some property of the cells at each coordinate - identity, 
    level of gene expression, etc. 

    Output: saves input array as an .xyz file
    """
    # Get the number of 'atoms'. This is the no. of particles in the .xyz file
    no_atoms = str(len(array))
    # Get the comment for the file (goes on line 2).
    comment = f"Snapshot for time = {time}."
    # initialise the output file
    file_name = f"snapshots/snapshot_t{time}.xyz"
    with open(file_name, "w") as f:
        # write the first two lines
        f.write(no_atoms + "\n" + comment + "\n")
        # TODO allow tracking of cells by identity 
        # perhaps let them take on an identity depending on simulation?
        # Now loop over each set of coords in the array and append to the .xyz file
        for index, row in enumerate(array):
            x, y, z, identity, _ = row

            f.write(str(x) + "\t" + str(y) + "\t" + str(z) + "\t" + str(identity) + "\n")

def clear_dir(path):
    """
    Removes all files in the directory at path. 
    Note - this will not remove directories.

    Input:
    path (str) - path to directory to clear. 

    Output:
    Deletes all files in the directory at {path}.
    """
    if not os.path.isdir(path):
        raise ValueError(f"'{path}' is not a directory or does not exist.")

    for file in os.listdir(path):
        if os.path.isfile(path + "/" + file):
            os.remove(path + "/" + file)

def get_times(data, dt, t_axis):
    """
    Returns list of time values, taken at even intervals from the time course.
    """
    n = int(1 / dt) 
    times = list(set(data[t_axis]))
    return times[::n]

def translate_phase(x):
    return 0.5 * (1.0 + np.sin(x))

def normalised_color(df, q):
    return [x / q for x in df.to_numpy()]



def convert_file(file, t_axis, p_axis, dt):
    """Convert .csv tracks to directory of .xyz files."""
    # clear the directory
    clear_dir("snapshots")
    # Now convert to .xyz
    # Subset by time 
    data = pd.read_csv(file) # TODO alter script to accept .txt or .xlsx files
    # remove NaN values
    data.dropna(subset=["Position X Reference Frame", "Position Y Reference Frame", "Position Z Reference Frame"], inplace=True)
    times = get_times(data, dt, t_axis)

    if p_axis in ["pher1", "pher7", "pdelta", "mher1", "mher7", "mdelta"]:

        # find maximum value and normalise by it
        q = data[p_axis].max()
        data[p_axis] = normalised_color(data[p_axis], q)

    # save a snapshot for each time
    for t in times:
        snapshot = data[data[t_axis] == t] # subset by that time
        if snapshot.empty:
            indx = data[t_axis].sub(t).abs().idxmin() # get next closest value. 
            t_close = float(data[t_axis].iloc[indx])
            snapshot = data[data[t_axis] == t_close]
        snapshot = snapshot[["Position X Reference Frame", 
                        "Position Y Reference Frame",
                        "Position Z Reference Frame", 
                        p_axis,
                        "TrackID"
                        ]] # subset this by relevant columns
        if p_axis == "TrackID":
            # max_val = data[p_axis].max()
            translate_discrete = lambda x : 1.0
            snapshot[p_axis] = snapshot[p_axis].apply(translate_discrete)
        elif p_axis == "color":
            pass
        elif p_axis in ["pher1", "pher7", "pdelta", "mher1", "mher7", "mdelta"]:
            pass
        else:
            snapshot[p_axis] = snapshot[p_axis].apply(translate_phase)
        snapshotnp = snapshot.to_numpy()
        to_xyz(snapshotnp, t)

if __name__ == "__main__":
    data_file = r"data\20221201_frequency_tracking_mitosing_simulation_seed400.csv"
    convert_file(data_file, t_axis="Time_Mins", p_axis="Phase_sim", dt=1)
