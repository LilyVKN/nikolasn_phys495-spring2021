import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

# create the program description
name = "Compile FORTRAN Cache"
desc = "Reads .fcache files which store particle position values and compiles \
        them into a *.mp4 video"
parser = argparse.ArgumentParser(prog=name,description=desc)

# create the argument list
parser.add_argument('--files','-f',required=True,
                    help='the formatted file template for the cache data of \
                        each frame (e.g. \"frame_%%04d.png\")')
parser.add_argument('--output','-o',default="video",
                    help='name of the output <output>.mp4 (default: video)')

# parse the arguments
args = parser.parse_args()
file_template = args.files
output = args.output

def read_FTCache(filename):
    """
    Reads the given .fcache file

    Parameters
    ----------
    filename : string
        The name of the .fcache file (the file extension is not checked)
    
    Returns
    -------
    x : array_like
        The collated x positions of each particle
    y : array_like
        The collated y positions of each particle
    z : array_like
        The collated z positions of each particle

    Notes
    -----
    The .fache file format must be one line per particle with x, y, z values
    entered using white space as a delimiter. This is done in FORTRAN using
    statements such as: `WRITE(*,*) pos(i,:)` or `WRITE(*,*) x, y, z`

    It might be the case that sufficiently large simulation sets will have too
    many lines of data for python to handle.
    """
    # Attempt to open the file as read-only
    with open(filename,"r") as f:
        # Read all the lines of data
        data = f.readlines()

        # Determine the number of lines in the file (i.e. the particle count)
        pnts = int(len(data))
        x = np.zeros((pnts))    # Create an empty float array for x values
        y = np.zeros((pnts))    # ~ y values
        z = np.zeros((pnts))    # ~ z values

        # Iterate through each line
        for i in range(0,pnts):
            # Retrieve the line as an array of delimited floats
            line = data[i].split()

            # Make sure there are only three values to match the format
            if len(line) < 3:
                print("Error: mismatch in position dimensions")
                print("\tLine {} of file {}",line,filename)
                quit()

            # Retrieve each of the particle's values
            x[i] = float(line[0])
            y[i] = float(line[1])
            z[i] = float(line[2])

    # Return the collated position values
    return x, y, z

# Make a temporary directory for each image frame
os.makedirs("./tmp/python",exist_ok=True)

# Initialize the figure and axes
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.set_xlim(-2., 2.)                        # View Limits
ax.set_ylim(-2., 2.)                        # TODO: Make it adaptable to domain
ax.set_zlim(-2., 2.)
ax.grid()
scatter = ax.scatter3D([], [], [], s=2)

# Start reading each frame and set up the data for plotting
frame = 1
filename = file_template % frame    # Insert the frame number into the template

# Keep reading the frame-organized files in order
while(os.path.exists(filename)):
    # Read the positions from the file and convert it to a graph
    x, y, z = read_FTCache(filename)
    scatter._offsets3d = (x, y, z)
    plt.savefig('./tmp/python/ft2py_{:04d}.png'.format(frame))

    # increment the file for the frame and update the filename
    frame += 1
    filename = file_template % frame

# Compile to mp4 using ffmpeg
os.system("ffmpeg -r 24 -f image2 -i ./tmp/python/ft2py_%04d.png -vcodec \
    libx264 -crf 25 -pix_fmt yuv420p {}.mp4".format(output))

# Remove all temporary files and the folder
shutil.rmtree("tmp/python")