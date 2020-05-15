import matplotlib.pyplot as plt;
import numpy as np;
import glob;
import os;
import sys;

if len(sys.argv) == 1:
  print("Usage -- python3 plotData.py <list of filenames>")

if os.path.isdir("./plots_data")==False:
  os.mkdir("./plots_data");

# Get all the data files
sim_files = [];
for i in range(1, len(sys.argv)):
  sim_files.append(glob.glob(sys.argv[i]));

# We need to flatten the array.
sim_flat = [];
for l in sim_files:
  for item in l:
    sim_flat.append(item);

n = len(sim_flat);

# Find the maximum and minimum values across all plots.
max_y = 0;
min_y = 0;
for i in range(0,n):
  data = np.genfromtxt(fname=sim_flat[i], delimiter=",");
  ydata = data[:,2];
  tmax = np.amax(ydata);
  tmin = np.amin(ydata);
  if tmax > max_y:
    max_y = tmax;
  if tmin > min_y:
    min_y = tmin;

# Now that we have a min and a max range, let's plot our data.

for i in range(0,n):
  data = np.genfromtxt(fname=sim_flat[i], delimiter=",");

  # Get the head of the file name.
  name = os.path.split(sim_flat[i])[1];
  name = os.path.splitext(name)[0];

  # Make the plots
  plt.plot(data[:,1], data[:,2],'.');
  plt.ylim(min_y,max_y);
  plt.xlabel('x');
  plt.ylabel(name);
  plt.tight_layout();

  plt.savefig('./plots_data/' + name + '.png');
  plt.clf();
