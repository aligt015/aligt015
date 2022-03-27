import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Instead of opening multiple files, it is better to consolidate the data into fewer files if possible.
plt.rcParams.update({'font.size': 22})
diff_coord = pd.read_fwf('diffusedata.txt')
det_coord = pd.read_fwf('detections.txt')
vel = pd.read_fwf('velocity.txt')
ang = pd.read_fwf('vlsr.txt')
diff_long = np.array(diff_coord[diff_coord.columns[0]])
diff_lat = np.array(diff_coord[diff_coord.columns[1]])
det_long = np.array(det_coord[det_coord.columns[0]])
det_lat = np.array(det_coord[det_coord.columns[1]])
velocities = np.array(vel[vel.columns[0]])
lsr = np.array(ang[ang.columns[0]])
angles = np.arcsin(lsr/15.)*100

def plot_detections(point, angle, length):
     '''
     point - array of coordinates based off the detections
     angle - angle from the velocity lsr
     length - linewidths from detections
     '''
     x, y = point

     # formulate the angle based off the lsr velocities
     begx = x - length/2 * np.cos(np.radians(angle))
     begy = y - length/2 * np.sin(np.radians(angle))
     endx = x + length/2 * np.cos(np.radians(angle))
     endy = y + length/2 * np.sin(np.radians(angle))

     # plot the detections
     fig = plt.figure(figsize=(75,55)) # very large figure size in order to easily see the detection points or else it will be cramped and hard to see.
     ax = plt.subplot(111)
     ax.set_xlim(360., 0.)
     ax.set_ylim(-90., 90.)
     ax.set_xlabel("Galactic Longitude")
     ax.set_ylabel("Galactic Latitude")
     ax.scatter(diff_long, diff_lat, s=3, color='black')
     ax.scatter(det_long, det_lat, s=20, color='blue')
     ax.plot([begx, endx], [begy, endy], color='blue')
     fig.show()

plot_detections([det_long, det_lat], angles, velocities)
