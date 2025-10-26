import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from matplotlib.animation import FuncAnimation

# Read data
def read_data(filename):
    data = np.loadtxt(filename)
    return data[:,0], data[:,1], data[:,2]  # x, y, z

# Initialize plot
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surface = None

# Update function for animation
def update(frame):
    global surface
    filename = f'sol.{frame}.dat'  # Construct filename based on frame index
    x, y, z = read_data(filename)
    
    # Create grid values
    xi = np.linspace(min(x), max(x), num=10)
    yi = np.linspace(min(y), max(y), num=10)
    xi, yi = np.meshgrid(xi, yi)
    
    # Interpolate
    zi = griddata((x, y), z, (xi, yi), method='cubic')

    # Clear the previous surface
    if surface:
        surface.remove()
    
    # Plot new surface
    surface = ax.plot_surface(xi, yi, zi, cmap='viridis')

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Quantity')
    ax.set_title(f'Time = {frame}')

# Animation
num_frames = 250  # Change this based on the number of files/time steps
ani = FuncAnimation(fig, update, frames=num_frames, repeat=False)

plt.show()
