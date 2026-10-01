import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import icemodel
from functools import partial

plt.close('all')

# Geometrical constants
L = 1.e5  # length of domain (m)
initialthickness=800
angle = 0.2  # bedrock slope in degrees
amplitudes=[-0.5,0.9]

# Numerical constants
m = 51  # number of grid points
x = np.linspace(0, L, m)  # horizontal distance
dx = L / (m - 1)   # dx *(m-1)=L, m-1 =L/dx, m= L/dx+1
dt = 0.1  # time step (years)  # Make sure you don't violate the CFL condition: dt < C/dx^2
time_end = 200  # duration (years)


# Physical constants  MPa - m - a
A = 300  # Ice flow parameter
rho = 9.138e-19  # Ice density
grav = 9.7692e15  # Gravitation constant 
n = 3.  # Glen index

# Bed shape setup
b=icemodel.CreateBedRock(angle,amplitudes,x,L)
h=icemodel.CreateInitialSurface(angle,initialthickness,x,L)


fig,ax2,ax2_uz,line1,line2,line3,line4,line5,time_text=icemodel.InitAnimation(L,b,h,h,x)

# Update function for animation
def update(t):
    global h #OBS, bad coding due to laziness
    dhdx, ux, uz = icemodel.VelocityAndSurfaceGradient(b, h, dx,rho,grav,A,n)
    h = h + dt * (-np.multiply(ux, dhdx) + uz)

    ax2.set_ylim([min(-0.01,np.min(ux)), max(0.01,np.max(ux))])
    ax2_uz.set_ylim([min(-0.01,np.min(uz)), max(0.01,np.max(uz))])  # Adjust this range based on uz values

    # Update the data in the first subplot (elevation)
    line2.set_ydata(h);
    
    # Update the data in the second subplot (u_x)
    line4.set_ydata(ux)
    line5.set_ydata(uz)

    
    time_text.set_text(f't={int(t)} years')

    # if t==time_end:
    #     with open('coordinates.npy', 'wb') as fcoords:
    #         np.save(fcoords,x)

    #     with open('surface.npy', 'wb') as fsurface:
    #         np.save(fsurface,h)

    #     with open('ux.npy', 'wb') as fux:
    #         np.save(fux,ux)

    #     with open('uz.npy', 'wb') as fuz:
    #         np.save(fuz,uz)

    return line2, line4, line5, time_text




# Create the animation
frames = np.arange(0, time_end + 1, dt)  # time steps
ani = FuncAnimation(fig, partial(update), frames=frames, blit=True, interval=10, repeat=False)

# Display the animation
plt.show()

