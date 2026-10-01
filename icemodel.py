import numpy as np
import matplotlib.pyplot as plt

def bump(amplitude,x,position,L,tan):
    bump=amplitude* np.exp(-(2 * (x - position))**2 / (L / 4.)**2) / (0.7 * tan)
    return bump

def CreateBedRock(angle,amplitudes,x,L):
    tan = np.tan(angle * np.pi / 180)
    bedbase=-tan * x
    position1=0.5*L/2.
    position2=1.3*L/2.

    bump1=bump(amplitudes[0],x,position1,L,tan)
    bump2=bump(amplitudes[1],x,position2,L,tan)
    
    b = bedbase + bump1+bump2
    return b

def CreateInitialSurface(angle,initialthickness,x,L):
    tan = np.tan(angle * np.pi / 180)
    bedbase=-tan * x
    h0 = bedbase + initialthickness  # initial ice thickness
    return h0


def VelocityAndSurfaceGradient(b, h, dx,rho,grav,A,n):
    H = h - b
    m= len(h)
    dhdx = np.zeros(m)
    d2hdx2 = np.zeros(m)
    dHdx = np.zeros(m)
    ux = np.zeros(m)
    uz = np.zeros(m)

    heightdiff=h[0]-h[m-1]  # L*tan

    # left point
    j = 0
    dhdx[j] = (h[1] - (h[m-2] + heightdiff)) / (2 * dx)
    d2hdx2[j] = (h[1] - 2 * h[j] + (h[m-2] + heightdiff)) / (dx**2)
    dHdx[j] = (H[1] - H[m-2]) / (2 * dx)

    # interior points
    for j in np.arange(1, m-1, 1):
        dhdx[j] = (h[j+1] - h[j-1]) / (2 * dx)
        d2hdx2[j] = (h[j+1] - 2 * h[j] + h[j-1]) / (dx**2)
        dHdx[j] = (H[j+1] - H[j-1]) / (2 * dx)

    # right point
    j = m-1
    dhdx[j] = ((h[1] - heightdiff) - h[j-1]) / (2 * dx)
    d2hdx2[j] = ((h[1] - heightdiff) - 2 * h[j] + h[j-1]) / (dx**2)
    dHdx[j] = (H[1] - H[j-1]) / (2 * dx)

    for j in np.arange(0, m, 1):
        ux[j] =-2. * A * (rho * grav)**n * dhdx[j]**n * H[j]**(n+1.) / (n+1.)
        uz[j] = 2. / (n+1) * A * (rho * grav)**n * (3. * dhdx[j]**2 * d2hdx2[j] * (4./5.) * H[j]**5. + 4. * dhdx[j]**3. * H[j]**4 * (dHdx[j] - 0.25 * dhdx[j]))

    return dhdx, ux, uz


def InitAnimation(L,b,h,h0,x):

    plt.rcParams['figure.figsize'] = (10, 8)  # Increase figure size for two subplots
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titleweight'] = 'bold'

    # Set up the figure and axes for subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})  # Two subplots, one below the other
    fig.subplots_adjust(hspace=0.4)  # Adjust the space between subplots

    # First subplot: Elevation plot (bedrock and ice surface)
    ax1.set_xlim([0, L])
    ax1.set_ylim([np.min(b)-100, np.max(h)+100])
    ax1.set_ylabel('Elevation (m)')
    ax1.set_title('Bed and Surface')

    line1, = ax1.plot(x, b, c='black', label='bedrock')  # bedrock plot
    line2, = ax1.plot(x, h, c='orangered', label='surface')  # ice surface plot
    line3, = ax1.plot(x, h0, c='black',linestyle='dotted', label='initial surface')  # initial ice surface plot

    # Second subplot: Horizontal velocity plot (ux)
    ax2.set_xlim([0, L])
    #ax2.set_ylim([min(-0.01,np.min(ux)), max(0.01,np.max(ux))])
    ax2.set_xlabel('Horizontal distance (m)')
    ax2.set_ylabel('Horizontal Velocity (m/s)')
    ax2.set_title('Surface Velocity')
    ax2_uz = ax2.twinx()
    #ax2_uz.set_ylim([min(-0.01,np.min(uz)), max(0.01,np.max(uz))])  # Adjust this range based on uz values
    ax2_uz.set_ylabel('Vertical velocity (m/s)')

    line4, = ax2.plot(x, np.zeros(len(x)), c='magenta', label='horizontal velocity')  # velocity plot
    line5, = ax2_uz.plot(x, np.zeros(len(x)), c='blue', label='vertical velocity')

    # Add legends to both axes
    ax2.legend(loc='upper left')  # Legend for u_x
    ax2_uz.legend(loc='upper right')  # Legend for u_z


    # Time text
    time_text = ax1.text(0.97, 0.95, '', transform=ax1.transAxes, va='top', ha='right')

    return fig,ax2,ax2_uz,line1,line2,line3,line4,line5,time_text

