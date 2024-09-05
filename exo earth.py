# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 15:53:18 2023

@author: ppytr13
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:51:06 2023

@author: super
"""



import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
plt.close('all')
def odes(y, t, G, masses):

    derivatives=np.zeros_like(y)
    vx=[]
    vy=[]
    vz=[]
    r=[]
    axx=[]
    ay=[]
    az=[]
    for i in range(len(masses)):
        #appending velocities for re-use in next iteration
        vx.append(y[(i*6)+3])
        vy.append(y[(i*6)+4])
        vz.append(y[(i*6)+5])
        

                              
    for i in range(len(masses)):
        #initialising accelerations 
        axi=0
        ayi=0
        azi=0
        for j in range(len(masses)): 
            if i != j:
            #calculation of accelerations and masses
               r = (np.sqrt((y[j*6]-y[i*6])**2 + (y[(j*6)+1]-y[(i*6)+1])**2 + (y[(j*6)+2]-y[(i*6)+2])**2))
               axi += (G*masses[j]*(y[j*6]-y[i*6])/(r**3))
               ayi += (G*masses[j]*(y[(j*6)+1]-y[(i*6)+1])/(r**3))
               azi += (G*masses[j]*(y[(j*6)+2]-y[(i*6)+2])/(r**3))
                   
        axx.append(axi)
        ay.append(ayi)
        az.append(azi)
        #setting up output in a form readable by odeint
        derivatives[i*6]=vx[i]
        derivatives[(i*6)+1]=vy[i]
        derivatives[(i*6)+2]=vz[i]
        derivatives[(i*6)+3]=axx[i]
        derivatives[(i*6)+4]=ay[i]
        derivatives[(i*6)+5]=az[i]
    
    #print(r,f)
    
    return derivatives

# Parameters
G = 6.6743e-11
msun = 1.989e30
mearth=5.97e24
au = 1.496e11

# N bodies

masses = [0.389*msun,2.64*mearth,4.1*mearth,14.2*mearth,mearth]

#approximating orbits as circular motion
v1=2*np.pi*(0.0788*au)/(12.94*86400)
v2=2*np.pi*(0.5141*au)/(215.62*86400)
v3=2*np.pi*(2.845*au)/(2806*86400)

# Initial conditions - from nasa horizons system 2023-Dec-05 00:00:00.0000
initial_state = ([0,0,0,0,0,0,0.0788*au,0,0,0,v1,0,0.5141*au,0,0,0,v2,0,2.845*au,0,0,0,v3,0,au,0,0,0,30e3,0])

# Time grid

t = np.linspace(0, 250000000, 50000)

# Solve the ODE
solution = odeint(odes,initial_state,t, args=(G, masses))

bodies=['Lalande 21185','Lalande 21185-b','Lalande 21185-d','Lalande 21185-c','Earth']

# Plot the trajectories
fig = plt.figure()
ax1=plt.axes([0.27,0.51,0.46,0.45])
ax2=plt.axes([0.27,0.09,0.46,0.3])


for i in range(len(masses)):
    ax1.plot(solution[:, i*6], solution[:, (i*6)+1], label=bodies[i])
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.3))

ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('Distance (m)')


ax1.set_xlim(-6e11,6e11)
ax1.set_ylim(-6e11,6e11)


for i in range(len(masses)):
    ax2.plot(t, solution[:, i*6])

ax2.set_xlabel('Time (s)')
ax2.set_ylabel('X Position (m)')
