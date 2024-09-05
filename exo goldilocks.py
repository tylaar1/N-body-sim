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




# Time grid

t = np.linspace(0, 10000000, 5000)

# Solve the ODE


initial_state_1 = ([0,0,0,0,0,0,0.0788*au,0,0,0,v1,0,0.5141*au,0,0,0,v2,0,2.845*au,0,0,0,v3,0,0.17*au,0,0,0,35e3,0])

solution = odeint(odes,initial_state_1,t, args=(G, masses))

goldilocks_star_dist_1=abs(np.sqrt((solution[:,0]-solution[:,24])**2+(solution[:,1]-solution[:,25])**2
                       +(solution[:,2]-solution[:,26])**2))

initial_state_2 = ([0,0,0,0,0,0,0.0788*au,0,0,0,v1,0,0.5141*au,0,0,0,v2,0,2.845*au,0,0,0,v3,0,0.17*au,0,0,0,40e3,0])

solution = odeint(odes,initial_state_2,t, args=(G, masses))

goldilocks_star_dist_2=abs(np.sqrt((solution[:,0]-solution[:,24])**2+(solution[:,1]-solution[:,25])**2
                       +(solution[:,2]-solution[:,26])**2))

initial_state_3 = ([0,0,0,0,0,0,0.0788*au,0,0,0,v1,0,0.5141*au,0,0,0,v2,0,2.845*au,0,0,0,v3,0,0.17*au,0,0,0,45e3,0])

solution = odeint(odes,initial_state_3,t, args=(G, masses))

goldilocks_star_dist_3=abs(np.sqrt((solution[:,0]-solution[:,24])**2+(solution[:,1]-solution[:,25])**2
                       +(solution[:,2]-solution[:,26])**2))

initial_state_4 = ([0,0,0,0,0,0,0.0788*au,0,0,0,v1,0,0.5141*au,0,0,0,v2,0,2.845*au,0,0,0,v3,0,0.17*au,0,0,0,50e3,0])

solution = odeint(odes,initial_state_4,t, args=(G, masses))

goldilocks_star_dist_4=abs(np.sqrt((solution[:,0]-solution[:,24])**2+(solution[:,1]-solution[:,25])**2
                       +(solution[:,2]-solution[:,26])**2))



plt.plot(t,goldilocks_star_dist_1,label='35km/s')
plt.plot(t,goldilocks_star_dist_2,label='40km/s')
plt.plot(t,goldilocks_star_dist_3,label='45km/s')
plt.plot(t,goldilocks_star_dist_4,label='50km/s')

plt.ylabel('distance (m)')
plt.xlabel('time (seconds)')

plt.legend(loc="upper left")

