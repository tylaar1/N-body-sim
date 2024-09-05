# -*- coding: utf-8 -*-

"""
Created on Tue Nov 28 16:42:22 2023

@author: ppytr13
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
    
    
    return derivatives

# Parameters
G = 6.6743e-11
msun = 1.989e30
au = 1.496e11

# N bodies

masses = [msun, 3.29e23, 4.87e24,5.97e24,6.41e23,1.90e27,5.68e26,8.69e25,1.02e26]

# Initial conditions - from nasa horizons system 2023-Dec-05 00:00:00.0000
initial_state = ([-1.211268576898178E+09, -4.057617460661397E+08, 3.163898964273166E+07, 8.072761040686311E+00, -1.264839836757296E+01 ,-7.222027978150917E-02,
                  5.266070322317566E+10, -9.789098706517318E+09, -5.676537933956491E+09, -1.065013371773722E+03, 5.015839676728469E+04 ,4.198385617078390E+03,
                  -8.613457506479302E+10, 6.523388841115820E+10, 5.833237666525669E+09,-2.154216553367116E+04,-2.789947444238368E+04,8.604141874895621E+02,
                  4.376081683153060E+10,1.400096748463262E+11,2.355651078157127E+07,-2.884599232857679E+04,8.961652677752641E+03,3.900157580347674E-01,
                  -1.005980420167290E+11,-2.041985387166618E+11,-1.801413551008239E+09,2.269769778598640E+04,-8.557296595341191E+03,-7.357015224301335E+02,
                  5.430111710839081E+11,5.081446334995103E+11,-1.425681890934390E+10,-9.072275623426602E+03,1.015860665390452E+04,1.608319680116543E+02,
                  1.337087088645709E+12,-5.770757051583652E+11,-4.320194250054500E+10,3.288271020485106E+03,8.849911704641903E+03,-2.848064582601899E+02,
                  1.847007858228591E+12,2.279227642566298E+12,-1.546323032258439E+10,-5.340863295440040E+03,3.970210761445738E+03,8.393562685451617E+01,
                  4.462559299006256E+12,-2.810810617798207E+11,-9.705588037391016E+10,3.057387830169898E+02,5.456761129464996E+03,-1.194176868401768E+02])

# Time grid

t = np.linspace(0, 6000000000, 50000)

# Solve the ODE
solution = odeint(odes,initial_state,t, args=(G, masses))

bodies=['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Passing Star']

# Plot the trajectories
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(len(masses)):
    ax.plot(solution[:, i*6], solution[:, (i*6)+1], solution[:, (i*6)+2], label=bodies[i])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.3))

ax.set_xlabel('X distance (Tm)')
ax.set_ylabel('Y distance (Tm)')
ax.set_zlabel('Z distance (Tm)')

ax.set_xlim(-6e12,6e12)
ax.set_ylim(-6e12,6e12)
ax.set_zlim(-6e12,6e12)