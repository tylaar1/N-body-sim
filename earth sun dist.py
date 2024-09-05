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
au = 1.496e11

# N bodies
#moon passing through solar system
masses = [msun, 3.29e23, 4.87e24,5.97e24,6.41e23,7.34e22]

# Initial conditions - from nasa horizons system 2023-Dec-05 00:00:00.0000
initial_state = ([-1.211268576898178E+09, -4.057617460661397E+08, 3.163898964273166E+07, 8.072761040686311E+00, -1.264839836757296E+01 ,-7.222027978150917E-02,
                  5.266070322317566E+10, -9.789098706517318E+09, -5.676537933956491E+09, -1.065013371773722E+03, 5.015839676728469E+04 ,4.198385617078390E+03,
                  -8.613457506479302E+10, 6.523388841115820E+10, 5.833237666525669E+09,-2.154216553367116E+04,-2.789947444238368E+04,8.604141874895621E+02,
                  4.376081683153060E+10,1.400096748463262E+11,2.355651078157127E+07,-2.884599232857679E+04,8.961652677752641E+03,3.900157580347674E-01,
                  -1.005980420167290E+11,-2.041985387166618E+11,-1.801413551008239E+09,2.269769778598640E+04,-8.557296595341191E+03,-7.357015224301335E+02,
                  1e12,1e12,0,-1e4,-0.9e4,0])

# Time grid

t = np.linspace(0, 1000000000, 1000000)

# Solve the ODE
solution = odeint(odes,initial_state,t, args=(G, masses))


earth_sun_dist_1=abs(np.sqrt((solution[:,0]-solution[:,18])**2+(solution[:,1]-solution[:,19])**2
                       +(solution[:,2]-solution[:,20])**2))
#jupiter mass passing through solar system
masses = [msun, 3.29e23, 4.87e24,5.97e24,6.41e23,1.90e27]

# Initial conditions - from nasa horizons system 2023-Dec-05 00:00:00.0000
initial_state = ([-1.211268576898178E+09, -4.057617460661397E+08, 3.163898964273166E+07, 8.072761040686311E+00, -1.264839836757296E+01 ,-7.222027978150917E-02,
                  5.266070322317566E+10, -9.789098706517318E+09, -5.676537933956491E+09, -1.065013371773722E+03, 5.015839676728469E+04 ,4.198385617078390E+03,
                  -8.613457506479302E+10, 6.523388841115820E+10, 5.833237666525669E+09,-2.154216553367116E+04,-2.789947444238368E+04,8.604141874895621E+02,
                  4.376081683153060E+10,1.400096748463262E+11,2.355651078157127E+07,-2.884599232857679E+04,8.961652677752641E+03,3.900157580347674E-01,
                  -1.005980420167290E+11,-2.041985387166618E+11,-1.801413551008239E+09,2.269769778598640E+04,-8.557296595341191E+03,-7.357015224301335E+02,
                  1e12,1e12,0,-1e4,-0.9e4,0])

# Time grid

t = np.linspace(0, 1000000000, 1000000)

# Solve the ODE
solution = odeint(odes,initial_state,t, args=(G, masses))


earth_sun_dist_2=abs(np.sqrt((solution[:,0]-solution[:,18])**2+(solution[:,1]-solution[:,19])**2
                       +(solution[:,2]-solution[:,20])**2))
#solar mass passing through solar system
masses = [msun, 3.29e23, 4.87e24,5.97e24,6.41e23,msun]

# Initial conditions - from nasa horizons system 2023-Dec-05 00:00:00.0000
initial_state = ([-1.211268576898178E+09, -4.057617460661397E+08, 3.163898964273166E+07, 8.072761040686311E+00, -1.264839836757296E+01 ,-7.222027978150917E-02,
                  5.266070322317566E+10, -9.789098706517318E+09, -5.676537933956491E+09, -1.065013371773722E+03, 5.015839676728469E+04 ,4.198385617078390E+03,
                  -8.613457506479302E+10, 6.523388841115820E+10, 5.833237666525669E+09,-2.154216553367116E+04,-2.789947444238368E+04,8.604141874895621E+02,
                  4.376081683153060E+10,1.400096748463262E+11,2.355651078157127E+07,-2.884599232857679E+04,8.961652677752641E+03,3.900157580347674E-01,
                  -1.005980420167290E+11,-2.041985387166618E+11,-1.801413551008239E+09,2.269769778598640E+04,-8.557296595341191E+03,-7.357015224301335E+02,
                  1e12,1e12,0,-1e4,-0.9e4,0])

# Time grid

t = np.linspace(0, 1000000000, 1000000)

# Solve the ODE
solution = odeint(odes,initial_state,t, args=(G, masses))


earth_sun_dist_3=abs(np.sqrt((solution[:,0]-solution[:,18])**2+(solution[:,1]-solution[:,19])**2
                       +(solution[:,2]-solution[:,20])**2))


#supergiant 10x msun passing through solar system
masses = [msun, 3.29e23, 4.87e24,5.97e24,6.41e23,10*msun]

# Initial conditions - from nasa horizons system 2023-Dec-05 00:00:00.0000
initial_state = ([-1.211268576898178E+09, -4.057617460661397E+08, 3.163898964273166E+07, 8.072761040686311E+00, -1.264839836757296E+01 ,-7.222027978150917E-02,
                  5.266070322317566E+10, -9.789098706517318E+09, -5.676537933956491E+09, -1.065013371773722E+03, 5.015839676728469E+04 ,4.198385617078390E+03,
                  -8.613457506479302E+10, 6.523388841115820E+10, 5.833237666525669E+09,-2.154216553367116E+04,-2.789947444238368E+04,8.604141874895621E+02,
                  4.376081683153060E+10,1.400096748463262E+11,2.355651078157127E+07,-2.884599232857679E+04,8.961652677752641E+03,3.900157580347674E-01,
                  -1.005980420167290E+11,-2.041985387166618E+11,-1.801413551008239E+09,2.269769778598640E+04,-8.557296595341191E+03,-7.357015224301335E+02,
                  1e12,1e12,0,-1e4,-0.9e4,0])

# Time grid

t = np.linspace(0, 1000000000, 1000000)

# Solve the ODE
solution = odeint(odes,initial_state,t, args=(G, masses))


earth_sun_dist_4=abs(np.sqrt((solution[:,0]-solution[:,18])**2+(solution[:,1]-solution[:,19])**2
                       +(solution[:,2]-solution[:,20])**2))

#plotting displacements against one another
plt.plot(t,np.log10(earth_sun_dist_1),label='moon')
plt.plot(t,np.log10(earth_sun_dist_2),label='jupiter')
plt.plot(t,np.log10(earth_sun_dist_3),label='sun')
plt.plot(t,np.log10(earth_sun_dist_4),label='supergiant')
plt.ylabel('log of distance')
plt.xlabel('time (seconds)')

plt.legend(loc="lower center")
