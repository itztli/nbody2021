#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from time import sleep
#from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches


ccolor=['k','r','g']
# Data for plotting

#t = np.arange(0.0, 2.0, 0.01)
#s = 1 + np.sin(2 * np.pi * t)

with open("view3D.dat") as f:
    data = f.readlines()

iterations, nparticle, x_As, x_Bs, y_As, y_Bs,z_As, z_Bs, r_circle = data[0].split("\t")
x_A = float(x_As)
x_B = float(x_Bs)
y_A = float(y_As)
y_B = float(y_Bs)
z_A = float(z_As)
z_B = float(z_Bs)
r_circle=float(r_circle)
skip = 0

#iterations = iterations[1:]

nparticle = int(nparticle)
iterations = int(iterations)

print(nparticle,iterations)

X = [[] for x in range(nparticle)]
Y = [[] for x in range(nparticle)]
Z = [[] for x in range(nparticle)]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

iterator = 0
new_sick_vector = []
bar_x=[]


#line in data[1:]:


for i in range(iterations):
    ax.cla()
    ax.set_xlim3d(x_A,x_B)
    ax.set_ylim3d(y_A,y_B)
    ax.set_zlim3d(z_A,z_B)
    
    for j in range(nparticle):
        position = data[i*nparticle+j+1]
        x,y,z = position.split("\t")
        x = float(x)
        y=float(y)
        z=float(z)
        print(j,x,y,z)
        ax.scatter(x, y, z)

    plt.draw()        
    fig.savefig("animation/"+("%08d"%i)+"-view3D.png")
    
