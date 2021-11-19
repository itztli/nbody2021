#!/usr/bin/env python

#author: vdelaluz
#email: vdelaluz

import math
import random

G = 6.67e-11 # Nm2/kg2

class particle():
    m = 1.0 #kg
    p = [0.0, 0.0, 0.0] #m
    v = [0.0, 0.0, 0.0] #m/s
    r = [0.0, 0.0, 0.0] #m
    F0 = [0.0,0.0,0.0]
    
    def __init__(self, m, p, v):
        self.m = m
        self.p = p
        self.v = v

    def print_position(self):
        print(self.p[0],"\t",self.p[1],"\t",self.p[2])

    def u(self,p_j):
        self.r[0] = p_j.p[0] - self.p[0] 
        self.r[1] = p_j.p[1] - self.p[1]
        self.r[2] = p_j.p[2] - self.p[2]
        return self.r
    
class container():
    particles = []
    time=0.0

    def __init__(self, particles):
        self.particles = particles
        self.time=0.0
        
    def print_position(self):
        for particle in self.particles:
            particle.print_position()
        
class integrator():
    t0 = 0.0
    dt = 0.01 #sec
    steps = 1
    M=0
    
    def __init__(self, t0, dt, steps, M):
        self.t0 = t0
        self.dt = dt
        self.steps = steps
        self.M = M

    def riemman(self, container):
        for j in range(self.M):
            p_j = container.particles[j]            
            for i in range(self.M):
                p_i = container.particles[i]
                if (i == j):
                    continue
                
                u = p_j.u(p_i)
                ux = u[0]
                uy = u[1]
                uz = u[2]
                r = math.sqrt(ux*ux+uy*uy+uz*uz)
                #print("r",r)
                C = (G*p_j.m*p_i.m/math.pow(r,3))
                Fx = C*ux
                Fy = C*uy
                Fz = C*uz

                Vx= (self.dt/p_j.m)*Fx 
                Vy= (self.dt/p_j.m)*Fy
                Vz= (self.dt/p_j.m)*Fz
                
                p_j.v[0] = p_j.v[0] + Vx
                p_j.v[1] = p_j.v[1] + Vy
                p_j.v[2] = p_j.v[2] + Vz

        # 2 on the spot        
        for j in range(self.M):
            p_j = container.particles[j]            
            #v=d/t
            # d = vt
            p_j.p[0] = p_j.p[0] + p_j.v[0]*self.dt
            p_j.p[1] = p_j.p[1] + p_j.v[1]*self.dt
            p_j.p[2] = p_j.p[2] + p_j.v[2]*self.dt
                
        return container
                
        
    def trapecio(self, container):
        self.t0 = self.t0 + self.dt
        container.time = self.t0
        # sum F
        #1 on the spot        
        for j in range(self.M):
            p_j = container.particles[j]            
            for i in range(self.M):
                p_i = container.particles[i]
                if (i != j):
                    u = p_j.u(p_i)
                    ux = u[0]
                    uy = u[1]
                    uz = u[2]
                    r = math.sqrt(ux*ux+uy*uy+uz*uz)
                    print("r",r)
                    C = (G*p_j.m*p_i.m/math.pow(r,3))
                    Fx = C*ux
                    Fy = C*uy
                    Fz = C*uz

                    if (i==0):
                        V1 = (self.dt*Fx - p_j.m*p_j.v[0])/p_j.m
                        a0 = (V1 - p_j.v[0])/self.dt
                        F0x = p_j.m * a0
                        V1 = (self.dt*Fy - p_j.m*p_j.v[1])/p_j.m
                        a0 = (V1 - p_j.v[1])/self.dt
                        F0y = p_j.m * a0
                        V1 = (self.dt*Fz - p_j.m*p_j.v[2])/p_j.m
                        a0 = (V1 - p_j.v[2])/self.dt
                        F0z = p_j.m * a0
                    else:
                        F0x = p_j.F0[0] 
                        F0y = p_j.F0[1] 
                        F0z = p_j.F0[2] 

                    
                    Vx= (self.dt/2*p_j.m) * (F0x + Fx)
                    Vy= (self.dt/2*p_j.m) * (F0y + Fy)
                    Vz= (self.dt/2*p_j.m) * (F0z + Fz)

                    p_j.F0[0] = p_j.F0[0] + Fx
                    p_j.F0[1] = p_j.F0[1] + Fy
                    p_j.F0[2] = p_j.F0[2] + Fz

                    #p_j.F0[0] = Fx
                    #p_j.F0[1] = Fy
                    #p_j.F0[2] = Fz
                
                    p_j.v[0] = p_j.v[0] + Vx
                    p_j.v[1] = p_j.v[1] + Vy
                    p_j.v[2] = p_j.v[2] + Vz

        # 2 on the spot        
        for j in range(self.M):
            p_j = container.particles[j]            
            #v=d/t
            # d = vt
            p_j.p[0] = p_j.p[0] + p_j.v[0]*self.dt
            p_j.p[1] = p_j.p[1] + p_j.v[1]*self.dt
            p_j.p[2] = p_j.p[2] + p_j.v[2]*self.dt
                
        return container
        
    def get_t(self):
        return t0

dt = 0.1 
M = 100 #number of particles
N = 100 #number of iterations

x_As = -0.01
x_Bs = 0.01
y_As = -0.01
y_Bs = 0.01
z_As = -0.01
z_Bs = 0.01
r_circle = 2

universe = []

for i in range(M):
    x = random.random()
    y = random.random()
    z = random.random()
    universe.append(particle(1000.0, [x,y,z], [0,0,0]))


#p0 = particle(100.0, [0,0,0], [0,0,0])
#p1 = particle(100.0, [0.01,0,0], [0,0,0])
#p2 = particle(100.0, [0.0,0.0,0.01], [0,0,0])
#p3 = particle(100.0, [-0.01,0.0,0.0], [0,0,0])
#p4 = particle(100.0, [0.0,0.0,-0.01], [0,0,0])

c1 = container(universe)
s1 = integrator(0.0, dt, 1, M)



print(N,"\t",M,"\t",x_As,"\t",x_Bs,"\t",y_As,"\t",y_Bs,"\t",z_As,"\t",z_Bs,"\t",r_circle)

for i in range(N):
    c1 = s1.riemman(c1)
    c1.print_position()

# 1 cm de distancia sin v0
#  3.335330384569905e-07 0.0 0.0
#  0.009999332933923093 0.0 0.0

# 1 mm de distancia sin v0
# 3.7336274064142836e-05 0.0 0.0
# 0.0009253274518717142 0.0 0.0

# con suma de fuerzas
#0.03045098884816577 0.0 0.0
#-0.0007796116278593226 0.0 0.0

#con riemman
# 0.06845554368520183 0.0 0.0
# -0.06745554368520183 0.0 0.0

# con M0 = 10kg
#-0.011107696749138017 0.0 0.0
#0.11207696749138017 0.0 0.0

# v0 = 0 a 1 cm.
#3.411114224407899e-05 0.0 0.0
#0.009658888577559211 0.0 0.0


