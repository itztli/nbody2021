#!/usr/bin/env python

#author: vdelaluz
#email: vdelaluz

class particle():
    m = 1.0 #kg
    p = [0.0, 0.0, 0.0] #m
    v = [0.0, 0.0, 0.0] #m/s
    r = [0.0, 0.0, 0.0] #m
    
    def __init__(self, m, p, v):
        self.m = m
        self.p = p
        self.v = v

    def print_position(self,time):
        print(time, self.p)

    def u(self,p_j):
        self.r[0] = p_j[0] - self.p[0] 
        self.r[1] = p_j[1] - self.p[1]
        self.r[1] = p_j[1] - self.p[1]
        return r
    
class container():
    particles = []
    time=0.0

    def __init__(self, particles):
        self.particles = particles
        self.time=0.0
        
    def print_position(self):
        for particle in self.particles:
            particle.print_position(self.time)
        
class integrator():
    t0 = 0.0
    dt = 1 #sec
    steps = 1
    M=0
    
    def __init__(self, t0, dt, steps, M):
        self.t0 = t0
        self.dt = dt
        self.steps = steps
        self.M = M
        
    def trapecio(self, container):
        self.t0 = self.t0 + self.dt
        container.time = self.t0

        # sum F
        for j in range(self.M):
            p_j = container.particles[j]            
            for i in range(self.M):
                p_i = container.particles[i]
                u = p_j.u(p_i)
                ux = u[0]
                uy = u[1]
                uz = u[2]
                r = math.sqrt(ux*ux+uy*uy+uz*uz)
                C = (G*p_j.m*p_i.m/math.pow(r,3))
                Fx = C*ux
                Fy = C*uy
                Fz = C*uz

                Vx=
                Vy=
                Vz=
                
        return container
        
    def get_t(self):
        return t0

M = 2 #number of particles
N = 100 #number of iterations

p0 = particle(1, [0,0,0], [0,0,0])
p1 = particle(1, [1,0,0], [0,0,0])

c1 = container([p0,p1])

s1 = integrator(0.0, 0.1, 1, M)

print(M,N)

for i in range(N):
    c1 = s1.trapecio(c1)
    c1.print_position()
