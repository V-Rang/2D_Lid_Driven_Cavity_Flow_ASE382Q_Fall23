import numpy as np
from utils import *
T = 1
dt = 0.1
dx = 0.1
Re = 100    

#initializing
n = 10 #sqaure mesh, no of points in pressure grid along each dimension
nx,ny = n,n #keeping it square for now.
phi = np.zeros(int(n**2),dtype=np.float32)
u_n = np.zeros(((ny+2)*(nx+1)),dtype=np.float32)
v_n = np.zeros(((ny+1)*(nx+2)),dtype=np.float32)

G_x = np.zeros(len(u_n),dtype=np.float32)
G_y = np.zeros(len(v_n),dtype=np.float32)

#dont know how to update velocities for first time, dont have N(u_n_1). 
# Prof said use ERK2, how?

N_u_n = np.zeros(len(u_n),dtype=np.float32) 
N_v_n = np.zeros(len(v_n),dtype=np.float32)

#need to make only once
lhs_u_star = u_star_solver_lhs_operator(n,dt,Re,dx)
lhs_v_star = v_star_solver_lhs_operator(n,dt,Re,dx)
laplacian_pressure = create_laplacian_for_pressure(n,dx,dt)

#to start
u_n_1 = u_n
v_n_1 = v_n 
 
for t in range(10):

    u_star_rhs = u_star_solver_rhs(n,dt,Re,dx,u_n,v_n,u_n_1,v_n_1)
    v_star_rhs = v_star_solver_rhs(n,dt,Re,dx,u_n,v_n,u_n_1,v_n_1)
    
    u_star = scipy.sparse.linalg.spsolve(lhs_u_star,u_star_rhs)
    v_star = scipy.sparse.linalg.spsolve(lhs_v_star,v_star_rhs)

    Duv = Duv_maker(u_star,v_star,n,dx)
    phi = scipy.sparse.linalg.spsolve(laplacian_pressure,Duv)

    Gx,Gy = G_update_u_v(phi,dx)

    u_n_1,v_n_1 = u_n,v_n
    
    u_n,v_n = u_star - dt*Gx, v_star - dt*Gy


plot_velocities_quiver_plot(u_n,v_n,[0,1],[0,1],n,n)
# plot_pressure_contour(phi,[0,1],[0,1],nx,ny)
