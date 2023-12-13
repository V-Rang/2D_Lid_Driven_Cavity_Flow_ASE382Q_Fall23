import numpy as np
import matplotlib.pyplot as plt
from utils import u_star_solver_lhs_operator,u_star_solver_rhs
from utils import v_star_solver_lhs_operator,v_star_solver_rhs
from utils import create_laplacian_for_pressure
from utils import Duv_maker, G_update_u_v
import scipy



n = 19 #sqaure mesh, no of points in pressure grid along each dimension
nx,ny = n,n #keeping it square for now. [0,1] x [0,1]

dx = 1/(nx-1)

#Specify time-step size and number of time iterations, vary to achieve convergence.
# dt = T/(nt-1) #time 

dt = 0.05
T =  100
nt = 1 + int(T/(dt))

Re = 100 # Reynolds Number
tol = 1e-3 
# tolerance for early termination, i.e. if difference x-comp
# of velocity from one t step to another is less than tol   

phi = np.zeros(int( (n+2)**2),dtype=np.float64)
u_n = np.zeros(((ny+2)*(nx+1)),dtype=np.float64)
v_n = np.zeros(((ny+1)*(nx+2)),dtype=np.float64)


G_x = np.zeros(len(u_n),dtype=np.float64)
G_y = np.zeros(len(v_n),dtype=np.float64)

N_u_n = np.zeros(len(u_n),dtype=np.float64) 
N_v_n = np.zeros(len(v_n),dtype=np.float64)

#to start the loop
u_n_1 = u_n
v_n_1 = v_n 

#operators that need to be made only once
lhs_u_star = u_star_solver_lhs_operator(n,dt,Re,dx)
lhs_v_star = v_star_solver_lhs_operator(n,dt,Re,dx)
laplacian_pressure = create_laplacian_for_pressure(n,dx,dt)


# err_vals = []
for t in range(nt):

    u_star_rhs = u_star_solver_rhs(n, dt, Re, dx, u_n, v_n, u_n_1, v_n_1)
    v_star_rhs = v_star_solver_rhs(n, dt, Re, dx, u_n, v_n, u_n_1, v_n_1)
    
    u_star = scipy.sparse.linalg.spsolve(lhs_u_star,u_star_rhs)
    v_star = scipy.sparse.linalg.spsolve(lhs_v_star,v_star_rhs)

    Duv = Duv_maker(u_star,v_star,n,dx)

    phi = scipy.sparse.linalg.spsolve(laplacian_pressure,Duv)

    Gx,Gy = G_update_u_v(phi,dx)

    u_n_1,v_n_1 = u_n,v_n
    
    u_n,v_n = u_star - dt*Gx, v_star - dt*Gy


    # err_vals.append(np.linalg.norm(u_n - u_n_1,2)/(np.linalg.norm(u_n,2)*dt)   )
    if(np.linalg.norm(u_n - u_n_1) < tol):
        print(f"Converged, no of iterations = {t}")
        break   


# plt.semilogy(err_vals)

plt.show()
u_n,v_n = u_n.reshape(ny+2,nx+1),v_n.reshape(nx+1,ny+2)

exp_val = (u_n[1:-1,1:] - u_n[1:-1,0:-1])/dx + (v_n[1:,1:-1] - v_n[0:-1,1:-1])/dx
print(np.linalg.norm(exp_val))


x = np.linspace(0,1,nx+1)
y = np.linspace(0,1,ny+1)
xcc = (x[0:-1] + x[1:])/2
ycc = (y[0:-1] + y[1:])/2
X,Y = np.meshgrid(xcc,ycc)

phi_pl1 = np.reshape(phi,(ny+2,nx+2))
phi_pl1 = phi_pl1[1:-1,1:-1]

ucc = (u_n[1:-1,0:-1] + u_n[1:-1,1:])/2
vcc = (v_n[0:-1,1:-1] + v_n[1:,1:-1])/2



plt.figure()
plt.contourf(X, Y, phi_pl1)
plt.colorbar()
# plt.quiver(X, Y, ucc, vcc, color="black")
plt.streamplot(X, Y, ucc, vcc, color="black")
plt.title(f"N = {n+1} points along each dimension and Re = {Re}")
plt.xlim(0,1)
plt.ylim(0,1)
savepath = f"{n+1}_{Re}.png"
plt.savefig(savepath)
plt.show()

