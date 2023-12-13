# We formulate our Laplacian that we use in computing our phi values (pressure)
# The order of convergence is compared by manufacturing a known vector:
# T(x,y) = sin(2pi x)cos(pi y)
#The expected rhs when this vector is multiplied by the Laplacian is 
# rhs(x,y) = -5 * pi**2 * sin(2pi x)cos(pi y)
#We compare the error of our numerical solution against this analytical solution
#The expected order of convergence is 2.


import numpy as np
import matplotlib.pyplot as plt
from utils import *

T = 1
dt = 0.1
dx = 0.1
Re = 5000
tol = 1e-3

#initializing
n = 15 #sqaure mesh, no of points in pressure grid along each dimension
nx,ny = n,n #keeping it square for now. [0,1] x [0,1]
dt = 0.1

# dx = 0.1
dx = 1/(nx-1)

phi = np.zeros(int( (n+2)**2),dtype=np.float64)
u_n = np.zeros(((ny+2)*(nx+1)),dtype=np.float64)
v_n = np.zeros(((ny+1)*(nx+2)),dtype=np.float64)



G_x = np.zeros(len(u_n),dtype=np.float64)
G_y = np.zeros(len(v_n),dtype=np.float64)

N_u_n = np.zeros(len(u_n),dtype=np.float64) 
N_v_n = np.zeros(len(v_n),dtype=np.float64)

#need to make only once
# dt,dx,Re = 1,1,1 #to test what the matrix looks like

lhs_u_star = u_star_solver_lhs_operator(n,dt,Re,dx)
lhs_v_star = v_star_solver_lhs_operator(n,dt,Re,dx)


laplacian_pressure = create_laplacian_for_pressure(n,dx,dt)


#validating the operator
# N_vals = [10,20,30,40,50,6070]


N_vals = np.linspace(10,100,10,dtype=np.int64)
err_vals = []
for i in range(len(N_vals)):
    
    n,dx,dt =  N_vals[i],1/(n-1),1

    laplacian_pressure = create_laplacian_for_pressure(n,dx,dt)
    lp_array = laplacian_pressure.toarray() 

    #shape = (n+2)*(n+2)
    # print(lp_array)

    T_analytical,rhs_analytical = np.zeros((n+2)**2), np.zeros((n+2)**2)

    for ind in range((n+2)**2):
        i  = ind//(n+2)
        j = ind - i * (n+2)
        x = j*dx
        y = i*dx
        T_analytical[ind] = np.sin(2*np.pi*x)*np.cos(np.pi*y)
        rhs_analytical[ind] = -5*np.pi**2 * np.sin(2*np.pi*x)*np.cos(np.pi*y)

    rhsa_mod = np.reshape(rhs_analytical,(n+2,n+2))
    rhsa_mod = rhsa_mod[1:-1,1:-1]

    rhs_numerical = lp_array@T_analytical
    rhsn_mod = np.reshape(rhs_numerical,(n+2,n+2))
    rhsn_mod = rhsn_mod[1:-1,1:-1]
    err_vals.append(np.linalg.norm(rhsa_mod.flatten() - rhsn_mod.flatten(),np.inf))



x_values,y_values = np.log(N_vals),np.log(err_vals)
slope, intercept = np.polyfit(x_values, y_values, 1)

best_fit_line = slope * x_values + intercept

plt.scatter(x_values, y_values, label='Data')

plt.plot(x_values, best_fit_line, label=f'Best Fit Line (Slope = {slope:.2f})', color='red')
plt.xlabel(r'$\log(1/N)$')
plt.ylabel(r'$\log(\|\epsilon\|_\infty)$')
plt.legend()
plt.title('Log plot for validating order of convergence of Laplacian - log(error) vs log(1/N)')


# Show the plot
plt.grid(True)

plt.savefig('laplace_validate_convergence.png')


plt.show()
