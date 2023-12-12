import numpy as np
import matplotlib.pyplot as plt
import scipy


def plot_velocities_quiver_plot(uvals,vvals,xlims,ylims,nx,ny):
    #uvals,vvals are 1d arrays that have all the x, y component of velocities,
    #they are assumed to be stored in row major order:
    #if physical domain looks like:
    #[u_6, u_7, u_8]
    # [u_3, u_4, u_5]
    # [u_0, u_1, u_2] 
    #the uvector is supplied as [u_0, u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8]

    #nx -> num of points in p grid along x
    #ny -> num of points in p grid along y

    xlower,xupper = xlims
    ylower,yupper = ylims
    
    uvals,vvals = uvals.reshape(ny+2,nx+1),vvals.reshape(nx+1,ny+2)

    ucc = (uvals[1:-1,0:-1] + uvals[1:-1,1:])/2
    vcc = (vvals[0:-1,1:-1] + vvals[1:,1:-1])/2
    

    #to create mesh
    xvals = np.linspace(xlower,xupper,nx+1)
    yvals = np.linspace(ylower,yupper,ny+1)

    xcc = (xvals[0:-1] + xvals[1:])/2
    ycc = (yvals[0:-1] + yvals[1:])/2

    X,Y = np.meshgrid(xcc,ycc)

    fig,ax = plt.subplots()
    ax.scatter(X,Y, s=50, c='red', label='Coordinates')
    #scale will have to be adjusted to get optimuum plot
    ax.quiver(X, Y, ucc, vcc, color="C0", angles='xy',scale_units=None, scale=0.8, width=.015,label='Velocity')
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    ax.set(xlim=(xlower, xupper), ylim=(ylower, yupper))
    plt.show()

#########test for plot_velocities_quiver_plot function################    
#p -> 10,10 => u,v => (11X12)
# u,v = np.zeros(shape=(132)), np.zeros(shape=(132))
# import random
# for i in range(132):
#     u[i],v[i] = random.random(),random.random()

# plot_velocities_quiver_plot(u,v,[0,1],[0,1],10,10)
#########test for plot_velocities_quiver_plot function################    



def plot_pressure_contour(p_c,xlims,ylims,nx,ny):
    #uvals,vvals are 1d arrays that have all the x, y component of velocities,
    #they are assumed to be stored in row major order:
    #if physical domain looks like:
    #[u_6, u_7, u_8]
    # [u_3, u_4, u_5]
    # [u_0, u_1, u_2] 
    #the uvector is supplied as [u_0, u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8]

    #nx -> num of points in p grid along x
    #ny -> num of points in p grid along y

    xlower,xupper = xlims
    ylower,yupper = ylims
    
    # p_c,p_a = p_c.reshape(ny+2,nx+2),p_a.reshape(nx+2,ny+2)
    # p_c = p_c[1:-1,1:-1]
    # p_a = p_a[1:-1,1:-1]

    p_c = p_c.reshape(ny+2,nx+2)
    p_c = p_c[1:-1,1:-1]


    
    # print(p_c.shape)

    # ucc = (uvals[1:-1,0:-1] + uvals[1:-1,1:])/2
    # vcc = (vvals[0:-1,1:-1] + vvals[1:,1:-1])/2
    
    #to create mesh
    xvals = np.linspace(xlower,xupper,nx+1)
    yvals = np.linspace(ylower,yupper,ny+1)

    xcc = (xvals[0:-1] + xvals[1:])/2
    ycc = (yvals[0:-1] + yvals[1:])/2

    X,Y = np.meshgrid(xcc,ycc)

    fig,ax = plt.subplots()

    # min_range = min(np.min(p_c),np.min(p_a))
    # max_range = max(np.max(p_c),np.max(p_a))

    # print(np.min(p_c),np.min(p_a))
    # print(np.max(p_c),np.max(p_a))

    min_range = np.min(p_c)
    max_range = np.max(p_c)

          
    # print(min_range)
    # print(max_range)

    plt.contourf(X,Y,p_c,vmin=min_range,vmax=max_range)
    plt.show()

    # plt.contourf(X,Y,p_a,vmin=min_range,vmax=max_range)
    # plt.show()
    
    # fig,ax = plt.subplots()
    # ax.scatter(X,Y, s=50, c='red', label='Coordinates')
    # #scale will have to be adjusted to get optimuum plot
    # ax.quiver(X, Y, ucc, vcc, color="C0", angles='xy',scale_units=None, scale=10, width=.015,label='Velocity')
    # ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    # ax.set(xlim=(xlower, xupper), ylim=(ylower, yupper))
    # plt.show()


# def Laplace_maker(n):

    
#     #

#making Laplacian for pressure
# n = 2
def create_laplacian_for_pressure(n,h,k):
    #n  = number of unknowns along each dimension (assumed same in x and y dimension)
    #h = mesh spacing

    #Introducing ghost nodes => effective system = (n+2) X (n+2)
    #=> Total number of unknowns = (n+2) ^2
    # Therefore, return square matrix of size (n+2)^2 X (n+2)^2, make it sparse

    row,col,data =[],[],[]

    for ind in range((n+2)**2):
        new_phy_size = n+2
        #loop through rows of sparse matrix to be created, 
        # determine if equation is being written for a ghost node or interior node
        #if ghost, apply BC, if interior -> use 5 pt stencil

        #location in physical domain of row for which equation is being written
        i = ind // new_phy_size #y coordinate
        j = ind - i * new_phy_size #x coordinate in phy domain
        # print(ind, (i,j))
        
        #book keepging for corner, non contributing nodes
        # if( (i == 0 and j==0) or (i==0 and j==new_phy_size-1) or (i==new_phy_size-1 and j ==0) or (i==new_phy_size-1 and j ==new_phy_size-1) ):
        #     data.append(1)
        #     row.append(ind)
        #     col.append(ind) #corresp RHS = 0
        
        # else:
        if(i == 0): #bottom
            if(j == 0 or j == new_phy_size-1) : #left or right corner
                data.append(k/h**2)
                row.append(ind)
                col.append(ind)
            else:
                data.append(k/h**2)
                row.append(ind)
                col.append(ind)

                data.append(-k/h**2)
                row.append(ind)
                col.append(ind+new_phy_size)

        elif( i== new_phy_size - 1): #top
            if(j ==0 or j == new_phy_size - 1): #left or right corner
        
                data.append(k/h**2)
                row.append(ind)
                col.append(ind)
            else:
                data.append(k/h**2)
                row.append(ind)
                col.append(ind)

                data.append(-k/h**2)
                row.append(ind)
                col.append(ind-new_phy_size)

        else:
            if(j == 0): #left
                data.append(k/h**2)
                row.append(ind)
                col.append(ind)

                data.append(-k/h**2)
                row.append(ind)
                col.append(ind+1)
            elif(j == new_phy_size-1): #right
                data.append(k/h**2)
                row.append(ind)
                col.append(ind)

                data.append(-k/h**2)
                row.append(ind)
                col.append(ind-1)
            else: #interior
                # b[ind] = -5*np.pi**2 * np.sin(2*np.pi*(j-0.5)*h) * np.cos(np.pi *(i-0.5)*h)
                # T_a[ind] = np.sin(2*np.pi*(j-0.5)*h) * np.cos(np.pi *(i-0.5)*h)
                
                # if(ind == new_phy_size+1):
                    # data.append(1)
                    # row.append(ind)
                    # col.append(ind)
                    # b[ind] = 123
                # else:

                data.append(-4*k/h**2)
                row.append(ind)
                col.append(ind)

                data.append(k/h**2)
                row.append(ind)
                col.append(ind-1)

                data.append(k/h**2)
                row.append(ind)
                col.append(ind+1)

                data.append(k/h**2)
                row.append(ind)
                col.append(ind-new_phy_size)

                data.append(k/h**2)
                row.append(ind)
                col.append(ind+new_phy_size)
        
    laplacian_for_pressure = scipy.sparse.csr_array((data,(row,col)),shape=((n+2)**2 ,(n+2)**2   ) )

    
    return laplacian_for_pressure
    # return b,T_a,laplacian_for_pressure


# row,col,data = np.array(row),np.array(col),np.array(data)
# print(row)
# print(col)
# print(data)

# h = 0.25
# n = 20
# b = np.zeros(shape=((n+2)**2,),dtype=np.float32)

# T_a = np.zeros(shape=((n+2)**2,),dtype=np.float32)

# b,T_a,pr_lap = create_laplacian_for_pressure(n,h,b,T_a)

# pr_np_arr= pr_lap.toarray()
# x =np.linalg.pinv(pr_np_arr) @b
# # for i in range(len(b)):
#     print(i,x[i],T_a[i])


#########********************######################
# plot_pressure_contour(x,T_a,[0,1],[0,1],20,20) #Imp line
#########********************######################


# pr_np_arr= pr_lap.toarray()
# print(pr_np_arr)

# pr_np_arr *= 1/(h**2)




# b[5] = -5*np.pi**2 * np.sin(2*np.pi*h/2) * np.cos(np.pi *h/2)
# b[6] = -5*np.pi**2 * np.sin(2*np.pi*3*h/2) * np.cos(np.pi * h/2)
# b[9] = -5*np.pi**2 * np.sin(2*np.pi*h/2) * np.cos(np.pi * 3*h/2)
# b[10] = -5*np.pi**2 * np.sin(2*np.pi*3*h/2) * np.cos(np.pi * 3*h/2)



# T_a[5] = np.sin(2*np.pi*h/2) * np.cos(np.pi *h/2)
# T_a[6] = np.sin(2*np.pi*3*h/2) * np.cos(np.pi * h/2)
# T_a[9] = np.sin(2*np.pi*h/2) * np.cos(np.pi * 3*h/2)
# T_a[10] = np.sin(2*np.pi*3*h/2) * np.cos(np.pi *3*h/2)


# # print(b)
# print(np.linalg.det(pr_np_arr))



    # print(i,x[i],b[i])
    # print(i,b[i])
# print(x)

# print(test_matrix.toarray()) #as expected

##########end of Laplacian for pressure##############################

######rhs of Laplacian of pressure#######D(u,v)

#function to create Duv(u,v) for RHS of laplacian of pressure equation

def Duv_maker(u,v,n,h):
    #expecting square system of n X n, LHS is Laplacian of pressure which is (n+2)**2 X (n+2)**2 system with an equation for each of the (n+2)**2 unknowns (including boundary nodes)
    #expecting u and v to supplied in matrix format.
    nx,ny = n,n
    u = u.reshape( ((ny+2),(nx+1)) )
    v = v.reshape( ((ny+1),(nx+2))  )
    
    Duv = np.zeros((n+2)**2,dtype=np.float32)
    Duv_int = (u[1:-1,1:] - u[1:-1,0:-1])/h + (v[1:,1:-1] - v[0:-1,1:-1])/h
    Duv_int = Duv_int.flatten() 

    int_index = 0
    new_phy_size = n+2
    Duv = np.zeros((n+2)**2,dtype=np.float32)

    for ind in range(Duv.shape[0]):
        i = ind // new_phy_size
        j = ind - i * new_phy_size
        if(i ==0 or i== new_phy_size -1 or j ==0 or j == new_phy_size-1):
            Duv[ind] = 0
        else:
            Duv[ind] = Duv_int[int_index]
            int_index +=1
        
    return Duv


# nx = 2
# ny  = 2
# n = 2

# import random
# u = np.zeros(shape=((nx+2)*(ny+1)),dtype=np.float32)
# v = np.zeros(shape=((nx+1)*(ny+2)),dtype=np.float32)

# for i in range(u.shape[0]):
#     # u[i],v[i] = random.random(),random.random()
#     u[i] = random.random()
#     v[i] = random.random()

# u = u.reshape(((ny+2),(nx+1)))
# v = v.reshape(((ny+1),(nx+2)))

# h = 1 #mesh spacing

# Duv_calc = Duv_maker(u,v,n,h)

# print(u,'\n')
# print(v,'\n')
# print(Duv_calc) #works as expected


# print(v,'\n')

##############end of rhs of laplacian##################################################

#########create Gx, Gy to update un+1######################
def G_update_u_v(p,h):
    #p ->vector of length (n+2)**2
    n = int(np.sqrt(len(p))) # = n+2
    p = p.reshape((n,n))
    Gx = (p[:,1:] - p[:,0:-1])/h
    Gy = (p[1:,:] - p[0:-1,:])/h
    Gx,Gy = Gx.flatten(),Gy.flatten()

    return Gx,Gy
    #Gx and Gy are of same shape as matrix u and v and so can be subtracted as is in u(n+1) = u* - Gx



# n = 2
# h = 1
# new_phy_size = n+2
# import random
# p = np.zeros(shape=(new_phy_size**2),dtype=np.float32)

# for i in range(len(p)):
#     p[i] = random.random()

# p = p.reshape((n+2,n+2))

# Gx = (p[:,1:] - p[:,0:-1])/h
# Gy = (p[1:,:] - p[0:-1,:])/h

# print(p,'\n')
# gx,gy = G_update_u_v(p,h)
# print(gy,'\n')

#########end of create Gx, Gy to update un+1######################


###create function for N(u) and N(v) ###################
#has to only be created for interior nodes. Will be used to solve for u*, v*. Equations for u*, v* are different for ghost and boundary points,
# so N (scalar value for each interior node) is only needed at interior nodes

#for u4*, N(u_4) will be called; so have 0 value for boundary and ghost nodes.


def N_u_maker(u,v,n,h):
    #if supplied as vectors and not matrices, then make matrices first:
    nx,ny = n,n
    u = u.reshape((ny+2,nx+1))
    v = v.reshape((ny+1,nx+2))

    vcc = (v[0:-1,1:-1] + v[1:,1:-1])/2
    vcc_u_displ = (vcc[:,0:-1] + vcc[:,1:])/2

    # print(v,'\n')
    # print(vcc,'\n')
    # print(vcc_u_displ)

    N_u = np.zeros(((ny+2)*(nx+1)),dtype=np.float32)
    # h = 1
    N_u_int = -np.multiply(u[1:-1,1:-1],  (1/(2*h))*(  u[1:-1,2:] - u[1:-1,0:-2] )  ) -np.multiply( vcc_u_displ, (1/(2*h))*(  u[2:,1:-1] - u[0:-2,1:-1]    )  )
    N_u_int = N_u_int.flatten()
    
    N_u_int_index = 0

    #Need N_u vals at all u vals (set = 0 for ghost and boundary values)
    for ind in range((ny+2)*(nx+1)):
        i = ind // (nx+1)
        j = ind  - i * (nx+1)
        if(i == 0 or i == ny+1 or j == 0 or j == nx):
            N_u[ind] = 0
        else:
            N_u[ind] = N_u_int[N_u_int_index]
            N_u_int_index += 1
    
    return N_u #vector for all u (ghost, boundary = 0, internal = -u.del(u))


def N_v_maker(u,v,n,h):
    #if supplied as vectors and not matrices, then make matrices first:
    nx,ny = n,n
    u = u.reshape((ny+2,nx+1))
    v = v.reshape((ny+1,nx+2))
    ucc = (u[1:-1,0:-1]+ u[1:-1,1:])/2
    ucc_v_displ = (ucc[0:-1,:] + ucc[1:,:])/2

    N_v = np.zeros(((ny+1)*(nx+2)),dtype=np.float32)
    # h = 1
# N_v_int = -np.multiply(u[1:-1,1:-1],  (1/(2*h))*(  u[1:-1,2:] - u[1:-1,0:-2] )  ) -np.multiply( vcc_u_displ, (1/(2*h))*(  u[2:,1:-1] - u[0:-2,1:-1]    )  )
    N_v_int = -np.multiply( ucc_v_displ,(1/(2*h))*(v[1:-1,2:] - v[1:-1,0:-2] )) - np.multiply( v[1:-1,1:-1], (1/(2*h))*(v[2:,1:-1] - v[0:-2,1:-1]) )
    N_v_int = N_v_int.flatten()

    N_v_int_index = 0

    #Need N_v vals at all v vals (set = 0 for ghost and boundary values)
    for ind in range((ny+1)*(nx+2)):
        i = ind // (nx+2)
        j = ind  - i * (nx+2)
        if(i == 0 or i == ny or j == 0 or j == nx+1):
            N_v[ind] = 0
        else:
            N_v[ind] = N_v_int[N_v_int_index]
            N_v_int_index += 1

    return N_v


###################################################################
#function for LHS of u* solver
def u_star_solver_lhs_operator(n,k,Re,h):
    nx,ny = n,n
    row_u,col_u,data_u_lhs = [],[],[]
    factor = k/( 2*Re*(h**2) )    
    for ind in range((ny+2)*(nx+1)):
        i = ind//(nx+1) #will change for when making for u* and v* if not square mesh
        j = ind - i * (nx+1) 

        if(i == 0): #bottom
            row_u.append(ind)
            col_u.append(ind)
            data_u_lhs.append(1)

            row_u.append(ind)
            col_u.append(ind+nx+1)
            data_u_lhs.append(-1)


        elif( i == ny+1): #top
            row_u.append(ind)
            col_u.append(ind)
            data_u_lhs.append(1)

            row_u.append(ind)
            col_u.append(ind-(nx+1) )
            data_u_lhs.append(1) #utop + u(ind - (nx+1)) = 2, rhs wont be zero


        else:
            if( j == 0 or j == nx): #left or right boundary
                row_u.append(ind)
                col_u.append(ind)
                data_u_lhs.append(1)

            else: #interior
                row_u.append(ind)
                col_u.append(ind)
                data_u_lhs.append(1+4*factor)

                row_u.append(ind)
                col_u.append(ind-1)
                data_u_lhs.append(-factor)

                row_u.append(ind)
                col_u.append(ind+1)
                data_u_lhs.append(-factor)

                row_u.append(ind)
                col_u.append(ind-(nx+1))
                data_u_lhs.append(-factor)

                row_u.append(ind)
                col_u.append(ind+(nx+1))
                data_u_lhs.append(-factor)

    
    lhs_u_star_mat = scipy.sparse.csr_array((data_u_lhs,(row_u,col_u)),shape=( (ny+2)*(nx+1) ,(ny+2)*(nx+1)   ) )

    return lhs_u_star_mat

def u_star_solver_rhs(n,k,Re,h,u_n,v_n,u_n_1,v_n_1):
    N_u_n = N_u_maker(u_n,v_n,n,h)
    N_u_n_1 = N_u_maker(u_n_1,v_n_1,n,h)
    
    nx,ny = n,n
    data_u_rhs = np.zeros(((ny+2)*(nx+1)),dtype=np.float32)

    for ind in range((ny+2)*(nx+1)):
        i = ind//(nx+1)
        j = ind - i * (nx+1)

        if(i == 0): #bottom
            data_u_rhs[ind] = 0
        elif (i == ny+1): #top
                data_u_rhs[ind] = 2.
        else: 
            if(j ==0 or j == nx): #left or right
                data_u_rhs[ind] = 0
            else: #interior
                data_u_rhs[ind] = u_n[ind] + (k/2)*(3*N_u_n[ind] - N_u_n_1[ind]) + (k/(2*Re))*(1/h**2)*(u_n[ind-1] + u_n[ind+1] + u_n[ind+(nx+1)] + u_n[ind-(nx+1)] - 4*u_n[ind] )

    return data_u_rhs

def v_star_solver_lhs_operator(n,k,Re,h):
    nx,ny = n,n
    row_v,col_v,data_v_lhs = [],[],[]
    factor = k/( 2*Re*(h**2) )    
    for ind in range((ny+1)*(nx+2)):
        i = ind//(nx+2) #will change for when making for u* and v* even if sqaure mesh b/c of indexing
        j = ind - i * (nx+2) 

        if(i == 0 or i == ny): #bottom or top
            if(j == 0): #left
                
                row_v.append(ind)
                col_v.append(ind)
                data_v_lhs.append(1)
                
                row_v.append(ind)
                col_v.append(ind+1)
                data_v_lhs.append(-1)

            elif( j == nx+1): #right
                
                row_v.append(ind)
                col_v.append(ind)
                data_v_lhs.append(1)
                
                row_v.append(ind)
                col_v.append(ind-1)
                data_v_lhs.append(-1)
            else:
                row_v.append(ind)
                col_v.append(ind)
                data_v_lhs.append(1)
        
        else:
            if( j == 0): #left
                row_v.append(ind)
                col_v.append(ind)
                data_v_lhs.append(1)
                
                row_v.append(ind)
                col_v.append(ind+1)
                data_v_lhs.append(-1)
            
            elif(j == nx+1): #right
                row_v.append(ind)
                col_v.append(ind)
                data_v_lhs.append(1)
                
                row_v.append(ind)
                col_v.append(ind-1)
                data_v_lhs.append(-1)

            else: #interior

                row_v.append(ind)
                col_v.append(ind)
                data_v_lhs.append(1+4*factor)

                row_v.append(ind)
                col_v.append(ind-1)
                data_v_lhs.append(-factor)

                row_v.append(ind)
                col_v.append(ind+1)
                data_v_lhs.append(-factor)

                row_v.append(ind)
                col_v.append(ind-(nx+2))
                data_v_lhs.append(-factor)

                row_v.append(ind)
                col_v.append(ind+(nx+2))
                data_v_lhs.append(-factor)
    
    lhs_v_star_mat = scipy.sparse.csr_array((data_v_lhs,(row_v,col_v)),shape=( (ny+2)*(nx+1) ,(ny+2)*(nx+1)   ) )

    return lhs_v_star_mat

def v_star_solver_rhs(n,k,Re,h,u_n,v_n,u_n_1,v_n_1):
    N_v_n = N_v_maker(u_n, v_n, n, h)
    N_v_n_1 = N_v_maker(u_n_1, v_n_1, n, h)
    
    nx,ny = n,n
    data_v_rhs = np.zeros(((ny+1)*(nx+2)),dtype=np.float32)

    for ind in range((ny+1)*(nx+2)):
        i = ind//(nx+2)
        j = ind - i * (nx+2)

        if(i == 0 or i== ny or j ==0 or j == nx+1): #ghost or boundary point

            data_v_rhs[ind] = 0

        else: #interior
            data_v_rhs[ind] = v_n[ind] + (k/2)*(3*N_v_n[ind] - N_v_n_1[ind]) + (k/(2*Re))*(1/h**2)*(v_n[ind-1] + v_n[ind+1] + v_n[ind+(nx+2)] + v_n[ind-(nx+2)] - 4*v_n[ind] )

    return data_v_rhs
