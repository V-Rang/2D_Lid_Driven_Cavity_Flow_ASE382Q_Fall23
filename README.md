We will be simulating 2D lid driven cavity flow using fractional step method of Kim and Moin [1985].

main.py contains code for initialization of variables and running of the loop, results for the norm of the divergence of velocity and quiver, streamplots are plotted along with pressure contours.

utils.py contains all the functions that are used to implement the steps as listed as by Kim and Moin.

validate_laplacian.py contains a test code to confirm the 2nd order convergence of the Laplacian that is used to solve the Poisson equation for pressure.
