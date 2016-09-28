import numpy as np
import scitools.std as sci
import os, shutil, sys

#xk = 10
#alfa = 0
#time_steps = 200
#dt = 0.5
#Nx = 200
#depth = 1
#amplitude_factor = 0.2

#amplitude = amplitude_factor * depth



#######################################################################
###                     DOMAIN                                      ###
#######################################################################


def domain(xk, alfa, Nx, depth):

    if alfa == 0.0:
        dx = float(xk) / Nx
        bpu = depth*np.ones(Nx+1)
        xu = np.linspace(0, xk, Nx+1)
        bpeta = depth*np.ones(Nx+2)
        xeta = np.linspace(-0.5*dx, xk+0.5*dx, Nx+2)
    else:
        bpu = np.zeros(Nx+1)
        bpeta = np.zeros(Nx+2)

        degree_incline = alfa / (90.0 - alfa)
        xs = xk + depth / degree_incline   

        def slope(coordinate):
            if coordinate <= xk:
                return depth
            else:
                return depth + degree_incline * (xk-coordinate)
        dx = float(xs) / Nx
        xu = np.linspace(0, xs, Nx+1)
        xeta = np.zeros(Nx+2)
        xeta[0:-1] = np.linspace(-0.5*dx, xs-0.5*dx, Nx+1)
        xeta[-1] = xs

        for i in range(len(bpu)):
            bpu[i] = slope(xu[i])
        for i in range(len(bpeta)):
            bpeta[i] = slope(xeta[i])

    bp_out = open('bp', 'w')
    for i in range(len(bpu)):
        bp_out.write('%g\t%g\n' % (xu[i], bpu[i]))

    return xu, xeta, dx, bpu, bpeta

    
def time_grid(time_steps, dt):
    ts = dt*time_steps
    
    t = np.linspace(0, ts, time_steps+1)
    
    return t


def arrays(xu, xeta, amplitude, dt):
    u0 = np.zeros(len(xu))
    u1 = np.zeros(len(xu))
    eta0 = np.zeros(len(xeta))
    eta1 = np.zeros(len(xeta))

    u0[:] = gaussian(xu[:], 0, amplitude, 'u')
    eta0[1:-1] = gaussian(xeta[1:-1], 0, amplitude, 'eta')
    eta0[0], eta0[-1] = eta0[1], eta0[-2]
    output(u0, xu, eta0, xeta, 0)

    return u0, u1, eta0, eta1


def gaussian(x, t, amplitude, out):

    start = 0
    var1 = x - t
    var2 = x + t

    if out == 'eta':
        g = amplitude * ( np.exp(-(var1 - start)**2) + np.exp(-(var2 - start)**2) )
    elif out == 'u':
        g = amplitude * ( np.exp(-(var1 - start)**2) - np.exp(-(var2 - start)**2) )

    return g


###################################################################
###                 CALCULATION                                 ###
###################################################################


def ghost_points(eta0, eta1, xeta, dx, dt, alfa):
    #eta0[-1] = eta1[-1] - dt/dx*(eta1[-1]- eta1[-2])

    if alfa == 0.0:
        eta1[0] = eta1[1]
        eta1[-1] = eta1[-2] 
    else:
        eta1[0] = eta1[1]
        eta1[-1] = eta1[-4] * (xeta[-1]-xeta[-3]) / (xeta[-4]-xeta[-3]) * (xeta[-1]-xeta[-2]) / (xeta[-4]-xeta[-2]) + \
                    eta1[-3] * (xeta[-1]-xeta[-4]) / (xeta[-3]-xeta[-4]) * (xeta[-1]-xeta[-2]) / (xeta[-3]-xeta[-2]) + \
                    eta1[-2] * (xeta[-1]-xeta[-4]) / (xeta[-2]-xeta[-4]) * (xeta[-1]-xeta[-3]) / (xeta[-2]-xeta[-3])
    return eta1[0], eta1[-1]



###########################################################
###                 FORMAT                              ###
###########################################################



def output(u1, xu, eta1, xeta, time_step):

    uout = open('u/u%g'%time_step, 'w')
    for i in range(len(u1)):
        uout.write('%g\t%g\n' % (xu[i], u1[i]))
    uout.close()

    etaout = open('eta/eta%g'%time_step, 'w')
    for i in range(len(eta1)):
        etaout.write('%g\t%g\n' % (xeta[i], eta1[i]))
    etaout.close()

    
    
def error(u1, xu, eta1, xeta, amplitude, time_step, dt):

    error_u = np.zeros(len(u1))
    error_u[:] = u1[:] - gaussian(xu[:], (time_step)*dt, amplitude, 'u')
    error_u_out = open('error_u/error_u%g'%time_step, 'w')
    for i in range(len(error_u)):
        error_u_out.write('%g\t%g\n' % (xu[i], error_u[i]))
    error_u_out.close()

    error_eta = np.zeros(len(eta1))
    error_eta[:] = eta1[:] - gaussian(xeta[:], (time_step+0.5)*dt, amplitude, 'eta')
    error_eta_out = open('error_eta/error_eta%g'%time_step, 'w')
    for i in range(len(error_eta)):
        error_eta_out.write('%g\t%g\n' % (xeta[i], error_eta[i]))
    error_eta_out.close()
    #print max(error_u)
    #print max(error_eta)



def plotting(xu, xeta, u1, eta1, bpu, bpeta):

    sci.figure(1)
    sci.plot(xeta, eta1, xeta, -bpeta, axes=(xu[0], xu[-1], -1.1*max(bpeta), 2.0*max(bpeta)))#, xeta, gaussian(xeta, (time_step+0.5)*dt, 'eta'))
    sci.figure(2)
    sci.plot(xu, u1, xu, -bpu, axes=(xu[0], xu[-1], -1.1*max(bpu), 2.0*max(bpu)))#, xu, gaussian(xu, (time_step)*dt, 'u'))


def check_folder():
    if os.path.exists('eta/'):
        shutil.rmtree('eta/')
    os.makedirs('eta/')
    if os.path.exists('u/'):
        shutil.rmtree('u/')
    os.makedirs('u/')
    if os.path.exists('error_u/'):
        shutil.rmtree('error_u/')
    os.makedirs('error_u/')
    if os.path.exists('error_eta/'):
        shutil.rmtree('error_eta')
    os.makedirs('error_eta/')


def read_file(file_name):
    parameters = []
    infile = open(file_name, 'r')
    for line in infile:
        variables = line.strip().split()
        parameters.append(variables[1])
    infile.close()
    return parameters[:]


############################################################################################################


def main():
    try:
        input_data = sys.argv[1]
        eq, profile, depth, alfa, xk, Nx, amplitude_factor, dt, ratio, time_steps = read_file(input_data)
    except IndexError:
        xk = 5
        alfa = 45
        time_steps = 100
        dt = 0.1
        Nx = 60
        depth = 1
        amplitude_factor = 0.2
        eq = 1
        ratio = 1
        profile = 1

    eq, profile, Nx, time_steps = int(eq), int(profile), int(Nx), int(time_steps)
    depth, alfa, xk, dt, ratio = float(depth), float(alfa), float(xk), float(dt), float(ratio)

    amplitude = float(amplitude_factor) * depth


    check_folder()
    xu, xeta, dx, bpu, bpeta = domain(xk, alfa, Nx, depth)
    t = time_grid(time_steps, dx)
    u0, u1, eta0, eta1 = arrays(xu, xeta, amplitude, dt)
    nu = len(u0)
    neta = len(eta0)-1

    for time_step in range(1, len(t)):

        u1[0:nu] = u0[0:nu] - dt/dx * (eta0[1:nu+1] - eta0[0:nu])
        eta1[1:neta] = eta0[1:neta] - dt/dx * (bpu[1:neta]*u1[1:neta] - bpu[0:neta-1]*u1[0:neta-1])

        eta1[0], eta1[-1] = ghost_points(eta0, eta1, xeta, dx, dt, alfa)
  

        eta0 = eta1
        u0 = u1
        

     
        output(u1, xu, eta1, xeta, time_step)
        error(u1, xu, eta1, xeta, amplitude, time_step, dt)
        #plotting(xu, xeta, u1, eta1, bpu, bpeta)


main()

