import numpy as np
import scitools.std as sci
import scipy.special as sp
import sys
import matplotlib.pyplot as plt
import os
import shutil


"""
    parameters:
    x0:     start of grid, Defalult=0.0
    xs:     end of grid
    xk:     point where inclination starts
    xb:     point where bottom profile breach water surface
    xu:     array for velocity
    xeta:   array for surface elevation
    Nx:     number of grid points
    dx:     distance between grid points
    t0:     start time, Default=0.0
    T:      end time
    dt:     time steps
    alfa:   angle of inclination
    depth:  initial depth of water
    bpu:    array with depth of bottom profile, points matches xu
    bpeta:  array with depth of bottom profile, points matches xeta 
    C:      dt / dx

    
    functions:
    main:           main loop
    parameters:     takes parameters from terminal
                    returns parameters; eq, depth, profile_par, alfa, 
                                        Nx, initial_condition, end_time, output

    calculation:    does the calculation
                    takes eq, initial_condition, u0, u1, eta0, eta1, xu, xeta, bpu, bpeta, t, C, output    
    
    grid:           creates grid
                    takes xs, Nx
                    returns x0, xu, xeta, dx
    bottom_profile: creates bottomprofile
                    takes profile_par, Nx, depth, alfa
                    returns xu, xeta, bpu, bpeta, dx, x0

    time_grid:      creates time grid
                    takes end_time, dx
                    returns t, C

    arrays:         creates arrays for velocity and surface elevation with initial cond
                    takes xu, xeta, initial_condition
                    returns u0, u1, eta0, eta1
    
    output_format:  what does the user from the program
                    takes output
                    returns out_format

"""



def main():

    eq, depth, profile_par, alfa, Nx, initial_condition, end_time, output = parameters()

    xu, xeta, bpu, bpeta, dx, x0 = bottom_profile(profile_par, Nx, depth, alfa)

    t, C = time_grid(end_time, dx)

    u0, u1, eta0, eta1 = arrays(xu, xeta, initial_condition)

    calculation(eq, initial_condition, u0, u1, eta0, eta1, xu, xeta, bpu, bpeta, t, C, output)



def calculation(eq, 
                initial_condition, 
                u0, 
                u1, 
                eta0, 
                eta1,
                xu,
                xeta, 
                bpu,
                bpeta, 
                t, 
                C, 
                output):
    """
        Performs the calculation
    """

    #if initial_condition[0] == 1:
    #    for j in range(len(eta1)):
    #        eta1[j] = eta0[j] - C * (bpu[j+1]*u0[j+1] - bpu[j]*u0[j])

    for time_step in range(1, len(t)+1):
        #print 'progress: ', float(time_step) / len(t) * 100, '%'

        #if eq == 'LSW':
        if initial_condition[0] == 2:
            # starting with a surface elevation
            u0[0] = u0[1]
            u0[-1] = u0[-2]

            for i in range(1, len(u0)-1):
                u1[i] = u0[i] - C * (eta0[i] - eta0[i-1])

            for j in range(0, len(eta0)):
                eta1[j] = eta0[j] - C * (bpu[j+1]*u1[j+1] - bpu[j]*u1[j])
        
            eta0 = eta1
            u0 = u1
            
        
        elif initial_condition[0] == 1:
            # starting with a velocity condition

            if time_step < 30:
                u0[0] = time_dependent(initial_condition[1], time_step)
            else:
                u0[0] = u0[1]
                u0[-1] = u0[-2]

            for j in range(0, len(eta1)):
                eta1[j] = eta0[j] - C * (bpu[j+1]*u0[j+1] - bpu[j]*u0[j])

            for i in range(1, len(u1)-1):
                u1[i] = u0[i] - C * (eta1[i] - eta1[i-1])

            #u1[0] = u1[1]
            #u1[-1] = u1[-2]

            eta0[:] = eta1
            u0[:] = u1


        output_format(output, eta1, xeta, bpu, bpeta, time_step)




def parameters():
    """
        Let users choose parameters from terminal or read from file
        Returns parameters
    """

    print 'Read input from file?'
    read = int(raw_input('No(1), Yes(2), Init_vel(3), Init_surf(4)'))
    if read == 1:
        pass
    elif read == 2:
        file_name = str(raw_input('File name: '))
        while os.path.isfile(file_name) == False:
            print 'Not a valid file name'
            file_name = str(raw_input('File name: '))
        infile = open(file_name, 'r')
        for line in infile:
            print line
    elif read == 3:
        eq, depth, profile_par, alfa, Nx, initial_condition, end_time, output = 1, 2, [2, 0.1], 88, 400, [1, 2], 10, [2,1]  
        if output[1] == 1:
            if os.path.exists('eta/'):
                shutil.rmtree('eta/')
            os.makedirs('eta/')
        elif output[1] == 2:
            if os.path.exists('u/'):
                shutil.rmtree('u/')
            os.makedirs('u/')
        elif output[1] == 3:
            if os.path.exists('eta/'):
                shutil.rmtree('eta/')
            os.makedirs('eta/')
            if os.path.exists('u/'):
                shutil.rmtree('u/')
            os.makedirs('u/')      
        return eq, depth, profile_par, alfa, Nx, initial_condition, end_time, output       
    elif read == 4:
        eq, depth, profile_par, alfa, Nx, initial_condition, end_time, output = 1, 1, [2, 1], 45, 100, [2, 0.3], 10, [1]
        return eq, depth, profile_par, alfa, Nx, initial_condition, end_time, output
        
    print ''

    # equation type
    print 'Choose equation solver'
    print 'LSW(1)'
    try:
        eq = int(raw_input('equation: '))
    except TypeError:
        print 'Not a valid entry, number must be 1'
        eq = int(raw_input('equation: '))
    if eq == 1:
        eq = 'LSW'
    print ''

    # depth
    print 'Choose depth of bottom profile'
    try:
        depth = abs(float(raw_input('depth: ')))
    except TypeError:
        print 'Depth must be a positive integer'
        depth = abs(float(raw_input('depth: ')))

    # type of bottom profile
    profile_par = []
    print 'Choose bottom profile'
    print 'Constant(1), Linear(2), Read from file(3)'
    try:
        profile = int(raw_input('profile: '))
    except TypeError:
        print 'Not a valid entry, number must be between 1-3'
        profile = int(raw_input('profile: '))
    profile_par.append(profile)
    if profile == 1:
        print 'Enter length of pool'
        xs = float(raw_input('length: '))
        profile_par.append(xs)
    elif profile == 2:
        print 'Enter point where inclination starts'
        xk = float(raw_input('starting point of inclination: '))
        profile_par.append(xk)
    elif profile == 3:
        print 'option not availble yet'
        return 0
    print ''

    # angle of inclination
    if profile == 2:
        print 'Choose the angle of inclination'
        try:
            alfa = float(raw_input('degree: '))
            if alfa >= 90.0 or alfa < 0.0:
                print 'The degree must be in the range of 0-90'
                alfa = float(raw_input('degree: '))
        except TypeError:
            print 'Angle must be in degrees between 0-90'
            alfa = float(raw_input('degree: '))
        print ''
    else:
        alfa = 0.0

    # grid points
    print 'Choose number of spatial grid points'
    try:
        Nx = int(raw_input('Number of grid points: '))
    except TypeError:
        print 'Number of spatial grid points must be an interger > 0'
        Nx = int(raw_input('Number of grid points: '))
    print ''
    
    # initial condition
    initial_condition = []
    print 'Choose initial condition'
    print 'velocity condition(1) or surface condition(2)'
    try:    
        IC = int(raw_input('Initial condition: '))
    except TypeError:
        print 'Initial condition must be 1 or 2'
        IC = int(raw_input('initial condition: '))
    initial_condition.append(IC)
    if IC == 1:
        print 'Choose inital velocity'
        initial_velocity = float(raw_input('Initial velocity: '))
        initial_condition.append(initial_velocity)
    elif IC == 2:
        print 'Choose initial surface elevation'
        initial_elevation = float(raw_input('Initial surface elevation: '))
        initial_condition.append(initial_elevation)
    print ''

    # time
    print 'Choose time length'
    try:
        end_time = float(raw_input('time: '))
    except TypeError:
        print 'Time must be a positive integer or float number'
        end_time = float(raw_input('time: '))
    print ''
    
    # output
    output = []
    print 'Choose output'
    print 'Live plot(1), Save to file(2)'
    try:
        out1 = int(raw_input('Output: '))
    except TypeError:
        print 'Output must be chosen as an integer 1-2'
        out1 = int(raw_input('Output: '))
    output.append(out1)
    if out1 == 2:
        print 'Surface elevation(1), wave velocity(2), both(3)'
        out2 = int(raw_input('Data to collect: '))
        output.append(out2)
        if out2 == 1:
            if os.path.exists('eta/'):
                shutil.rmtree('eta/')
            os.makedirs('eta/')
        elif out2 == 2:
            if os.path.exists('u/'):
                shutil.rmtree('u/')
            os.makedirs('u/')
        elif out2 == 3:
            if os.path.exists('eta/'):
                shutil.rmtree('eta/')
            os.makedirs('eta/')
            if os.path.exists('u/'):
                shutil.rmtree('u/')
            os.makedirs('u/')

        
    

    return eq, depth, profile_par, alfa, Nx, initial_condition, end_time, output



def grid(xs, 
         Nx):
    """
        creates the grid with staring point x0 and end
        point xs with one ghost point on each side. 
        Nx is number of points and dx is the
        distance between two points
    """
        
    x0 = 0.0
    dx = (xs - x0) / (Nx)
    xu = np.linspace(x0-dx/2, xs+dx/2, Nx+2)
    xeta = np.linspace(x0, xs, Nx+1)


    return x0, xu, xeta, dx



def bottom_profile(profile_par,
                   Nx,
                   depth,
                   alfa):
    """
        creates a bottom profile 
        must fix for alfa=90
    """

    if profile_par[0] == 1:
        xs = profile_par[1]
        x0, xu, xeta, dx = grid(xs, Nx)
        bpu = depth*np.ones(len(xu))
        bpeta = depth*np.ones(len(xeta))

    elif profile_par[0] == 2:

        xk = profile_par[1]
        degree_incline = alfa / (90.0 - alfa)
        xs = depth / degree_incline + xk    
        x0, xu, xeta, dx = grid(xs, Nx)

        def slope(coordinate):
            if coordinate <= xk:
                return depth
            else:
                return depth + degree_incline * (xk-coordinate)

        bpu = np.zeros(len(xu))
        for i in range(len(bpu)):
            bpu[i] = slope(xu[i])
        bpeta = np.zeros(len(xeta))
        for i in range(len(bpeta)):
            bpeta[i] = slope(xeta[i])


    return xu, xeta, bpu, bpeta, dx, x0


def slope(coordinate):
    if coordinate <= start:
        return depth
    else:
        return depth + degree_incline * (start-coordinate)


def time_grid(end_time, dx):
    """
        creates a "time grid"
    """
    
    t0 = 0.0
    #dt = 0.05 * (end_time - t0) * dx
    #Nt = (end_time - t0) / (dt)
    Nt = (end_time - t0) / (0.5*dx)
    t = np.linspace(t0, end_time, Nt+1)
    dt = t[1]-t[0]

    C = dt / dx

    return t, C


def arrays(xu, xeta, initial_condition):
    """
        returns all arrays with initial condition
    """
    
    u0 = np.zeros(len(xu))
    u1 = np.zeros(len(xu))
    eta0 = np.zeros(len(xeta))
    eta1 = np.zeros(len(xeta))

    if initial_condition[0] == 1:
    #    u0[0] = initial_condition[1]
        u0[0] = 0
    
    elif initial_condition[0] == 2:
        for i in range(len(eta0)):
            eta0[i] = gaussian(initial_condition[1], xeta[i])

    return u0, u1, eta0, eta1


def gaussian(amplitude, x, time_step):
    f = amplitude * np.exp(-(x**2) / (2*0.1**2))
    if f < 0.0001:
        f = 0
    return f


def time_dependent(velocity, time_step):
    if time_step > 10 and time_step < 20:
        return 0.01*time_step*velocity
    else:
        return 0
        


def output_format(output, eta1, xeta, bpu, bpeta, time_step, max_runup=None):
    """
        What to return. 
        Data or plots
    """

    if output[0] == 1:
        plotting(eta1, xeta, bpeta)
    elif output[0] == 2:
        if output[1] == 1:
            write_eta(eta1, bpeta, xeta, time_step)
        elif output[1] == 2:
            write_u(u1, bpu, xu, time_step)
        elif output[1] == 3:
            write_eta(eta1, bpeta, xeta, time_step)
            write_u(u1, bpu, xu, time_step)

            
        
def plotting(eta1, xeta, bpeta):
    sci.plot(xeta, eta1, xeta, -bpeta, axis=(xeta[0], xeta[-1], -1.1*max(bpeta), 1.1*max(bpeta)))


def write_eta(eta1, bpeta, xeta, time_step):
    """
        stores data of surface elevation and bottom profile in a .txt file
    """
    outfile = open('eta/eta%g.txt'%time_step, 'w')
    for value in range(len(eta1)):
        outfile.write('%g\t%g\t%g\n' % (eta1[value], xeta[value], bpeta[value]))
    outfile.close()

def write_u(u1, bpu, xu, time_step):
    """
        stores data of surface elevation and bottom profile in a .txt file
    """
    outfile = open('u/u%g.txt'%time_step, 'w')
    for value in range(len(u1)):
        outfile.write('%g\t%g\t%g\n' % (u1[value], xu[value], bu[value]))
    outfile.close()
    
        

"""
def save_plot(eta1, xeta, bp, time_step):
    if not os.path.exists("figures/"):
        os.makedirs("figures/")
    #plt.ioff()
    #fig = plt.figure()
    plt.axis([xeta[0], xeta[-1], -1.1*max(bp), 1.1*max(bp)])
    plt.plot(xeta, eta1, xeta, -bp[1:])
    plt.savefig('figures/eta%g.png' % time_step)
    #plt.close()
"""

def convergence_test(profile_par,
         Nx,
         depth,
         alfa):

    J = lambda v, z: sp.jv(v, z)
    xu, xeta, bpu, bpeta, dx, x0 = bottom_profile(profile_par, Nx, depth, alfa)
    Jv = J(0, bpeta)
    plt.plot(xeta, Jv)
    plt.show()
    raw_input()


       


                
                           



if __name__ == '__main__':
    main()
    #convergence_test(profile_par=[2,1], Nx=100, depth=2, alfa=45)
    #runup_height([1, 2, 3, 4], 2)
    #write_eta([0,1,2,3,4,5],2, 3)









