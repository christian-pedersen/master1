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
    IC:     Initial condition (surface elevation or velocity)
    IV:     Initial value

    
    functions:
    main:           main loop
    parameters:     takes parameters from terminal
                    returns parameters; eq, depth, profile_par, alfa, 
                                        Nx, IC, IV, end_time, output

    calculation:    does the calculation
                    takes eq, IC, IV, u0, u1, eta0, eta1, xu, xeta, bpu, bpeta, t, C, output    
    
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
                    takes xu, xeta, IC, IV
                    returns u0, u1, eta0, eta1
    
    output_format:  what does the user from the program
                    takes output
                    returns out_format


    IF READ FROM FILE:
    print data in input file with one variable per line
    eq=                 LSW(1)
    depth=              a positive float number(e.g. 1.0)
    profile=            linear(1), read from file(2)
    xs=                 a positive float or integer
    alfa=               float number between 0.0<=alfa<90.0
    Nx=                 a positive integer
    IC=                 velocity(1), surface elevation(2)
    IV=                 a positive float
    end_time=           positive float
    ratio=              a positive float
    output=             live plot(1), save to file(2)
    if_save=            eta(1), u(2), both(3)   

"""



def main():

    try:
        input_data = sys.argv[1]
        check_folder('error')
        eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output = read_file(input_data)
    except IndexError:
        eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output = parameters()

    xu, xeta, bpu, bpeta, dx, x0 = bottom_profile(profile, xs, Nx, depth, alfa)

    t, dt, C = time_grid(end_time, ratio, dx)

    u0, u1, eta0, eta1 = arrays(xu, xeta, IC, IV)

    calculation(eq, IC, IV, u0, u1, eta0, eta1, xu, xeta, xs, bpu, bpeta, alfa, t, dt, C, ratio, output)



def calculation(eq, 
                IC,
                IV, 
                u0, 
                u1, 
                eta0, 
                eta1,
                xu,
                xeta,
                xs, 
                bpu,
                bpeta,
                alfa, 
                t, 
                dt,
                C,
                ratio, 
                output):
    """
        Performs the calculation
    """

    for time_step in range(1, len(t)+1):

        if IC == 1:
            # starting with a velocity condition

            u0[0], u0[-1] = ghost_points(u0, u1, C, alfa, IC, xu[0], xu[-1], IV, time_step, ratio, dt)

            for i in range(1, len(u1)-1):
                u1[i] = u0[i] - C * (eta0[i] - eta0[i-1])

            for j in range(0, len(eta1)):
                eta1[j] = eta0[j] - C * (bpu[j+1]*u1[j+1] - bpu[j]*u1[j])

            #u0[0], u0[-1] = ghost_points(u0, u1, C, alfa, IC, xu[0], xu[-1], IV, time_step, ratio, dt)

            eta0 = eta1
            u0 = u1


        elif IC == 2:
            # starting with a surface elevation
            u0[0], u0[-1] = ghost_points(u0, u1, C, alfa, IC)

            for i in range(1, len(u0)-1):
                u1[i] = u0[i] - C * (eta0[i] - eta0[i-1])

            for j in range(0, len(eta0)):
                eta1[j] = eta0[j] - C * (bpu[j+1]*u1[j+1] - bpu[j]*u1[j])

            eta0 = eta1
            u0 = u1

    
        output_format(output, u1, eta1, xu, xeta, bpu, bpeta, time_step)
        flat_error(eta1, u1, IV, bpeta, xeta, bpu, xu, xs, time_step, dt, ratio, IC)

def ghost_points(u0, u1, C, alfa, IC, x_start=None, x_end=None, IV=None, time_step=None, ratio=None, dt=None):

    if alfa == 0.0 and IC == 2:
        u0[0] = u1[0] + C*(u1[1] - u1[0])
        u0[-1] = u1[-1] - C*(u1[-1]- u1[-2])
    elif alfa == 0.0 and IC == 1:
        u0[0] = time_dependent(IV, x_start, x_end, time_step, ratio, dt)
        if u0[0] == False:
            u0[0] = u1[0] + C*(u1[1] - u1[0])
        u0[-1] = u1[-1] - C*(u1[-1]- u1[-2])
    else:
        u0[0] = u1[0] + C*(u1[1] - u1[0])
        u0[-1] = u0[-2]
    return u0[0], u0[-1]

####################################################################################
###                                PARAMETERS                                    ###
####################################################################################


def parameters():
    """
        Let users choose parameters from terminal or read from file
        Returns parameters
    """

    print 'Read input from file?'
    read = int(raw_input('No(1), Yes(2), Init_vel(3), Init_surf(4)'))
    while read not in (1,2,3,4):
        print 'User must choose an integer between 1-4'
        read = int(raw_input('No(1), Yes(2), Init_vel(3), Init_surf(4)'))
    if read == 1:
        pass
    elif read == 2:
        file_name = str(raw_input('File name: '))
        while os.path.isfile(file_name) == False:
            print 'Not a valid file name'
            file_name = str(raw_input('File name: '))
        par = read_file(file_name)
        eq, depth, profile_par, alfa, Nx, IC, IV, end_time, output = par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]
        return eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, output        

    elif read == 3:
        eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output = 1, 1, 1, 10, 0, 50, 1, 1, 20, 1, [1]  
        if len(output) == 2:
            check_folder(output[1])
            check_folder('error')
        return eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output       
    elif read == 4:
        eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output = 1, 1, 1, 10, 0, 50, 2, 1, 10, 0.1, [1]
        if len(output) == 2:
            check_folder(output[1])
        check_folder('error')
        return eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output
        
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


    # type of bottom profile
    print 'Choose bottom profile'
    print 'Linear(1), Read from file(2)'
    profile = int(raw_input('profile: '))
    while profile not in (1,2):
        print 'Not a valid entry, number must be between 1-2'
        profile = int(raw_input('profile: '))
    print ''
    if profile == 1:
        # pool length
        print 'Enter length of pool'
        xs = float(raw_input('length: '))
        print ''

        # depth
        print 'Choose depth of bottom profile'
        try:
            depth = abs(float(raw_input('depth: ')))
        except TypeError:
            print 'Depth must be a positive integer'
            depth = abs(float(raw_input('depth: ')))
        print ''

        # angle of inclination
        print 'Choose the angle of inclination'
        alfa = float(raw_input('degree: '))
        while alfa >= 90.0 or alfa < 0.0:
            print 'The degree must be in the range of 0-90'
            alfa = float(raw_input('degree: '))
        print ''

        # grid points
        print 'Choose number of spatial grid points'
        Nx = int(raw_input('Number of grid points: '))
        print ''

    elif profile == 2:
        print 'option not availble yet'
        return 0
    

    # initial condition
    print 'Choose initial condition'
    print 'velocity condition(1) or surface condition(2)'
    try:    
        IC = int(raw_input('Initial condition: '))
    except TypeError:
        print 'Initial condition must be 1 or 2'
        IC = int(raw_input('initial condition: '))
    if IC == 1:
        print 'Choose inital velocity'
        IV = float(raw_input('Initial velocity: '))
    elif IC == 2:
        print 'Choose initial surface elevation'
        IV = float(raw_input('Initial surface elevation: '))
    print ''

    # time
    print 'Choose time length'
    try:
        end_time = float(raw_input('time: '))
    except TypeError:
        print 'Time must be a positive integer or float number'
        end_time = float(raw_input('time: '))
    print ''

    # time / space ratio
    print 'Choose dt / dx ratio, e.g. 0.5'
    try:
        ratio = float(raw_input('ratio: '))
    except TypeError:
        print 'ratio must be float number'   
        ratio = float(raw_input('ratio: ')) 
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
        check_folder(out2)


    return eq, depth, profile, xs, alfa, Nx, IC, IV, end_time, ratio, output





#################################################################################
###                         GRID POINTS                                       ###
#################################################################################


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



def bottom_profile(profile,
                   xs,
                   Nx,
                   depth,
                   alfa):
    """
        creates a bottom profile 
        must fix for alfa=90
    """

    x0, xu, xeta, dx = grid(xs, Nx)
    bpu = np.zeros(len(xu))
    bpeta = np.zeros(len(xeta))

    if alfa == 0.0:
        bpu = depth*np.ones(len(xu))
        bpeta = depth*np.ones(len(xeta))
    
    else:
        degree_incline = alfa / (90.0 - alfa)
        xk = xs - depth / degree_incline   

        def slope(coordinate):
            if coordinate <= xk:
                return depth
            else:
                return depth + degree_incline * (xk-coordinate)

        for i in range(len(bpu)):
            bpu[i] = slope(xu[i])
        for i in range(len(bpeta)):
            bpeta[i] = slope(xeta[i])


    return xu, xeta, bpu, bpeta, dx, x0


def time_grid(end_time, ratio, dx):
    """
        creates a "time grid"
    """
    
    t0 = 0.0
    #dt = 0.05 * (end_time - t0) * dx
    #Nt = (end_time - t0) / (dt)
    Nt = (end_time - t0) / (ratio * dx)
    t = np.linspace(t0, end_time, Nt+1)
    dt = t[1]-t[0]

    C = dt / dx

    return t, dt, C


def arrays(xu, xeta, IC, IV):
    """
        returns all arrays with initial condition
    """
    
    u0 = np.zeros(len(xu))
    u1 = np.zeros(len(xu))
    eta0 = np.zeros(len(xeta))
    eta1 = np.zeros(len(xeta))

    if IC == 1:
        u0[0] = np.exp(-(0.3*xeta[-1])**2)
    
    elif IC == 2:
        xs = xeta[-1]
        dx = xu[2] - xu[1]
        dt = dx
        for i in range(len(eta0)):
            eta0[i] = gaussian(IV, xeta[i], xs)
        for i in range(len(u0)):
            u0[i] = 0.5 * (np.exp(-(xu[i] - (0-0.5)*dt - 0.5*xs)**2) - np.exp(-(xu[i] + (0-0.5)*dt - 0.5*xs)**2))

    return u0, u1, eta0, eta1





###########################################################################################
###                             CONDITIONS                                              ###                                      
###########################################################################################


def gaussian(amplitude, x, xs, time_step=None, dt=None):
    variable = x 
    a = 0.5*xs
    f = amplitude * np.exp((-(variable-a)**2))
    return f


def time_dependent(velocity, x, xs, time_step, ratio, dt):
    variable = (time_step)*dt
    starting_point = -0.3*xs
    f = velocity * np.exp(-(-variable-starting_point)**2)  

    return f
    
        



###########################################################################################
###                             FORMATING                                               ###
###########################################################################################



def output_format(output, u1, eta1, xu, xeta, bpu, bpeta, time_step, max_runup=None):
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
    sci.figure(1)
    sci.plot(xeta, eta1, xeta, -bpeta, axis=(xeta[0], xeta[-1], -1.1*max(bpeta), 1.1*max(bpeta)))


def read_file(file_name):
    parameter = []
    infile = open(file_name)
    for line in infile:
        variables = line.strip().split()
        parameter.append(variables[1]) 
    if len(parameter) == 12:
        check_folder(int(parameter[11]))
        return int(parameter[0]), float(parameter[1]), int(parameter[2]), \
               float(parameter[3]), float(parameter[4]), int(parameter[5]), \
               int(parameter[6]), float(parameter[7]), float(parameter[8]), \
               float(parameter[9]), [int(parameter[10]), int(parameter[11])]
    elif len(parameter) == 11:
        return int(parameter[0]), float(parameter[1]), int(parameter[2]), \
               float(parameter[3]), float(parameter[4]), int(parameter[5]), \
               int(parameter[6]), float(parameter[7]), float(parameter[8]), \
               float(parameter[9]), int(parameter[10])
       

def check_folder(folder_name):
    if folder_name == 1:
        if os.path.exists('eta/'):
            shutil.rmtree('eta/')
        os.makedirs('eta/')
    elif folder_name == 2:
        if os.path.exists('u/'):
            shutil.rmtree('u/')
        os.makedirs('u/')
    elif folder_name == 3:
        if os.path.exists('eta/'):
            shutil.rmtree('eta/')
        os.makedirs('eta/')
        if os.path.exists('u/'):
            shutil.rmtree('u/')
        os.makedirs('u/')
    elif folder_name == 'error':
        if os.path.exists('error_u/'):
            shutil.rmtree('error_u/')
        os.makedirs('error_u/')
        if os.path.exists('error_eta/'):
            shutil.rmtree('error_eta')
        os.makedirs('error_eta/')
    


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
        outfile.write('%g\t%g\t%g\n' % (u1[value], xu[value], bpu[value]))
    outfile.close()

def write_error_eta(error_eta, xeta, bpeta, time_step):
    """
        stores data of error in a .txt file
    """
    outfile = open('error_eta/error_eta%g.txt'%time_step, 'w')
    for value in range(len(error_eta)):
        outfile.write('%g\t%g\t%g\n' % (error_eta[value], xeta[value], bpeta[value]))
    outfile.close()
    
def write_error_u(error_u, xu, bpu, time_step):
    outfile = open('error_u/error_u%g.txt'%time_step, 'w')
    for value in range(len(error_u)):
        outfile.write('%g\t%g\t%g\n' % (error_u[value], xu[value], bpu[value]))
    outfile.close()



##################################################################################################
###                                     TESTING                                                ###
##################################################################################################




def exact_gaussian(amplitude, depth, x, xs, time_step, dt, ratio, IC, out):
    c = np.sqrt(depth)
    variable1u = x - (time_step-0.5)*(dt) / c
    variable2u = x + (time_step-0.5)*(dt) / c
    variable1eta = x - (time_step)*(dt) / c
    variable2eta = x + (time_step)*(dt) / c
    if IC == 1:
        starting_point = -0.3*xs
        eta = amplitude * np.exp(-(variable1eta - starting_point)**2) 
        u = amplitude / depth * np.exp(-(variable1u - starting_point)**2)
        if out == 'eta':
            return eta
        elif out == 'u':
            return u

    elif IC == 2:
        starting_point = 0.5*xs
        eta =   0.5 * amplitude * ( np.exp(-(variable1eta - starting_point)**2) \
                                    + np.exp(-(variable2eta - starting_point)**2) )
        u = 0.5 * amplitude / c * ( np.exp(-(variable1u - starting_point)**2) \
                                    - np.exp(-(variable2u - starting_point)**2) )
        if out == 'eta':
            return eta
        elif out == 'u':
            return u


def flat_error(eta1,
               u1,
               amplitude,
               bpeta,
               xeta,
               bpu,
               xu,               
               xs,
               time_step,
               dt,
               ratio,
               IC):
    depth = bpeta[0]
    eta_exact = np.zeros(len(xeta))
    for i in range(len(xeta)):
        eta_exact[i] = exact_gaussian(amplitude, depth, xeta[i], xs, time_step, dt, ratio, IC, 'eta')
    #sci.figure(2)
    #sci.plot(xeta, eta_exact, axis=(0, xs, -1.1, 1.1))
    error_eta = [eta1[i] - eta_exact[i] for i in range(len(xeta))]
    write_error_eta(error_eta, xeta, bpeta, time_step)

    u_exact = np.zeros(len(xu))
    for i in range(len(xu)):
        u_exact[i] = exact_gaussian(amplitude, depth, xu[i], xs, time_step, dt, ratio, IC, 'u')
    error_u = [u1[i] - u_exact[i] for i in range(len(xu))]
    write_error_u(error_u, xu, bpu, time_step)
    #sci.figure(3)
    #sci.plot(xu, u_exact, axis=(0,xs,-1.1,1.1))
    


       


                
                           



if __name__ == '__main__':
    main()
    #convergence_test(profile=1, xs=6, Nx=100, depth=2, alfa=45)
    #runup_height([1, 2, 3, 4], 2)
    #write_eta([0,1,2,3,4,5],2, 3)









