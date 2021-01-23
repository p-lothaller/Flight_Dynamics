# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:25:48 2020

@author: vladg, phillipe
@version: 2.2
"""

import numpy as np
from control.matlab import *
import matplotlib.pyplot as plt
import numpy.linalg


# +++++++++++++++++++++++++++++++++ Helper Functions ++++++++++++++++++++++++++++++++++++++++++++++

def plotting(x,y,name,variable,unit):
    """Use this for plotting."""
    ax = plt.figure(str(name))
    # ax.legend("best")

    plt.plot(x,y,label=name)
    plt.xlabel("t [s]")

    lab = str(str(variable)+" "+"["+unit+"]")
    plt.ylabel(lab)
    plt.grid(True)
    plt.show()


#+++++++++++++++++++++++++++++++++++ Global variables+++++++++++++++++++++++++++++++++++++++++++++++

# Citation 550 - Linear simulation
    # Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * np.pi / 180   # stabiliser angle of incidence [rad]

oew = 4157.174              #Operational Empty Weight [kg]
m_payload = 765             # Payload mass [kg]
    # Aerodynamic properties
e      = 0.5087   #0.7334          # Oswald factor [ ]
CD0    = 0.00227  #0.02180            # Zero lift drag coefficient [ ]
CLa    = 4.811    #4.383            # Slope of CL-alpha curve [ ]

    # Constant values concerning atmosphere and gravity
rho0   = 1.2250          # air density at sea level [kg/m^3]
lam    = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15       # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)


    #Simulation parameters:
nsteps = 10**3

#+++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def main(t0,deltat,input_type):
    """Input type: elevator
                    rudder
                    airleron"""

    #C.G location
    xcg = 0.25 * c

    # Stationary flight condition
    m_fuel = 1197.484           # CHANGE M_fule(t0) Total fuel mass [kg]
    gamma  = 0                  # !!!!!! flight path angle -
    hp0    = 1527.048      	    # CHANGE pressure altitude in the stationary flight condition [m]
    V0     = 127.067            # CHANGE true airspeed in the stationary flight condition [m/sec]
    alpha0 = np.radians(1.4)    # CHANGE angle of attack in the stationary flight condition [rad]
    th0    = alpha0 + gamma     # CHANGE pitch angle in the stationary flight condition [rad]

    # Aircraft mass
    m      =  4989.516 + m_payload         # CHANGE mass [kg]

    # Longitudinal stability
    Cma    = -0.4433 # -0.4934        # CHANGE longitudinal stabilty [ ]
    Cmde   = -1.001 #-1.031           # CHANGE elevator effectiveness [ ]

    # air density [kg/m^3]
    rho    = rho0 * pow( ((1+(lam * hp0 / Temp0))), (-((g / (lam*R)) + 1)))
    W      = m * g            # [N]       (aircraft weight)

    # Aircraft inertia (depend on t0):
    muc    = m / (rho * S * c) #CHANGE
    mub    = m / (rho * S * b) #CHANGE
    KX2    = 0.019
    KZ2    = 0.042
    KXZ    = 0.002
    KY2    = 1.25 * 1.114

    # Aerodynamic constants:

    Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
    CNwa   = CLa                    # Wing normal force slope [1/rad]
    CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
    depsda = 4 / (A + 2)            # Downwash gradient [ ]

    # Lift and drag coefficient (depend on t0):

    CL = 2 * W / (rho * V0 ** 2 * S)     # CHANGES Lift coefficient [ ]
    CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]

    # Stabiblity derivatives
    CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu    = -0.095         #corrected
    CXa    = +0.47966		# Positive! (has been erroneously negative since 1993)
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728

    CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = -0.37616
    CZa    = -5.74340
    CZadot = -0.00350
    CZq    = -5.66290
    CZde   = -0.69612

    Cmu    = +0.06990   #positive!
    Cmadot = +0.17800   #positive!
    Cmq    = -8.79415

    CYb    = -0.7500
    CYbdot =  0
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300

    Clb    = -0.10260
    Clp    = -0.71085
    Clr    = +0.23760
    Clda   = -0.23088
    Cldr   = +0.03440

    Cnb    =  +0.1348
    Cnbdot =   0
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    Cnda   =  -0.0120
    Cndr   =  -0.0939

    #c-matrix dimensions
    s1 = (4,4)
    s2 = (4,1)
    s3 = (4,2)

    #Creating the different c-matrices (c1, c2 &c3) for symmetrical flight
    #c1 matrix
    c1 = np.zeros(s1)
    c1[0,0] = -2*muc*(c/V0)
    c1[1,1] = (CZadot - 2*muc)*(c/V0)
    c1[2,2] = -(c/V0)
    c1[3,1] = Cmadot*(c/V0)
    c1[3,3] = -2*muc*KY2*((c/V0)**2)

    #c2 matrix
    c2 = np.zeros(s1)
    c2[0,0] = -CXu
    c2[0,1] = -CXa
    c2[0,2] = -CZ0
    c2[0,3] = -CXq
    c2[1,0] = -CZu
    c2[1,1] = -CZa
    c2[1,2] = -CX0
    c2[1,3] = -(CZq + 2*muc)*(c/V0)
    c2[2,3] = -(c/V0)
    c2[3,0] = -Cmu
    c2[3,1] = -Cma
    c2[3,3] = -Cmq*(c/V0)

    #c3 matrix
    c3 = np.zeros(s2)
    c3[0,0] = -CXde
    c3[1,0] = -CZde
    c3[3,0] = -Cmde


    #Creating the different c-matrices (c4, c5 &c6) for asymmetrical flight

    #c4 matrix
    c4 = np.zeros(s1)
    c4[0,0] = (CYbdot - 2*mub)*(b/V0)
    c4[1,1] = (-0.5)*(b/V0)
    c4[2,2] = -4*mub*KX2*(b/V0)*(b/(2*V0))
    c4[2,3] = 4*mub*KXZ*(b/V0)*(b/(2*V0))
    c4[3,0] = Cnb*(b/V0)
    c4[3,2] = 4*mub*KXZ*(b/V0)*(b/(2*V0))
    c4[3,3] = -4*mub*KZ2*(b/V0)*(b/(2*V0))

    #c5 matrix
    c5 = np.zeros(s1)
    c5[0,0] = CYb
    c5[0,1] = CL
    c5[0,2] = CYp*(b/(2*V0))
    c5[0,3] = (CYr - 4*mub)*(b/(2*V0))
    c5[1,2] = (b/(2*V0))
    c5[2,0] = Clb
    c5[2,2] = Clp*(b/(2*V0))
    c5[2,3] = Clr*(b/(2*V0))
    c5[3,0] = Cnb
    c5[3,2] = Cnp*(b/(2*V0))
    c5[3,3] = Cnr*(b/(2*V0))

    #c6 matrix
    c6 = np.zeros(s3)
    c6[0,0] = -CYda
    c6[0,1] = -CYdr
    c6[2,0] = -Clda
    c6[2,1] = -Cldr
    c6[3,0] = -Cnda
    c6[3,1] = -Cndr

    # Time responses for unit steps:
    t = np.linspace(t0, deltat, nsteps)

    #Now, we distinct between inputs:

    if input_type=="elevator":
        #Symmetric system is triggered:

        #Creating the state matrix(A) and the input matrix(B) for symmetrical flight - xdot = c1^-1*c2*x c1^-1*c3*u = Ax + Bu
        A_s = np.dot(np.linalg.inv(c1), c2)
        B_s = np.dot(np.linalg.inv(c1), c3)
        C_s = np.identity(4)
        D_s = np.zeros((4, 1))

        #System in state-space
        sys_s = StateSpace(A_s, B_s, C_s, D_s)
        poles_s = pole(sys_s)
        damp(sys_s)
        #print("Eigen values of the symmetric system: ", poles_s) #verified

        # Time responses for unit steps:

        # u = input_array
        u = np.ones(t.shape)
        yout,t,u = lsim(sys_s,u,t)   #general time response

        u_out_s =     yout[:,0]
        alpha_out_s = yout[:,1]
        theta_out_s = yout[:,2]
        q_out_s =     yout[:,3]

        #Plotting....
        plotting(t,u_out_s,str("u Response for " +input_type+ " input"),r"$u$","m/s")
        plotting(t,alpha_out_s,str("Alpha Response for " +input_type+ " input"),r"$\alpha$","-")
        plotting(t,theta_out_s,str("Theta Response for " +input_type+ " input"),r"$\theta$","-")
        plotting(t,q_out_s,str("q Response for " +input_type+ " input"),"$q$",r"1/s")


    else:
        #Creating the state matrix(A) and the input matrix(B) for asymmetrical flight - y = c4^-1*c5*x c4^-1*c5*u = Ax + Bu
        A_a = -np.dot(np.linalg.inv(c4), c5)
        B_a = np.dot(np.linalg.inv(c4), c6)
        C_a = np.identity(4)
        #D_a depends on the input

        if input_type =="rudder":
            D_a = np.zeros((4, 2))
            D_a[:,0] = 1   #we should check this...
            u = np.ones((len(t),2)) * -0.804 #step input
            u[:,0]=1
            print(u.shape)

        elif input_type=="aileron":
            D_a = np.zeros((4, 2))
            D_a[:,1] = 1
            u = np.ones((len(t),2)) #step input
            u[:,1]=1

        #System in state-space
        sys_a = StateSpace(A_a, B_a, C_a, D_a)
        poles_a = pole(sys_a)
        #print("Eigen values of the asymmetric system: ", poles_a) #verified



        yout,t,u = lsim(sys_a,u,t)   #general time response for the input u

        u_out_a =     yout[:,0]
        alpha_out_a = yout[:,1]
        theta_out_a = yout[:,2]
        q_out_a =     yout[:,3]

        #Plotting...
        plotting(t,u_out_a,str("Beta Response for " + input_type +" input"), r"$beta$","-")
        plotting(t,alpha_out_a,str("Phi Response for " +input_type + " input"), r"$\phi$","-")
        plotting(t,theta_out_a,str("p Response for " +input_type + " input") , r"$p$" ,"1/s")
        plotting(t,q_out_a,str("r Response for " +input_type + " input"),  "$r$" ,r"1/s")

    return 1

if __name__=="__main__":

    main(0,140,"elevator")
