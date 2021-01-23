# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:25:48 2020

@author: Group B07
Worked on plotting
@version: 3.5 22/03/2020
"""

import numpy as np
import control.matlab as cm
import math
import matplotlib.pyplot as plt
import numpy.linalg
from massbalance import *
from data_generator_complete import *
import warnings
import scipy.optimize as sp
warnings.filterwarnings("ignore")

plt.close('all')

# +++++++++++++++++++++++++++++++++ Helper Functions ++++++++++++++++++++++++++++++++++++++++++++++

#---------------------------------- fitting ------------------------------


def double(t,A1,eps1,omega_d1,A2,eps2,omega_d2):
    return 1.2 + A1 * np.exp(eps1*t) * np.sin(omega_d1 * t + 0) + A2 * np.exp(eps2*t) * np.cos(omega_d2 * t + 0)

def doubleSP(t,A1,eps1,omega_d1,A2,eps2,omega_d2,phi01,phi02):
     return   A1 * np.exp(eps1*t) * np.sin(omega_d1 * t + phi01) + A2 * np.exp(eps2*t) * np.cos(omega_d2 * t + phi02)

def simple(t,A1,eps1,omega_d1,phi01):
    return   A1 * np.exp(eps1*t) * np.sin(omega_d1 * t + phi01)

def short(t,A,eps,omega,B,C):
    return A*np.exp(t*eps)*np.cos(omega*t+b)

def calcT(eps):
    return -np.log(1/2)/abs(eps)

def calcP(omega):
    return 2*np.pi/abs(omega)
# ----------------------------------- plots ---------------------------------
def plotting(x,y,name,variable,unit,label_name="Simulation",title=None,mins=False):
    """Use this for plotting. It returns the figure so we can add more to it"""
    ax = plt.figure(str(name))
    # ax.legend("best")

    if mins:
        x/=60  #change time to mins from secs
        plt.xlabel("t [min]")
    else:
        plt.xlabel("t [s]")

    if title!= None:
        plt.title(str(title))

    if label_name=="Simulation":
        plt.plot(x-x[0],y,'--', label=label_name)
    else:
        plt.plot(x - x[0], y, label=label_name)


    lab = str(str(variable)+" "+"["+unit+"]")
    plt.legend(loc='best')
    plt.ylabel(lab)
    plt.grid(True)
    # plt.savefig(title)
    plt.show()

    return ax

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
e      = 0.5466             # Oswald factor [ ]
CD0    = 0.01980           # Zero lift drag coefficient [ ]
CLa    = 4.547            # Slope of CL-alpha curve [ ]

    # Constant values concerning atmosphere and gravity
rho0   = 1.2250          # air density at sea level [kg/m^3]
lam    = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15       # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

    #Loading in data

altlst = np.genfromtxt("Data_SI_correct/Dadc1_altSI.txt",skip_header=2)  #reading in the altitude values
taslst = np.genfromtxt("Data_SI_correct/Dadc1_tasSI.txt",skip_header=2)  #reading in the true airspeed values
aoalst = np.genfromtxt("Data_SI_correct/vane_AOASI.txt",skip_header=2)  #reading in the angle of attack
pitchlst = np.genfromtxt("Data_SI_correct/Ahrs1_bPitchRateSI.txt",skip_header=2)  #reading in the pitch
    #Simulation parameters:
nsteps = 10**3


#+++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++
def main(t0,deltat,t,input_type,input_u, CZa = -5.74340 , CZadot = -0.00350, CZq = -5.66290, Cmadot = 0.17800, Cmq = -8.79415, CXu=-0.095, CZu=-0.37616):
    """Input type: elevator
                    rudder
                    airleron"""


    #Find time
    idx = np.where(timelst == t0)[0]

    #Flight condition
    #  m_fuel    = 1197.484        # CHANGE Total fuel mass [kg]

    hp    = altlst[idx]      	    # CHANGE pressure altitude in the stationary flight condition [m]
    V     = taslst[idx]            # CHANGE true airspeed in the stationary flight condition [m/sec]
    alpha = np.radians(aoalst[idx])        # angle of attack in the stationary flight condition [rad]
    theta = np.radians(pitchlst[idx])    # pitch angle in the stationary flight condition [rad]
    gamma  = theta - alpha  # CHANGE flight path angle -

    # Aircraft mass
    m      =  mass(t0)        # mass [kg]

    # Longitudinal stability
    Cma    = -0.4435           # CHANGE longitudinal stabilty [ ]
    Cmde   = -1.001            # CHANGE elevator effectiveness [ ]

    # air density [kg/m^3]
    rho    = rho0 * pow( ((1+(lam * hp / Temp0))), (-((g / (lam*R)) + 1)))
    W      = m * g            # [N]       (aircraft weight)

    # Aircraft inertia (depend on t0):
    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
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

    CL = 2 * W / (rho * V ** 2 * S)                 # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha) ** 2 / (np.pi * A * e) # Drag coefficient [ ]

    # Stabiblity derivatives
    CX0    = W * np.sin(theta) / (0.5 * rho * V ** 2 * S)
    #CXu    = -0.095 #corrected original
    CXa    = +0.47966		# Positive! (has been erroneously negative since 1993)
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728
    print("CX0 = ", CX0)
    CZ0    = -W * np.cos(theta) / (0.5 * rho * V ** 2 * S)
    #CZu    = -0.37616 #original
    #CZa    = -5.74340
    #CZadot = -0.00350
    #CZq    = -5.66290
    CZde   = -0.69612

    Cmu    = +0.06990   #positive!
    #Cmadot = +0.17800   #positive!
    #Cmq    = -8.79415

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
    c1[0,0] = -2*muc*(c/V)
    c1[1,1] = (CZadot - 2*muc)*(c/V)
    c1[2,2] = -(c/V)
    c1[3,1] = Cmadot*(c/V)
    c1[3,3] = -2*muc*KY2*((c/V)**2)

    #c2 matrix
    c2 = np.zeros(s1)
    c2[0,0] = -CXu
    c2[0,1] = -CXa
    c2[0,2] = -CZ0
    c2[0,3] = -CXq*(c/V)
    c2[1,0] = -CZu
    c2[1,1] = -CZa
    c2[1,2] = -CX0
    c2[1,3] = -(CZq + 2*muc)*(c/V)
    c2[2,3] = -(c/V)
    c2[3,0] = -Cmu
    c2[3,1] = -Cma
    c2[3,3] = -Cmq*(c/V)

    #c3 matrix
    c3 = np.zeros(s2)
    c3[0,0] = -CXde
    c3[1,0] = -CZde
    c3[3,0] = -Cmde


    #Creating the different c-matrices (c4, c5 &c6) for asymmetrical flight

    #c4 matrix
    c4 = np.zeros(s1)
    c4[0,0] = (CYbdot - 2*mub)*(b/V)
    c4[1,1] = (-0.5)*(b/V)
    c4[2,2] = -4*mub*KX2*(b/V)*(b/(2*V))
    c4[2,3] = 4*mub*KXZ*(b/V)*(b/(2*V))
    c4[3,0] = Cnb*(b/V)
    c4[3,2] = 4*mub*KXZ*(b/V)*(b/(2*V))
    c4[3,3] = -4*mub*KZ2*(b/V)*(b/(2*V))

    #c5 matrix
    c5 = np.zeros(s1)
    c5[0,0] = CYb
    c5[0,1] = CL
    c5[0,2] = CYp*(b/(2*V))
    c5[0,3] = (CYr - 4*mub)*(b/(2*V))
    c5[1,2] = (b/(2*V))
    c5[2,0] = Clb
    c5[2,2] = Clp*(b/(2*V))
    c5[2,3] = Clr*(b/(2*V))
    c5[3,0] = Cnb
    c5[3,2] = Cnp*(b/(2*V))
    c5[3,3] = Cnr*(b/(2*V))

    #c6 matrix
    c6 = np.zeros(s3)
    c6[0,0] = -CYda
    c6[0,1] = -CYdr
    c6[2,0] = -Clda
    c6[2,1] = -Cldr
    c6[3,0] = -Cnda
    c6[3,1] = -Cndr

    # Time responses for unit steps:
    # t = np.linspace(t0,t0+ deltat, nsteps) -t0
    u = input_u

    # print("u and t:",u,t,sep='\n')
    # print("u shape:",u.shape)
    # print("t shape:",t.shape)
    if t.shape!=u.shape:
        print("Wrong slicing for input and time!\n")
        return -1

    #Now, we distinct between inputs:

    if input_type=="elevator":
        print("Calculating for elevator input...")
        #Symmetric system is triggered:

        #Creating the state matrix(A) and the input matrix(B) for symmetrical flight - xdot = c1^-1*c2*x c1^-1*c3*u = Ax + Bu
        A_s = np.dot(np.linalg.inv(c1), c2)
        B_s = np.dot(np.linalg.inv(c1), c3)
        C_s = np.identity(4)
        D_s = np.zeros((4, 1))

        #System in state-space
        sys_s = cm.StateSpace(A_s, B_s, C_s, D_s)
        poles_s = cm.pole(sys_s)
        # print("Eigenvalues of the symmetric system: ", poles_s,sep='\n') #verified

        # Time responses for unit steps:
        yout,tout,uout = cm.lsim(sys_s,u,t)   #general time response

        u_out_s =     yout[:,0]
        alpha_out_s = yout[:,1] + alpha*180/math.pi
        theta_out_s = yout[:,2] + theta*180/math.pi + 2
        q_out_s =     yout[:,3]

        if CXu != -0.095 or CZu != -0.37616 or CZa != -5.74340 or CZadot != -0.00350 or CZq != -5.66290 or Cmadot != 0.17800 or Cmq!= -8.79415:
            #Plotting....
            plotting(t,u_out_s,str("u Response for " +input_type+ " input, t0= "+ str(t0)),"u","m/s", "Modified Model")
            plotting(t,alpha_out_s,str("Alpha Response for " +input_type+ " input, t0= "+ str(t0)),r"$\alpha$","deg","Modified Model")
            plotting(t,theta_out_s,str("Theta Response for " +input_type+ " input, t0= "+ str(t0)),r"$\theta$","deg", "Modified Model")
            plotting(t,q_out_s,str("q Response for " +input_type+ " input, t0= "+ str(t0)),"$q$",r"deg/s", "Modified Model")
            print("\tPlotted!")

        else:
            plotting(t, u_out_s, str("u Response for " + input_type + " input, t0= " + str(t0)), "u", "m/s")
            plotting(t, alpha_out_s, str("Alpha Response for " + input_type + " input, t0= " + str(t0)), r"$\alpha$",
                     "deg")
            plotting(t, theta_out_s, str("Theta Response for " + input_type + " input, t0= " + str(t0)), r"$\theta$",
                     "deg")
            plotting(t, q_out_s, str("q Response for " + input_type + " input, t0= " + str(t0)), "$q$", r"deg/s")
            print("\tPlotted!")
        return q_out_s, theta_out_s, poles_s

    else:
        #Creating the state matrix(A) and the input matrix(B) for asymmetrical flight - y = c4^-1*c5*x c4^-1*c5*u = Ax + Bu
        A_a = -np.dot(np.linalg.inv(c4), c5)
        B_a = np.dot(np.linalg.inv(c4), c6)
        C_a = np.identity(4)
        #D_a depends on the input

        if input_type =="rudder":
            print("Calculating for rudder input...")
            D_a = np.zeros((4, 2))
            D_a[:,0] = 0   #we should check this...
            uarray = np.ones((len(t),2)) #step input
            uarray[:,1] = -u        #ADDED MINUS!!!!!
            uarray[:,0] = 0

        elif input_type=="aileron":
            print("Calculating for aileron input...")
            D_a = np.zeros((4, 2))
            D_a[:,1] = 1
            uarray = np.ones((len(t),2)) #step input
            uarray[:,0] = -u        #ADDED MINUS!!!!!
            uarray[:,1] = 0

        #System in state-space
        sys_a = cm.StateSpace(A_a, B_a, C_a, D_a)
        poles_a = cm.pole(sys_a)
        # print("Eigenvalues of the asymmetric system: ", poles_a) #verified


        yout,tout,uout = cm.lsim(sys_a,uarray,t)   #general time response for the input uarray

        beta_out_a = yout[:,0]
        phi_out_a = yout[:,1]
        p_out_a = yout[:,2]
        r_out_a = yout[:,3]

        if CXu!=-0.095 or CZu!=-0.37616 or CZa!= -5.74340 or CZadot!= -0.00350 or CZq!= -5.66290 or Cmadot!= 0.17800 or Cmq!= -8.79415:
        # #Plotting...
            plotting(t,beta_out_a,str("Beta Response for " + input_type +" input, t0= "+ str(t0)), r"$\beta$","deg", "Modified Model")
            plotting(t,phi_out_a,str("Phi Response for " +input_type + " input, t0= "+ str(t0)), r"$\phi$","deg", "Modified Model")
            plotting(t,p_out_a,str("p Response for " +input_type + " input, t0= "+ str(t0)) , r"$p$" ,"deg/s", "Modified Model")
            plotting(t,r_out_a,str("r Response for " +input_type + " input, t0= "+ str(t0)),  "$r$" ,r"deg/s", "Modified Model")
            print("\tPlotted!")
        else:
            plotting(t, beta_out_a, str("Beta Response for " + input_type + " input, t0= " + str(t0)), r"$\beta$","deg")
            plotting(t, phi_out_a, str("Phi Response for " + input_type + " input, t0= " + str(t0)), r"$\phi$", "deg")
            plotting(t, p_out_a, str("p Response for " + input_type + " input, t0= " + str(t0)), r"$p$", "deg/s")
            plotting(t, r_out_a, str("r Response for " + input_type + " input, t0= " + str(t0)), "$r$", r"deg/s")
            print("\tPlotted!")


        return poles_a

    return 1

#++++++++++++++++++++++++++++++++++++++ RMSE calculation +++++++++++++++++++++++++++++++++++++++++++++++++++

def RMSE(flight_data,simulation_data):
    return np.sqrt(((flight_data - simulation_data) ** 2).mean())

#++++++++++++++++++++++++++++++++++++++ Input & Output +++++++++++++++++++++++++++++++++++++++++++++++++++

# Simulation parameters for dynamic measurements:
# input: ph -> elevator def
#       short period -> elevator def
#       dutch roll -> rudder def
#       dutch roll_yd -> rudder def
#       ar -> aileron def
#       spi -> rudder def (pulse-like input)

# output: ph -> pitch, pitch rate
#       shp -> pitch, pitch rate
#       dr -> yaw rate, roll rate
#       dr_yd -> yaw rate, roll rate
#       ar -> roll, roll_rate
#       spi -> roll, yaw_rate



if __name__=="__main__":

    print("Collecting data...")

    t0_lst         = [53.5*60,58.6*60+3,60.1*60+5,60.95*60+5,57.0*60,3746]           #s
    deltat_lst     = [148, 5, 28 ,19 ,60 ,50]                                 #s -- these should match data_generator.py values (at the end)
    input_type_lst = ["elevator","elevator","rudder","rudder","aileron","aileron"]


  ################################# PHUGOID ###############################################
    # print("Phugoid")
    # CXu_modified= -0.15 #-0.165 # best one so far -0.14500000000000002
    # CZu_modified= -0.44616 # best one so far: -0.42616
    #
    #t0, deltat, utime_ph, u_ph, u_ph_p, u_ph_p_rate = phugoid()
    # fig_q = plotting(utime_ph ,u_ph_p_rate,str("q Response for " +input_type_lst[0]+ " input, t0= "+ str(t0)),"$q$",r"deg/s",label_name="Flight Test",title="Phugoid -> Pitch Rate")
    # fig_theta = plotting(utime_ph ,u_ph_p,str("Theta Response for " +input_type_lst[0]+ " input, t0= "+ str(t0)),r"$\theta$",r"deg",label_name="Flight Test",title="Phugoid -> Pitch")
    # q_out_s_p_original, theta_outs_s_p_original, eig_p_original =
    # q_out_s_p, theta_outs_s_p, eig_p = main(t0, deltat, utime_ph, input_type_lst[0], u_ph, CXu_modified, CZu_modified)
    # RMSE_original_phugoid_p_rate = RMSE(u_ph_p_rate,q_out_s_p)
    # RMSE_original_phugoid_p = RMSE(u_ph_p,theta_outs_s_p)
    # print(RMSE_original_phugoid_p, RMSE_original_phugoid_p_rate)
    # #...debugging....working?!
    # utime = utime_ph-utime_ph[0]                                            # translate the interval for better fitting
    # coeffs,cov = sp.curve_fit(double,utime,u_ph_p, p0=[5,0.01,-0.12,0.1,-1,1.7])  #initial guess is IMPORTANT
    # eig_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
    # print("Eigenvalues phugoid model: ", eig_p[3])
    # print("Eigenvalues phugoid period test: %r + j %r" %(coeffs[1],coeffs[2]))
    #
    # # print(utime)
    # # print(coeffs,cov,sep="\n")
    # # plt.figure("Testing")
    # # plt.plot(utime,double(utime,*coeffs),'r')
    # print("Phugoid relative error [%]: ", abs((abs(eig_p[3])-eig_test))*100/abs(eig_p[3]))   #first two are short period (large omega), last two are phugoid
    # print("\nTest:T1/2 = %r [s] and P= %r [s]" %(round(calcT(coeffs[1]),3),round(calcP(coeffs[2]),3)))

  ######################################## SHORT PERIOD ############################################

    print("Shord period")
    t0, deltat, utime_shp, u_shp, u_shp_p, u_shp_p_rate = short_period()
    # CZa_final_m = -5.74340
    # CZadot_final_m = -0.00350
    # CZq_final_m = -5.66290
    # Cmadot_final_m = 0.17800
    # Cmq_final_m = -8.79415

    CZa_final_m = -5.7434 #yes
    CZadot_final_m = 0.0035 # definitely no
    CZq_final_m = -5.6629 #no
    Cmadot_final_m = 0.178 #probably not
    Cmq_final_m = -0.01#8.79415 #yes



    q_out_s_shp_original, theta_outs_s_shp_original, eig_shp_original = main(t0, deltat, utime_shp, input_type_lst[0], u_shp)
    q_out_s_shp, theta_outs_s_shp, eig_shp = main(t0, deltat, utime_shp, input_type_lst[0], u_shp,CZa_final_m, CZadot_final_m, CZq_final_m, Cmadot_final_m, Cmq_final_m)
    RMSE_original_short_period_p_rate = RMSE(u_shp_p_rate,q_out_s_shp_original)
    RMSE_original_short_period_p = RMSE(u_shp_p, theta_outs_s_shp_original)
    print(RMSE_original_short_period_p, RMSE_original_short_period_p_rate)
    plotting(utime_shp,u_shp_p_rate,str("q Response for " +input_type_lst[1]+ " input, t0= "+ str(t0)),"$q$",r"deg/s",label_name="Flight Test", title="Short Period -> Pitch Rate")
    plotting(utime_shp,u_shp_p,str("Theta Response for " +input_type_lst[1]+ " input, t0= "+ str(t0)),r"$\theta$",r"deg",label_name="Flight Test", title="Short Period -> Pitch")
    #eig_s = main(t0,deltat,utime_shp,input_type_lst[1],u_shp)
    # print(eig_s)

    # #...debugging....
    utime = utime_shp-utime_shp[0]                                            # translate the interval for better fitting
    coeffs,cov = sp.curve_fit(doubleSP,utime,u_shp_p_rate)  #initial guess is IMPORTANT
    eig_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
    print("Eigenvalues short period model: ", eig_shp[1])
    print("Eigenvalues  short period test: %r + j %r" %(coeffs[1],coeffs[2]))

    # print(coeffs,cov,sep="\n")
    # plt.figure("Testing")
    # plt.plot(utime,doubleSP(utime,*coeffs),'r')
    # plt.plot(utime,u_shp_p_rate,'b')

    print("Eigenvalue error [%]: ", abs((abs(eig_shp[1])-eig_test))*100/abs(eig_shp[1]))   #first two are short period (large omega), last two are phugoid
    print("\nTest:T1/2 = %r [s] and P= %r [s]" %(round(calcT(coeffs[1]),3),round(calcP(coeffs[2]),3)))

  ######################################## DUTCH ROLL ##############################################

    # print("Dutch roll")
    # t0, deltat, utime_dr, u_dr, u_dr_y, u_dr_r = dutch_roll()
    #
    # plotting(utime_dr,u_dr_y,str("r Response for " +input_type_lst[2]+ " input, t0= "+ str(t0)),"$r$",r"1/s",label_name="Flight Test", title="Dutch Roll -> Yaw Rate")
    # plotting(utime_dr,u_dr_r,str("p Response for " +input_type_lst[2]+ " input, t0= "+ str(t0)),"$p$",r"1/s",label_name="Flight Test", title="Dutch Roll -> Roll Rate")
    # eig_dr = main(t0,deltat,utime_dr,input_type_lst[2],u_dr)
    # print("Eigenvalues dutch roll: ",eig_dr)
    #
    #
    # #...debugging....working?!
    # utime = utime_dr-utime_dr[0]                                            # translate the interval for better fitting
    # coeffs,cov = sp.curve_fit(simple,utime,u_dr_y, p0=[6,-0.6,2,np.pi/4])  #initial guess is IMPORTANT
    # eig_dr_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
    # print("Eigenvalues dutch roll test: %r + j %r" %(coeffs[1],coeffs[2]))
    #
    # # print(utime)
    # # print(coeffs,cov,sep="\n")
    # # plt.figure("Testing")
    # # plt.plot(utime,simple(utime,*coeffs),'r')
    # # plt.plot(utime,u_dr_y,'b')
    # print("Dutch roll relative error [%]: ", abs((abs(eig_dr[1])-eig_dr_test))*100/abs(eig_dr[1]))
    # print("T1/2 = %r [s] and P= %r [s]" %(round(calcT(coeffs[1]),3),round(calcP(coeffs[2]),3)))

    ######################################## DUTCH ROLL YD ###########################################

    # print("Dutch roll YD")
    # t0, deltat, utime_dr_yd, u_dr_yd, u_dr_yd_y, u_dr_yd_r= dutch_roll_yd()
    # plotting(utime_dr_yd,u_dr_yd_y,str("r Response for " +input_type_lst[3]+ " input, t0= "+ str(t0)),"$r$",r"1/s",label_name="Flight Test", title="Dutch Roll YD -> Yaw Rate")
    # plotting(utime_dr_yd,u_dr_yd_r,str("p Response for " +input_type_lst[3]+ " input, t0= "+ str(t0)),"$p$",r"1/s",label_name="Flight Test", title="Dutch Roll YD -> Roll Rate")
    # eig_dr_yd = main(t0,deltat,utime_dr_yd,input_type_lst[3],u_dr_yd)
    #
    # print("Numerical model YD eigenvalue:",eig_dr_yd[1])
    # #...debugging....working?!
    # utime = utime_dr_yd - utime_dr_yd[0]                                            # translate the interval for better fitting
    # coeffs,cov = sp.curve_fit(simple,utime,u_dr_yd_y, p0=[6,-0.6,2,np.pi/4])  #initial guess is IMPORTANT
    # eig_dr_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
    #
    # print("Eigenvalues dutch roll test: %r + j %r" %(coeffs[1],coeffs[2]))
    #
    # # print(utime)
    # # print(coeffs,cov,sep="\n")
    # # plt.figure("Testing")
    # # plt.plot(utime,simple(utime,*coeffs),'r')
    # # plt.plot(utime,u_dr_yd_y,'b')
    # print("Dutch roll relative error [%]: ", abs((abs(eig_dr_yd[1])-eig_dr_test))*100/abs(eig_dr_yd[1]))
    # print("T1/2 = %r [s] and P= %r [s]" %(round(calcT(coeffs[1]),3),round(calcP(coeffs[2]),3)))


   ######################################## APERIODIC ROLL ##########################################

    # print("Aperiodic roll")
    # t0, deltat, utime_ar, u_ar, u_ar_r, u_ar_r_rate = aperiodic_roll()
    # plotting(utime_ar,u_ar_r,str("Roll Response for " +input_type_lst[4]+ " input, t0= "+ str(t0)),"$\phi$",r"-",label_name="Flight Test",title="Aperiodic Roll -> Roll")
    # plotting(utime_ar,u_ar_r_rate,str("p Response for " +input_type_lst[4]+ " input, t0= "+ str(t0)),"$p$",r"1/s",label_name="Flight Test",title="Aperiodic Roll -> Roll Rate")
    # main(t0,deltat,utime_ar,input_type_lst[4],u_ar)
    # coeffs, cov = sp.curve_fit(simple, utime, u_dr_yd_y, p0=[6, -0.6, 2, np.pi / 4])  # initial guess is IMPORTANT
    # eig_dr_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
    # #
    #  print("Eigenvalues dutch roll test: %r + j %r" %(coeffs[1],coeffs[2]))
   # ####################################### SPIRAL ###############################################
   #  print("Spiral stability")
   #  t0, deltat, utime_spi, u_spi, u_spi_r, u_spi_y = spiral()
   #  plotting(utime_spi,u_spi_r,str("Phi Response for " +input_type_lst[5] + " input, t0= "+ str(t0)),"$\phi$",r"-",label_name="Flight Test",title="Spiral -> Roll")
   #  plotting(utime_spi,u_spi_y,str("r Response for " +input_type_lst[5]+ " input, t0= "+ str(t0)),"$r$",r"1/s",label_name="Flight Test",title="Spiral -> Yaw Rate")
   #  main(t0,deltat,utime_spi,input_type_lst[5],u_spi)

# sorry for using the same variable names...
