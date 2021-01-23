# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:25:48 2020

@author: Group B07
Worked on plotting
@version: 3.5 22/03/2020
"""

import numpy as np
import control.matlab as cm
import matplotlib.pyplot as plt
import numpy.linalg
from massbalance import *
from data_generator_complete import *
import warnings
import scipy.optimize as sp
import time
import sys
from progress.bar import Bar
#import sublime, sublime_plugin
#import time


warnings.filterwarnings("ignore")

plt.close('all')


# +++++++++++++++++++++++++++++++++ Helper Functions ++++++++++++++++++++++++++++++++++++++++++++++

# ---------------------------------- fitting ------------------------------


def double(t, A1, eps1, omega_d1, A2, eps2, omega_d2):
    return 1.2 + A1 * np.exp(eps1 * t) * np.sin(omega_d1 * t + 0) + A2 * np.exp(eps2 * t) * np.cos(omega_d2 * t + 0)


def doubleSP(t, A1, eps1, omega_d1, A2, eps2, omega_d2, phi01, phi02):
    return A1 * np.exp(eps1 * t) * np.sin(omega_d1 * t + phi01) + A2 * np.exp(eps2 * t) * np.cos(omega_d2 * t + phi02)


def simple(t, A1, eps1, omega_d1, phi01):
    return A1 * np.exp(eps1 * t) * np.sin(omega_d1 * t + phi01)


def short(t, A, eps, omega, B, C):
    return A * np.exp(t * eps) * np.cos(omega * t + b)


def calcT(eps):
    return -np.log(1 / 2) / abs(eps)


def calcP(omega):
    return 2 * np.pi / abs(omega)


# ----------------------------------- plots ---------------------------------
#def plotting(x, y, name, variable, unit, label_name="Simulation", title=None, mins=False):
    #"""Use this for plotting. It returns the figure so we can add more to it"""
    #ax = plt.figure(str(name))
    # ax.legend("best")

    ##if mins:
       # x /= 60  # change time to mins from secs
       # plt.xlabel("t [min]")
    #else:
        #plt.xlabel("t [s]")

    #if title != None:
        #plt.title(str(title))

    ##if label_name == "Simulation":
       # plt.plot(x - x[0], y, '--', label=label_name)
    #else:
       # plt.plot(x - x[0], y, label=label_name)

   # lab = str(str(variable) + " " + "[" + unit + "]")
   # plt.legend(loc='best')
   # plt.ylabel(lab)
   # plt.grid(True)
    # plt.savefig(title)
   # plt.show()

   # return ax


# +++++++++++++++++++++++++++++++++++ Global variables+++++++++++++++++++++++++++++++++++++++++++++++

# Citation 550 - Linear simulation
# Aircraft geometry

S = 30.00  # wing area [m^2]
Sh = 0.2 * S  # stabiliser area [m^2]
Sh_S = Sh / S  # [ ]
lh = 0.71 * 5.968  # tail length [m]
c = 2.0569  # mean aerodynamic cord [m]
lh_c = lh / c  # [ ]
b = 15.911  # wing span [m]
bh = 5.791  # stabilser span [m]
A = b ** 2 / S  # wing aspect ratio [ ]
Ah = bh ** 2 / Sh  # stabilser aspect ratio [ ]
Vh_V = 1  # [ ]
ih = -2 * np.pi / 180  # stabiliser angle of incidence [rad]

oew = 4157.174  # Operational Empty Weight [kg]
m_payload = 765  # Payload mass [kg]

# Aerodynamic properties
e = 0.5466  # Oswald factor [ ]
CD0 = 0.01980  # Zero lift drag coefficient [ ]
CLa = 4.547  # Slope of CL-alpha curve [ ]

# Constant values concerning atmosphere and gravity
rho0 = 1.2250  # air density at sea level [kg/m^3]
lam = -0.0065  # temperature gradient in ISA [K/m]
Temp0 = 288.15  # temperature at sea level in ISA [K]
R = 287.05  # specific gas constant [m^2/sec^2K]
g = 9.81  # [m/sec^2] (gravity constant)

# Loading in data

altlst = np.genfromtxt("Data_SI_correct/Dadc1_altSI.txt", skip_header=2)  # reading in the altitude values
taslst = np.genfromtxt("Data_SI_correct/Dadc1_tasSI.txt", skip_header=2)  # reading in the true airspeed values
aoalst = np.genfromtxt("Data_SI_correct/vane_AOASI.txt", skip_header=2)  # reading in the angle of attack
pitchlst = np.genfromtxt("Data_SI_correct/Ahrs1_bPitchRateSI.txt", skip_header=2)  # reading in the pitch
# Simulation parameters:
nsteps = 10 ** 3


# +++++++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++++++++++++++++++
def main(t0, deltat, t, input_type, input_u, CZu, CXu, CZa, Cmq):
    """Input type: elevator
                    rudder
                    airleron"""

    # Find time
    idx = np.where(timelst == t0)[0]

    # Flight condition
    #  m_fuel    = 1197.484        # CHANGE Total fuel mass [kg]

    hp = altlst[idx]  # CHANGE pressure altitude in the stationary flight condition [m]
    V = taslst[idx]  # CHANGE true airspeed in the stationary flight condition [m/sec]
    alpha = np.radians(aoalst[idx])  # angle of attack in the stationary flight condition [rad]
    theta = np.radians(pitchlst[idx])  # pitch angle in the stationary flight condition [rad]
    gamma = theta - alpha  # CHANGE flight path angle -

    # Aircraft mass
    m = mass(t0)  # mass [kg]

    # Longitudinal stability
    Cma = -0.4435  # CHANGE longitudinal stabilty [ ]
    Cmde = -1.001  # CHANGE elevator effectiveness [ ]

    # air density [kg/m^3]
    rho = rho0 * pow(((1 + (lam * hp / Temp0))), (-((g / (lam * R)) + 1)))
    W = m * g  # [N]       (aircraft weight)

    # Aircraft inertia (depend on t0):
    muc = m / (rho * S * c)
    mub = m / (rho * S * b)
    KX2 = 0.019
    KZ2 = 0.042
    KXZ = 0.002
    KY2 = 1.25 * 1.114

    # Aerodynamic constants:

    Cmac = 0  # Moment coefficient about the aerodynamic centre [ ]
    CNwa = CLa  # Wing normal force slope [1/rad]
    CNha = 2 * np.pi * Ah / (Ah + 2)  # Stabiliser normal force slope [ ]
    depsda = 4 / (A + 2)  # Downwash gradient [ ]

    # Lift and drag coefficient (depend on t0):

    CL = 2 * W / (rho * V ** 2 * S)  # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha) ** 2 / (np.pi * A * e)  # Drag coefficient [ ]

    # Stabiblity derivatives
    CX0 = W * np.sin(theta) / (0.5 * rho * V ** 2 * S)
    # CXu    = -0.095         #corrected
    CXa = +0.47966  # Positive! (has been erroneously negative since 1993)
    CXadot = +0.08330
    CXq = -0.28170
    CXde = -0.03728

    CZ0 = -W * np.cos(theta) / (0.5 * rho * V ** 2 * S)
    # CZu    = -0.37616
    CZa = -5.74340
    CZadot = -0.00350
    CZq = -5.66290
    CZde = -0.69612

    Cmu = +0.06990  # positive!
    Cmadot = +0.17800  # positive!
    Cmq = -8.79415

    #CYb = -0.7500
    #CYbdot = 0
    #CYp = -0.0304
    #CYr = +0.8495
    #CYda = -0.0400
    #CYdr = +0.2300

    #Clb = -0.10260
    #Clp = -0.71085
    #Clr = +0.23760
    #Clda = -0.23088
    #Cldr = +0.03440

    #Cnb = +0.1348
    #Cnbdot = 0
    #Cnp = -0.0602
    #Cnr = -0.2061
    #Cnda = -0.0120
    #Cndr = -0.0939

    # c-matrix dimensions
    s1 = (4, 4)
    s2 = (4, 1)
    s3 = (4, 2)

    # Creating the different c-matrices (c1, c2 &c3) for symmetrical flight
    # c1 matrix
    c1 = np.zeros(s1)
    c1[0, 0] = -2 * muc * (c / V)
    c1[1, 1] = (CZadot - 2 * muc) * (c / V)
    c1[2, 2] = -(c / V)
    c1[3, 1] = Cmadot * (c / V)
    c1[3, 3] = -2 * muc * KY2 * ((c / V) ** 2)

    # c2 matrix
    c2 = np.zeros(s1)
    c2[0, 0] = -CXu
    c2[0, 1] = -CXa
    c2[0, 2] = -CZ0
    c2[0, 3] = -CXq * (c / V)
    c2[1, 0] = -CZu
    c2[1, 1] = -CZa
    c2[1, 2] = -CX0
    c2[1, 3] = -(CZq + 2 * muc) * (c / V)
    c2[2, 3] = -(c / V)
    c2[3, 0] = -Cmu
    c2[3, 1] = -Cma
    c2[3, 3] = -Cmq * (c / V)

    # c3 matrix
    c3 = np.zeros(s2)
    c3[0, 0] = -CXde
    c3[1, 0] = -CZde
    c3[3, 0] = -Cmde

    # Creating the different c-matrices (c4, c5 &c6) for asymmetrical flight

    # Time responses for unit steps:
    # t = np.linspace(t0,t0+ deltat, nsteps) -t0
    u = input_u

    # print("u and t:",u,t,sep='\n')
    # print("u shape:",u.shape)
    # print("t shape:",t.shape)
    if t.shape != u.shape:
        print("Wrong slicing for input and time!\n")
        return -1

    # Now, we distinct between inputs:

    if input_type == "elevator":
        #print("Calculating for elevator input...")
        # Symmetric system is triggered:

        # Creating the state matrix(A) and the input matrix(B) for symmetrical flight - xdot = c1^-1*c2*x c1^-1*c3*u = Ax + Bu
        A_s = np.dot(np.linalg.inv(c1), c2)
        B_s = np.dot(np.linalg.inv(c1), c3)
        C_s = np.identity(4)
        D_s = np.zeros((4, 1))

        # System in state-space
        sys_s = cm.StateSpace(A_s, B_s, C_s, D_s)
        poles_s = cm.pole(sys_s)
        # print("Eigenvalues of the symmetric system: ", poles_s,sep='\n') #verified

        # Time responses for unit steps:
        yout, tout, uout = cm.lsim(sys_s, u, t)  # general time response

        u_out_s = yout[:, 0]
        alpha_out_s = yout[:, 1] + alpha
        theta_out_s = yout[:, 2] + theta
        q_out_s = yout[:, 3]

        # Plotting....
        # plotting(t,u_out_s,str("u Response for " +input_type+ " input, t0= "+ str(t0)),"u","m/s")
        # plotting(t,alpha_out_s,str("Alpha Response for " +input_type+ " input, t0= "+ str(t0)),r"$\alpha$","deg")
        # plotting(t,theta_out_s,str("Theta Response for " +input_type+ " input, t0= "+ str(t0)),r"$\theta$","deg")
        # plotting(t,q_out_s,str("q Response for " +input_type+ " input, t0= "+ str(t0)),"$q$",r"deg/s")
        # print("\tPlotted!")
        return theta_out_s, q_out_s, poles_s

    return 1


# ++++++++++++++++++++++++++++++++++++++ Input & Output +++++++++++++++++++++++++++++++++++++++++++++++++++

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


#if __name__ == "__main__":

# print("Collecting data...")

# t0_lst         = [53.5*60,58.6*60+3,60.1*60+5,60.95*60+5,57.0*60,3746]           #s
# deltat_lst     = [148, 5, 28 ,19 ,60 ,50]                                 #s -- these should match data_generator.py values (at the end)
input_type_lst = ["elevator","elevator","rudder","rudder","aileron","aileron"]


################################# PHUGOID ###############################################
# print("Phugoid")
# t0, deltat, utime_ph, u_ph, u_ph_p, u_ph_p_rate = phugoid()
# fig_q = plotting(utime_ph ,u_ph_p_rate,str("q Response for " +input_type_lst[0]+ " input, t0= "+ str(t0)),"$q$",r"deg/s",label_name="Flight Test",title="Phugoid -> Pitch Rate")
# fig_theta = plotting(utime_ph ,u_ph_p,str("Theta Response for " +input_type_lst[0]+ " input, t0= "+ str(t0)),r"$\theta$",r"deg",label_name="Flight Test",title="Phugoid -> Pitch")
# eig_p = main(t0,deltat,utime_ph,input_type_lst[0],u_ph)
#
# #...debugging....working?!
# utime = utime_ph-utime_ph[0]                                            # translate the interval for better fitting
# coeffs,cov = sp.curve_fit(double,utime,u_ph_p, p0=[5,0.01,-0.12,0.1,-1,1.7])  #initial guess is IMPORTANT
# eig_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
# print("Eigenvalues phugoid model: ", eig_p[3])
# print("Eigenvalues phugoid period test: %r + j %r" %(coeffs[1],coeffs[2]))
#
##print(utime)
##print(coeffs,cov,sep="\n")
# plt.figure("Testing")
# plt.plot(utime,double(utime,*coeffs),'r')
# print("Phugoid relative error [%]: ", abs((abs(eig_p[3])-eig_test))*100/abs(eig_p[3]))   #first two are short period (large omega), last two are phugoid
# print("\nTest:T1/2 = %r [s] and P= %r [s]" %(round(calcT(coeffs[1]),3),round(calcP(coeffs[2]),3)))

######################################## SHORT PERIOD ############################################

# print("Shord period")
# t0, deltat, utime_shp, u_shp, u_shp_p, u_shp_p_rate = short_period()
# plotting(utime_shp,u_shp_p_rate,str("q Response for " +input_type_lst[1]+ " input, t0= "+ str(t0)),"$q$",r"deg/s",label_name="Flight Test", title="Short Period -> Pitch Rate")
# plotting(utime_shp,u_shp_p,str("Theta Response for " +input_type_lst[1]+ " input, t0= "+ str(t0)),r"$\theta$",r"deg",label_name="Flight Test", title="Short Period -> Pitch")
# eig_s = main(t0,deltat,utime_shp,input_type_lst[1],u_shp)
# # print(eig_s)
#
# # #...debugging....
# utime = utime_shp-utime_shp[0]                                            # translate the interval for better fitting
# coeffs,cov = sp.curve_fit(doubleSP,utime,u_shp_p_rate)  #initial guess is IMPORTANT
# eig_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)                      #absolute value
# print("Eigenvalues short period model: ", eig_s[1])
# print("Eigenvalues  short period test: %r + j %r" %(coeffs[1],coeffs[2]))
#
# # print(coeffs,cov,sep="\n")
# # plt.figure("Testing")
# # plt.plot(utime,doubleSP(utime,*coeffs),'r')
# # plt.plot(utime,u_shp_p_rate,'b')
#
# print("Eigenvalue error [%]: ", abs((abs(eig_s[1])-eig_test))*100/abs(eig_s[1]))   #first two are short period (large omega), last two are phugoid
# print("\nTest:T1/2 = %r [s] and P= %r [s]" %(round(calcT(coeffs[1]),3),round(calcP(coeffs[2]),3)))
#++++++++++++++++++++++++++++++++++++++ RMSE calculation +++++++++++++++++++++++++++++++++++++++++++++++++++

def RMSE(flight_data,simulation_data):
    return np.sqrt(((flight_data - simulation_data) ** 2).mean())

  #++++++++++++++++++++++++++++++++++++++ Phugoid Model Improvement +++++++++++++++++++++++++++++++++++++++++++++++++++
#RMSE_original_phugoid_p_rate = 0.18537027365097794
#RMSE_original_phugoid_p = 3.248328645062731

#t0, deltat, utime_ph, u_ph, u_ph_p, u_ph_p_rate = phugoid()
#utime = utime_ph-utime_ph[0]                                            # translate the interval for better fitting
#coeffs,cov = sp.curve_fit(double,utime,u_ph_p, p0=[5,0.01,-0.12,0.1,-1,1.7])  #initial guess is IMPORTANT
#eig_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)
#rel_error_final = 100
#CZu_final = -0.37616
#CXu_final = -0.095
#interval = 0.01
#step = 30
#lst_CXu =np.linspace(CXu_final - interval, CXu_final + interval, step)
#lst_CZu =np.linspace(CZu_final - interval, CZu_final + interval, step)
#from progress.bar import Bar
#bar = Bar('Processing', max=step**2)
#error_real_final = 100
#error_imag_final = 100
#min_RMSE_phugoid_p_rate = RMSE_original_phugoid_p_rate
#min_RMSE_phugoid_p = RMSE_original_phugoid_p
#for i in lst_CZu:
#  for j in lst_CXu:
#    bar.next()
#    theta_out_s_p, q_out_s_p, eig_p = main(t0,deltat,utime_ph,input_type_lst[0],u_ph, i, j, -5.74340, -0.00350, -5.66290, 0.17800, -8.79415) #i = CZu, j = CXu CZa, CZadot, CZq, Cmadot, Cmq
    # #...debugging....working?!
#    real_num = eig_p[3].real
#    imag_num = eig_p[3].imag
#    real_ft  = eig_test.real
#    imag_ft = eig_test.imag
#    error_real = abs(real_num - real_ft)
#    error_imag = abs(imag_num - real_ft)
#    rel_error = abs((abs(eig_p[3])-eig_test))*100/abs(eig_p[3])
#    if RMSE(q_out_s_p,u_ph_p_rate) < min_RMSE_phugoid_p_rate and  RMSE(theta_out_s_p,u_ph_p) < min_RMSE_phugoid_p:
#      if error_real < error_real_final and error_imag < error_imag_final:
#          print('x')
#          error_real_final = error_real
#          error_imag_final = error_imag
#          rel_error_final = rel_error
#          CZu_final = i
#          CXu_final = j
#   min_RMSE_phugoid_p_rate=RMSE(q_out_s_p,u_ph_p_rate)
#    min_RMSE_phugoid_p=RMSE(theta_out_s_p,u_ph_p)
#bar.finish()

#print("The relative error is:", rel_error_final)
#print("CZu = ", CZu_final)
#print("CXu = ", CXu_final)

#CZu =  -0.44616 -> period
#CXu =  -0.165 -> amplitude


#CXu_modified= -0.15 # best one so far -0.14500000000000002
#CZu_modified= -0.44616 # best one so far: -0.42616
########################################### Short Period Motion Model Improvement #############################
RMSE_original_short_period_p_rate = 100#0.39683309816093915
RMSE_original_short_period_p = 100#3.087676046949498
min_RMSE_short_period_p_rate = RMSE_original_short_period_p_rate
min_RMSE_short_period_p = RMSE_original_short_period_p
t0, deltat, utime_shp, u_shp, u_shp_p, u_shp_p_rate = short_period()
utime = utime_shp-utime_shp[0]                                            # translate the interval for better fitting
coeffs,cov = sp.curve_fit(doubleSP,utime,u_shp_p_rate)  #initial guess is IMPORTANT
eig_test = np.sqrt(coeffs[1]**2 + coeffs[2]**2)

rel_error_final = 100
CZa_final = -5.74340
#CZadot_final = -0.00350
#CZq_final = -5.66290
#Cmadot_final = 0.17800
Cmq_final = -8.79415
interval = 30
step = 30
lst_CZa =np.linspace(CZa_final - interval, CZa_final + interval, step)
#lst_CZadot =np.linspace(CZadot_final - interval, CZadot_final + interval, step)
#lst_CZq =np.linspace(CZq_final - interval, CZq_final + interval, step)
#lst_Cmadot =np.linspace(Cmadot_final - interval, Cmadot_final + interval, step)
lst_Cmq =np.linspace(Cmq_final - interval, Cmq_final + interval, step)
from progress.bar import Bar
bar = Bar('Processing', max=100)
error_real_final = 100
error_imag_final = 100
for i in lst_CZa:
#  for j in lst_Czadot:
    #for k in lst_CZq:
    #for l in lst_Cmadot:
  for j in lst_Cmq:
      bar.next()
      theta_out_s_shp, q_out_s_shp, eig_s = main(t0,deltat,utime_shp,input_type_lst[1],u_shp, -0.37616,-0.095 , i, j) #CZu, CXu, i=CZa, j=Cmq
     # #...debugging....working?!
      real_num = eig_s[1].real
      imag_num = eig_s[1].imag
      real_ft  = eig_test.real
      imag_ft = eig_test.imag
      error_real = abs(real_num - real_ft)
      error_imag = abs(imag_num - real_ft)
      rel_error = abs((abs(eig_s[1])-eig_test))*100/abs(eig_s[1])
      if RMSE(q_out_s_shp,u_shp_p_rate) < min_RMSE_short_period_p_rate and  RMSE(theta_out_s_shp,u_shp_p) < min_RMSE_short_period_p:
         if error_real < error_real_final and error_imag < error_imag_final:
            print('x')
            error_real_final = error_real
            error_imag_final = error_imag
            rel_error_final = rel_error
            CZa_final = i
            #CZadot_final = j
            #CZq_final = k
            #Cmadot_final = l
            Cmq_final = j
      min_RMSE_short_period_p_rate=RMSE(q_out_s_shp,u_shp_p_rate)
      min_RMSE_short_period_p=RMSE(theta_out_s_shp,u_shp_p)
bar.finish()

print("The relative error is:", rel_error_final)
print("CZa = ", CZa_final)
#print("CZadot = ", CZadot_final)
#print("CZq = ", CZq_final)
#print("Cmadot = ", Cmadot_final)
print("Cmq = ", Cmq_final)


"""Ok so these are the initial RMSE for the phugoid:
RMSE_original_phugoid_p_rate = 0.18537027365097794
RMSE_original_phugoid_p = 3.248328645062731"""
