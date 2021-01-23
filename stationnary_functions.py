import numpy as np
import math as m
import os
import time as tt

p0=101325 #sea-level pressure
g0=9.81 #sea-level gravity acceleration
R=287.05  #air gas constant
T0=288.15  #sea-level temperature
delta=-0.0065  #lapse rate
gamma=1.4  #specific heat ratio
rho0=1.225  #sea-level density
c=2.0569  #Cessna Citation II chord
mu0=1.789E-05  #sea-level air dynamic viscosity
C=110.4 #Sutherland's constant
S=30 #Cessna Citation II wing area
b=15.911 #Cessna Citation II wing span
AR=b**2/S #Cessna Citation II wing aspect ratio
Ws=60500 #Cessna Citation II standard weight
d=0.686 #Cessna Citation II JT15D-4 engine diameter
Tcs= 0.2 #Standard thrust coefficient #CHANGE!!!!!
Tc= 0.2 #Thrust coefficient #CHANGE!!!!!
CmTc=-0.0064 #Dimensionless thrust moment arm  #CHECK if need CHANGE!!!!!

def speed_correction(VIAS,hp,TAT):
    #param VIAS: Indicated airspeed
    #param hp: Pressure altitude
    #param TAT: Total air temperature
    #return: Equivalent velocity, density and static temperature
    Vc=(VIAS-2)*0.514444
    p=p0*(1+delta*hp*0.3048/T0)**(-g0/delta/R)
    M=m.sqrt(2/(gamma-1)*((1+p0/p*((1+(gamma-1)*rho0*Vc**2/(2*gamma*p0))**(gamma/(gamma-1))-1))**((gamma-1)/gamma)-1))
    T=(TAT+273.15)/(1+(gamma-1)/2*M**2)
    a=m.sqrt(gamma*R*T)
    rho=p/(R*T)
    Vt=M*a
    Ve=Vt*m.sqrt(rho/rho0)

    return Ve, Vt, rho, T, M

#1 alittude-->many TAT
def Reynolds(T_avg,rho,Vt):
    #param T_avg: Average temperature
    #param rho: Pressure altitude
    #param Vt: Total air temperature
    #return: Reynolds number

    mu=mu0*(T0+C)/(T_avg+C)*(T_avg/T0)**1.5
    Re=rho*Vt*c/mu

    return Re

def aero_coefficient(W,T,rho,Vt,alpha):
    #param W: Aircraft weight
    #param T: Thrust
    #param rho: Density
    #param V: True airspeed
    #return: Lift coefficient, Drag coefficient

    CL=2*W/(rho*S*Vt**2)
    CD=2*T/(rho*S*Vt**2)


    return CL, CD


def elevator_effectiveness(CN,delta_xcg,delta_de):
    #param Ve: Equivalent airspeed
    #param W: Aircraft weight
    #return: Elevator effectiveness

    #Assume no change in CN
    Cmd=-CN*delta_xcg/delta_de/c

    return Cmd


def weight_correction(Ve,Fe,W):
    #param Ve: Equivalent airspeed
    #param W: Aircraft weight
    #return: Weight corrected equivalent speed and elevetor control force

    Ve_tilde=Ve*m.sqrt(Ws/W)
    Fe_star=Fe*Ws/W

    return Ve_tilde, Fe_star

#FUnction that obtains Cmd

def thrust_correction(de,Cmd,Tcs,Tc):
    #param deq: Elevator deflection
    #param deq: Elevator effectiveness
    #param W: Standard thrust coefficient
    #param W: Thrust coefficient
    #return: Thrust corrected elevetor deflection

    de_star=de-CmTc/Cmd*(Tcs-Tc)

    return de_star

def longitudial_stability(alpha,de,Cmd):
    #param alpha: Angle of attack
    #param de: Elevator deflection
    #param Cmd: Elevator effectiveness
    #return: Longitudial stability

    m,c=np.polyfit(alpha,de,1)
    Cma=-m*Cmd

    return Cma


def thrust(hp,VIAS,FFL,FFR,TAT):
    #param hp: Pressure altitude
    #param VIAS: Indicated airspeed
    #param FFL: Fuel flow left
    #param FFR: Fuel flow right
    #param TAT: Total temperature
    #return: Creates matlab.dat file and executes Thrust.exe

    hp=[i*0.3048 for i in hp]
    FFL=[i*0.453582/3600  for i in FFL]
    FFR=[i*0.453582/3600  for i in FFR]

    M=[]
    deltaT=[]
    for i in range(0,len(VIAS)):
        Vc=(VIAS[i]-2)*0.514444
        p=p0*(1+delta*hp[i]*0.3048/T0)**(-g0/delta/R)
        a=m.sqrt(2/(gamma-1)*((1+p0/p*((1+(gamma-1)*rho0*Vc**2/(2*gamma*p0))**(gamma/(gamma-1))-1))**((gamma-1)/gamma)-1))
        b=(TAT[i]+273.15)/(1+(gamma-1)/2*a**2)-(T0+delta*hp[i])
        M.append(a)
        deltaT.append(b)

    root = r"C:\Users\Mathieu Pellé\Documents\AE\3\Q3\AE3212-II\FD"
    os.chdir(root)
    file = str("matlab.dat")

    try:
        fid = open(file,"w")
    except:
        print("Fail cannot be found/opne!\n")

    for i in range(0,len(M)):
        fid.write(str(hp[i]) + " ")
        fid.write(str(M[i]) + " ")
        fid.write(str(deltaT[i]) + " ")
        fid.write(str(FFL[i])+ " ")
        fid.write(str(FFR[i])+ " ")

        #For next data point
        fid.write("\n")
    os.startfile("thrust(1).exe")

def thrust_coefficient(Th,rho,Vt):
    #param Th: Thrust
    #param rho: Densitu
    #param V: Velocity
    #return: Thrust coefficient

    cT=Th/(0.5*rho*Vt**2*d**2) #CHECK Formula

    return cT
