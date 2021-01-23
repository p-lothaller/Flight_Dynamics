# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:25:48 2020

@author: vladg
@version: main - 3.4 End of 21/03/2020
"""

import numpy as np
import pandas as pd
import os
import scipy.io as sp
import matplotlib.pyplot as plt


#Change this for your computer!
root = r"C:\Users\lotha\Flight_Dyanmics_SVV"
mydir = r"C:\Users\lotha\Flight_Dyanmics_SVV\Data_SI_correct"


def calcT(h): #ISA tmeperature
    return 15.0 - 0.0065*h

def getData(mydir):
    os.chdir(mydir)

    lst = []
    for root, dirs, files in os.walk(mydir):
        for file in files:
            # print(file)
            fullData = pd.read_csv(file,"\n")
            data = fullData[2:]
            name = fullData.columns[0]
            # units = str(fullData[0])
            local = {"name":name,"data":data}
            lst.append(local)

    return lst

def genThrust(fid,h,M,deltaT,FMFr,FMFl):
    """Generates the tthrust.txt file for the thrust calculations
 Input: altitude                      h [m]
         mach number                  M [-]
         delta real-ISA               T [C]
         Fuel Mass Flow right engine  [kg/s]
         Fuel Mass flow left engine   [kg/s]
         All from Reader/pg32
 Output: matlab.dat file for thrust.exe app"""



    # Open the file

    #Create the command file
    fid.write( str(h) + " ")
    fid.write(str(M) + " ")
    fid.write(str(deltaT) + " ")
    fid.write(str(FMFr)+ " ")
    fid.write(str(FMFl)+ " ")

    #For next data point
    fid.write("\n")


    return 1

def plottingData(x,y,name,title,variable,unit,mins=False):
    """Use this for plottingData."""
    ax = plt.figure(str(name))
    # ax.legend("best")

    if mins:
        x/=60  #change time to mins from secs
        plt.xlabel("t [min]")
    else:
        plt.xlabel("t [s]")

    # if title!= None:
    #     plt.title(str(title))

    plt.plot(x,y,label=name)


    lab = str(str(variable)+" "+"["+unit+"]")
    plt.ylabel(lab)
    plt.grid(True)
    # plt.savefig(title)
    plt.show()



lst = getData(mydir)

# for item in lst:
#     if "Angle of attack" in item["name"]:
#         alphatab = np.array(item["data"])
#         print("Alpha:",alphatab)
#         lst.remove(item)

    # if "Time" in item["name"]:
    #     timetab = np.array(item["data"])
    #     print("Time:",timetab)




#Data aquisition & Unit conversion (when needed)
timetab = np.genfromtxt("timeSI.txt",dtype=float,skip_header=2,delimiter='\n')  #s DONT CHANGE!


Htab = np.genfromtxt("Dadc1_altSI.txt",dtype=float,skip_header=2,delimiter='\n')  #m
Ttab = np.genfromtxt("Dadc1_tatSI.txt",dtype=float,skip_header=2,delimiter='\n')    #C
Mtab = np.genfromtxt("Dadc1_machSI.txt",dtype=float,skip_header=2,delimiter='\n')  #-
Tisa = calcT(Htab)   #C

FMFr_tab = np.genfromtxt("rh_engine_FMFSI.txt",dtype=float,skip_header=2,delimiter='\n')  #kg/s
FMFl_tab = np.genfromtxt("lh_engine_FMFSI.txt",dtype=float,skip_header=2,delimiter='\n') #kg/s

pitch_tab = np.genfromtxt("Ahrs1_PitchSI.txt",dtype=float,skip_header=2,delimiter='\n')
pitch_rate_tab =np.genfromtxt("Ahrs1_bPitchRateSI.txt",dtype=float,skip_header=2,delimiter='\n')

roll_tab =  np.genfromtxt("Ahrs1_RollSI.txt",dtype=float,skip_header=2,delimiter='\n')
roll_rate_tab =  np.genfromtxt("Ahrs1_bRollRateSI.txt",dtype=float,skip_header=2,delimiter='\n')

yaw_rate_tab = np.genfromtxt("Ahrs1_bYawRateSI.txt",dtype=float,skip_header=2,delimiter='\n')

alphatab = np.genfromtxt("vane_AOASI.txt",dtype=float,skip_header=2,delimiter='\n')  #deg
gtab = np.genfromtxt("Ahrs1_VertAccSI.txt",dtype=float,skip_header=2,delimiter='\n')
elevtab =  np.genfromtxt("delta_eSI.txt",dtype=float,skip_header=2,delimiter='\n')
ruddertab = np.genfromtxt("delta_rSI.txt",dtype=float,skip_header=2,delimiter='\n')
ailerontab = np.genfromtxt("delta_aSI.txt",dtype=float,skip_header=2,delimiter='\n')

elevtrimtab = np.genfromtxt("elevator_dteSI.txt",dtype=float,skip_header=2,delimiter='\n')
phitab = np.genfromtxt("Ahrs1_RollSI.txt",dtype=float,skip_header=2,delimiter='\n')

plt.plot(timetab,ailerontab,'r-')
# plt.plot(timetab,phitab,'g-')
# plt.plot(timetab,ruddertab,'b-')
plt.grid(True)
os.chdir(root)
file = str("matlab.dat")

try:
    fid = open(file,"w")
except:
    print("Fail cannot be found/opne!\n")

for i in range(10):
    genThrust(fid,Htab[i],Mtab[i],abs(Ttab[i]-Tisa[i]),FMFl_tab[i],FMFr_tab[i])

fid.close()

#!!!!!!!!!!!!!!!!!!!!!  ADD BELOW !!!!!!!!!!!

def getInput(tab,timetab,t0,deltat): #to be exported; will be used for output as well
    """Returns the sliced array from tab when
    the time values (in seconds!) are contained in the interval (t0,t0+deltat)
    Output: tab[slices],time[slices]"""

    return tab[np.where((t0+deltat>timetab) & (t0<timetab))], timetab[np.where((t0+deltat>timetab) & (t0<timetab))]


######################################## PHUGOID ###############################################

def phugoid(t0=53.5*60, deltat=148, plot_input=False, plot_output=False):

    #input -> elevator deflection
    u_ph, utime_ph = getInput(elevtab,timetab,t0,deltat)
    # print("Local u and t:",u_ph,utime_ph,sep="\n")

    #output -> pitch
    u_ph_p, utime_ph_p = getInput(pitch_tab,timetab,t0,deltat)

    #       -> pitch rate
    u_ph_p_rate, utime_ph_p_rate = getInput(pitch_rate_tab, timetab, t0, deltat)

    ####plottingData####
    if plot_input==True:
        plottingData(utime_ph, u_ph, name="elevator_def_Phugoid", title="Phugoid Input", variable="${\delta}_e$", unit="deg", mins=False)

    if plot_output==True:
        plottingData(utime_ph_p, u_ph_p, name="pitch_Phugoid", title="Phugoid -> Pitch", variable="${\Theta}$", unit="deg", mins=False)
        plottingData(utime_ph_p_rate, u_ph_p_rate, name="pitch_rate_Phugoid", title="Phugoid -> Pitch Rate", variable="$q$", unit="deg/s",mins=False)

    return  t0, deltat, utime_ph, u_ph, u_ph_p, u_ph_p_rate


######################################## SHORT PERIOD ############################################

def short_period(t0=58.6*60+3, deltat=5, plot_input=False, plot_output=False):
    # input -> elevator deflection
    u_shp, utime_shp = getInput(elevtab, timetab, t0, deltat)

    # output -> pitch
    u_shp_p, utime_shp_p = getInput(pitch_tab, timetab, t0, deltat)

    #       -> pitch rate
    u_shp_p_rate, utime_shp_p_rate = getInput(pitch_rate_tab, timetab, t0, deltat)

    ####plottingData####
    if plot_input == True:
        plottingData(utime_shp, u_shp, name="elevator_def_Short_Period", title="Short Period Input", variable="${\delta}_e$",
                 unit="deg", mins=False)

    if plot_output == True:
        plottingData(utime_shp_p, u_shp_p, name="pitch_Short_Period", title="Short Period -> Pitch", variable="${\Theta}$", unit="deg",
                 mins=False)
        plottingData(utime_shp_p_rate, u_shp_p_rate, name="pitch_rate_Short_Period", title="Short Period -> Pitch Rate",
                 variable="$q$", unit="deg/s", mins=False)

    return t0, deltat, utime_shp, u_shp, u_shp_p, u_shp_p_rate


######################################## DUTCH ROLL ##############################################

def dutch_roll(t0=60.1*60+5, deltat=28, plot_input=False, plot_output=False):
    # input -> rudder deflection
    u_dr, utime_dr = getInput(ruddertab, timetab, t0, deltat)

    # output -> yaw rate
    u_dr_y, utime_dr_y = getInput(yaw_rate_tab, timetab, t0, deltat)

    #       -> roll rate
    u_dr_r, utime_dr_r = getInput(roll_rate_tab, timetab, t0, deltat)

    ####plottingData####
    if plot_input == True:
        plottingData(utime_dr, u_dr, name="rudder_def_Dutch_Roll", title="Dutch Roll Input", variable="${\delta}_r$",
                 unit="deg", mins=False)

    if plot_output == True:
        plottingData(utime_dr_y, u_dr_y, name="yaw_rate_Dutch_Roll", title="Dutch Roll -> Yaw Rate", variable="r", unit="deg",
                 mins=False)
        plottingData(utime_dr_r, u_dr_r, name="roll_rate_Dutch_Roll", title="Dutch Roll -> Roll Rate",
                 variable="$p$", unit="deg/s", mins=False)

    return t0, deltat, utime_dr, u_dr, u_dr_y, u_dr_r


######################################## DUTCH ROLL YD ###########################################

def dutch_roll_yd(t0=60.95*60+5, deltat=19, plot_input=False, plot_output=False):
    # input -> rudder deflection
    u_dr_yd, utime_dr_yd = getInput(ruddertab, timetab, t0, deltat)

    # output -> yaw rate
    u_dr_yd_y, utime_dr_yd_y = getInput(yaw_rate_tab, timetab, t0, deltat)

    #       -> roll rate
    u_dr_yd_r, utime_dr_yd_r = getInput(roll_rate_tab, timetab, t0, deltat)

    ####plottingData####
    if plot_input == True:
        plottingData(utime_dr_yd, u_dr_yd, name="rudder_def_Dutch_Roll_YD", title="Dutch Roll YD Input", variable="${\delta}_r$",
                 unit="deg", mins=False)

    if plot_output == True:
        plottingData(utime_dr_yd_y, u_dr_yd_y, name="yaw_rate_Dutch_Roll_YD", title="Dutch Roll YD -> Yaw Rate", variable="r", unit="deg",
                 mins=False)
        plottingData(utime_dr_yd_r, u_dr_yd_r, name="roll_rate_Dutch_Roll_YD", title="Dutch Roll YD -> Roll Rate",
                 variable="$p$", unit="deg/s", mins=False)

    return t0, deltat, utime_dr_yd, u_dr_yd, u_dr_yd_y, u_dr_yd_r


######################################## APERIODIC ROLL ##########################################

def aperiodic_roll(t0=57.0*60, deltat=60, plot_input=False, plot_output=False):
    # input -> aileron deflection
    u_ar, utime_ar = getInput(ailerontab, timetab, t0, deltat)

    # output -> roll
    u_ar_r, utime_ar_r = getInput(roll_tab, timetab, t0, deltat)

    #       -> roll rate
    u_ar_r_rate, utime_ar_r_rate = getInput(roll_rate_tab, timetab, t0, deltat)

    ####plottingData####
    if plot_input == True:
        plottingData(utime_ar, u_ar, name="aileron_def_Aperiodic_Roll", title="Aperiodic Roll Input", variable="${\delta}_a$",
                 unit="deg", mins=False)

    if plot_output == True:
        plottingData(utime_ar_r, u_ar_r, name="roll_Aperiodic_Roll", title="Aperiodic Roll -> Roll", variable="$\phi$", unit="deg",
                 mins=False)
        plottingData(utime_ar_r_rate, u_ar_r_rate, name="roll_rate_Aperiodic_Roll", title="Aperiodic Roll -> Roll Rate",
                 variable="$p$", unit="deg/s", mins=False)

    return t0, deltat, utime_ar, u_ar, u_ar_r, u_ar_r_rate


######################################## SPIRAL ###############################################

def spiral(t0=3746, deltat=50, plot_input=False, plot_output=False):
    # input -> rudder deflection
    u_spi, utime_spi = getInput(ruddertab, timetab, t0, deltat)

    # output -> roll
    u_spi_r, utime_spi_r = getInput(roll_tab, timetab, t0, deltat)

    #       -> yaw rate
    u_spi_y, utime_spi_y = getInput(yaw_rate_tab, timetab, t0, deltat)

    ####plottingData####
    if plot_input == True:
        plottingData(utime_spi, u_spi, name="rudder_def_Spiral", title="Spiral Input", variable="${\delta}_r$",
                 unit="deg", mins=False)

    if plot_output == True:
        plottingData(utime_spi_r, u_spi_r, name="roll_Spiral", title="Spiral -> Roll", variable="$\phi$", unit="deg",
                 mins=False)
        plottingData(utime_spi_y, u_spi_y, name="yaw_rate_Spiral", title="Spiral -> Yaw Rate",
                 variable="$r$", unit="deg/s", mins=False)

    return t0, deltat, utime_spi, u_spi, u_spi_r, u_spi_y

#++++++++++++++++++++++++++++++++++++++ Plotting +++++++++++++++++++++++++++++++++++++++++++++++++++

# print(phugoid(plot_input=True, plot_output=True))
# print(short_period(plot_input=True, plot_output=True))
# print(dutch_roll(plot_input=True, plot_output=True))
# print(dutch_roll_yd(plot_input=True, plot_output=True))
# print(aperiodic_roll(plot_input=True, plot_output=True))
# print(spiral(plot_input=True, plot_output=True))


