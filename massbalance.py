import numpy as np
import matplotlib.pyplot as plt



#VARIABLE parameters
blockfuel = 2640 #Fuel at t=0      check consistent with others

#fuelusedlst = np.array([666, 555, 444, 333, 222, 111])  # Wf (Fuel Used)                USE OUR DATA
#timelst = np.array([1, 2, 3, 4, 5, 6])  # t values for Wf (Fuel Used)   USE OUR DATA

timelst = np.genfromtxt("Data_SI_correct/timeSI.txt",skip_header=2)
lhfusi = np.genfromtxt("Data_SI_correct/lh_engine_FUSI.txt",skip_header=2)
rhfusi = np.genfromtxt("Data_SI_correct/rh_engine_FUSI.txt",skip_header=2)
fuelusedlst = (lhfusi+rhfusi)/0.453592   #given in SI and converted for consistency


#Pax weights   [kg]    [these are SI units which are converted later to lbs for consistency]
Arun = 95
Hans = 102
Paul = 89
Vlad = 82
Max = 66
Mat = 81
Flori = 69
Phil = 85
Jack = 96



#NON-VARIABLE parameters
BEW = 9165 #Basic Empty Weight    lbs    (BS Mass and Balance Report)
BEWarm = 291.65 #Basic Empty Weight Moment Arm   inches from datum  (BS Mass and Balance Report)
seatsarm = np.array([131,131,214,214,251,251,288,288,170,170]) #Seat Moment Arms   inches from datum   (Note: discrepancy between BS and Assignment 214 vs 216)

#Fuel Moment Arm Initial Parameters
fuelmass = np.linspace(100,2800,28)
fuelmomentall = 100*np.array([298.16,591.18,879.08,1165.42,1448.4,1732.53,2014.8,2298.84,2581.92,2866.3,3150.18,3434.52,3718.52,4003.23,4287.76,4572.24,4856.56,5141.16,5425.64,5709.9,5994.04,6278.47,6562.82,6846.96,7131,7415.33,7699.6,7984.34])

#print(len(fuelmass),len(fuelmomentall)) #verification

foo = 1  # Fuel on (1) or off (0)
soo = 1  # Seats on (1) or off (0)


def interpolate(x1,x2,y1,y2,x):
    return (y2-y1)*(x-x1)/(x2-x1) + y1



def fuelonboard(t):
    for i in range(len(timelst)-1):
        if timelst[i] <= t <= timelst[i+1]:
            f_used = interpolate(timelst[i],timelst[i+1],fuelusedlst[i],fuelusedlst[i+1],t)
            return blockfuel - f_used
    if t <= timelst[0]:
        #print(t)
        print("t less than first recorded value. No fuel used yet")
        return blockfuel
    if t >= timelst[-1]:
        print("t specified greater than last recorded t-value (constant extrapolation)")
        return blockfuel - fuelusedlst[-1]

def fuelmoment(t):
    for l in range(len(fuelmomentall)-1):
        if fuelmass[l] <= fuelonboard(t) <= fuelmass[l+1]:
            return interpolate(fuelmass[l],fuelmass[l+1],fuelmomentall[l],fuelmomentall[l+1],fuelonboard(t))
    if fuelonboard(t) < fuelmass[0]:
        print("Fuel too low - Moment arm inaccurate")
    if fuelonboard(t) > fuelmass[-1]:
        print("Fuel beyond moment arm data")


def cg(t,cgmove): #cgmove = False if Jack is in seat 8, set cgmove = True if Jack has moved for cg shift
    seat1 = Arun
    seat2 = Hans
    seat3 = Vlad
    seat4 = Max
    seat5 = Mat
    seat6 = Flori
    seat7 = Phil
    seat8 = Jack
    seat9 = 0
    seat10 = Paul

    if cgmove == True:
        seat8 = 0
        seat9 = Jack


    seats = np.array([seat1,seat2,seat3,seat4,seat5,seat6,seat7,seat8,seat9,seat10])
    seats = seats/0.453592  #convert to lbs !!

    seatmoment = 0
    seatmass = 0
    for s in range(len(seats)):
        seatmoment += seats[s]*seatsarm[s]
        seatmass += seats[s]



    momentsum = soo*seatmoment + BEW*BEWarm + foo*fuelmoment(t)
    totalmass = soo*seatmass + BEW + foo*fuelonboard(t)
    x_cg_inch = momentsum/totalmass
    x_cg = (x_cg_inch - 261.56)*2.54/100  #x_cg from LEMAC   in [m] !!
    return x_cg #This is in [m] from LEMAC !! (not as frac/MAC)

def mass(t):
    seat1 = Arun
    seat2 = Hans
    seat3 = Vlad
    seat4 = Max
    seat5 = Mat
    seat6 = Flori
    seat7 = Phil
    seat8 = Jack
    seat9 = 0
    seat10 = Paul

    seats = np.array([seat1,seat2,seat3,seat4,seat5,seat6,seat7,seat8,seat9,seat10])
    seats = seats/0.453592  #convert to lbs !!

    seatmass = 0
    for s in range(len(seats)):
        seatmass += seats[s]

    totalmass = soo*seatmass + BEW + foo*fuelonboard(t)
    return totalmass*0.453592 #Total mass in [kg] !!!


#print(cg(52.5*60,False)-cg(52.5*60,True))
#print(100*(cg(10,False)*39.3701+261.56-291.65)/291.65,mass(10)/0.453592)
xt = np.array([1358,1520,1682,1800,1890,2130,2450,2555,2652,2748,2805,2855,2970,3220,3435,3519,3610,3760]) #time list, note: no cg shift


# #finding Report Requirements Bullet nr 1:
# for l in xt:
#     print(round((cg(l,False)*100/2.54+261.56),3))
#
# for l in xt:
#     print(round((cg(l,False)),4))
#
# for l in xt:
#     print(round(mass(l),3))
#
# #cg shift
# print(round((cg(3157,True)*100/2.54+261.56),3),round((cg(3157,True)),3),round(mass(3157),3))
# print(round(cg(3157,True),4))