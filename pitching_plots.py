from stationnary_functions import*
import matplotlib.pyplot as plt
import numpy as np
import math as m

#Choose Data
file='ref' #reference
#file=fd #flight test
#=======Data extraction======#
data=np.genfromtxt("stationnary2_"+str(file)+".csv", delimiter=",")
data=data[:,1:]

time=data[:,0].T
hp=data[:,1].T
VIAS=data[:,2].T
alpha=data[:,3].T
de=data[:,4].T
Fe=data[:,6].T
FFL=data[:,7].T
FFR=data[:,8].T
fuel=data[:,9].T
TAT=data[:,-1].T

#=======Weight and Thrust extraction======#
if file=='ref':
    m_start=4050+9165 #[lbs]
    m_pax=95+92+74+66+61+75+78+86+68 #[kg]
else:
    m_start=2640+9165 #[lbs]
    m_pax=95+102+89+82+66+81+69+85+96 #[kg]
W=[]
for i in range(0,len(fuel)):
    a=((m_start-fuel[i])*0.453582+m_pax)*g0
    W.append(a)

thrust(hp,VIAS,FFL,FFR,TAT)
tt.sleep(3)
tdata=np.loadtxt( 'thrust.dat' )
Th=[]
for i in range(0,len(tdata[:,0])):
    t=tdata[i,0]+tdata[i,1]
    Th.append(t)

FF=np.ones(len(VIAS))*0.048/0.453582*3600
thrust(hp,VIAS,FF,FF,TAT)
tt.sleep(3)
tdata2=np.loadtxt( 'thrust_standard.dat' )
Ths=[]
for i in range(0,len(tdata2[:,0])):
    t=tdata2[i,0]+tdata2[i,1]
    Ths.append(t)

#========Elevator effectiveness========#
#CHECK assumptions
if file=='ref':
    Ve, Vt, rho, T, M=speed_correction(161,5730,5)
    CN,CT=aero_coefficient((m_start-881)*0.453582*g0,0,rho,Vt,5.3)
    delta_xcg=-0.04130226605496423
    Cmd=elevator_effectiveness(CN,delta_xcg,-0.5*m.pi/180)
else:
    Ve, Vt, rho, T, M=speed_correction(156,18360,-11.2)
    CN,CT=aero_coefficient((m_start-940)*0.453582*g0,0,rho,Vt,5.2)
    delta_xcg=-0.05073 #flight test
    Cmd=elevator_effectiveness(CN,delta_xcg,-0.6*m.pi/180)

print('Elevator effectiveness: '+str(Cmd))

#========Elevatro trim and contorl force curves========#

#=======Airspeed and Atmospheric properties======#
Ve=[]
Vt=[]
rho=[]
T=[]
M=[]
for i in range(0,len(VIAS)):
    a,b,c,d,e=speed_correction(VIAS[i],hp[i],TAT[i])
    Ve.append(a)
    Vt.append(b)
    rho.append(c)
    T.append(d)
    M.append(e)

#=======Weight correction======#
Ve_tilde=[]
Fe_star=[]
for i in range(0,len(Ve)):
    a,b=weight_correction(Ve[i],Fe[i],W[i])
    Ve_tilde.append(a)
    Fe_star.append(b)

#=======Thrust correction======#
de_star=[]
for i in range(0,len(Ve)):
    Tc=thrust_coefficient(Th[i],rho[i],Vt[i])
    Tcs=thrust_coefficient(Ths[i],rho0,Vt[i])
    a=thrust_correction(de[i],Cmd,Tcs,Tc)
    de_star.append(a)

print('Longitudinal Stability: '+str(longitudial_stability(alpha,de_star,Cmd)))
#CHECK if it should be corrected de or not. same for Cmd

#=======Fitting data with theoretical curves======#
Ve_inv=[]
for i in range(0,len(Ve_tilde)):
    a=1/(Ve_tilde[i]**2)
    Ve_inv.append(a)
grad, y_int=np.polyfit(Ve_inv,de_star,1)
Ve_theory=np.linspace(60,105,500)
de_theory=[]
for i in range(0,len(Ve_theory)):
    a=grad/(Ve_theory[i]**2)+y_int
    de_theory.append(a)

c1,c2,c3=np.polyfit(Ve_tilde,Fe_star,2)
Fe_theory=[]
for i in range(0,len(Ve_theory)):
    a=c1*Ve_theory[i]**2+c2*Ve_theory[i]+c3
    Fe_theory.append(a)

#=======Plotting======#
fig = plt.figure()
ax = fig.gca()

#=======Elevator trim curve======# #Mirror or not???, Radians???

# plt.scatter(Ve_tilde,de_star,color='black', marker='o')
# plt.plot(Ve_theory,de_theory,color='#0080FF')
# plt.xlim(60,100)
# plt.ylim(-2,1.5)
# ax.set_xticks(np.arange(60, 100, 5))
# ax.set_yticks(np.arange(-2, 2.0, 0.5))
# plt.gca().invert_yaxis()
# ax.grid(which='both')
# plt.xlabel('$V_e$ [m/s]',fontsize=16)
# plt.ylabel('$\\delta_e$ [deg]',fontsize=14)
# plt.title('Cessna Citation II Elevator Trim Curve', fontsize=14)
# plt.savefig('Plots\Trimcurve.jpeg')

#=======Elevator trim curve======# #Mirror or not???

plt.scatter(Ve_tilde,Fe_star,color='black', marker='o')
plt.plot(Ve_theory,Fe_theory,color='#0080FF')
plt.xlim(60,105)
plt.ylim(-60,100)
ax.set_xticks(np.arange(60, 105, 5))
ax.set_yticks(np.arange(-60, 120, 20))
plt.gca().invert_yaxis()
ax.grid(which='both')
plt.xlabel('$V_e$ [m/s]',fontsize=16)
plt.ylabel('$F_e$ [N]',fontsize=14)
plt.title('Cessna Citation II Elevator Control Force Curve', fontsize=14)
plt.savefig('Plots\controlforcecurve.jpeg')


plt.show()
