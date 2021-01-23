from stationnary_functions import*
import matplotlib.pyplot as plt
import numpy as np
import math as m
import time as tt
from decimal import Decimal

#Choose Data
file='ref' #reference
#file=fd #flight test
#=======Data extraction======#
data=np.genfromtxt("stationnary1_"+str(file)+".csv", delimiter=",")
data=data[:,1:]

time=data[:,0].T
hp=data[:,1].T
VIAS=data[:,2].T
alpha=data[:,3].T
FFL=data[:,4].T
FFR=data[:,5].T
fuel=data[:,6].T
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

Re=Reynolds(np.average(T),np.average(rho),np.average(Vt)) #CHECK assumption
#=======Aerodynamic coefficients======#
CL=[]
CL2=[]
CD=[]
for i in range(0,len(Vt)):
    a,b=aero_coefficient(W[i],Th[i],rho[i],Vt[i],alpha[i])
    CL.append(a)
    CL2.append(a**2)
    CD.append(b)


#=======Fitting data with theoretical curves======#
CLa, CL0=np.polyfit(alpha,CL,1)
grad, CD0=np.polyfit(CL2,CD,1)
CDtheory=np.linspace(CD0,0.08,500)
CL2theory=[]
for i in range(0,len(CDtheory)):
    a=m.sqrt((CDtheory[i]-CD0)/grad)
    CL2theory.append(a)
e=1/(grad*m.pi*AR)

print('Lift slope: '+str(CLa*180/m.pi))
print('Oswald factor: '+str(e))
print('Zero lift drag: '+str(CD0))

#=======Plotting======#
fig = plt.figure()
ax = fig.gca()

#=======Lift/alpha plot======#
# plt.scatter(alpha,CL,color='black', marker='o') #CHECK for linear part
# plt.plot([-CL0/CLa,12],[0,CLa*12+CL0],color='#0080FF')
# plt.xlim(-1,12)
# plt.ylim(0,1)
# ax.set_xticks(np.arange(-1, 13, 1))
# ax.set_yticks(np.arange(0, 1, 0.1))
# ax.grid(which='both')
# plt.xlabel('$\\alpha$ [deg]',fontsize=16)
# plt.ylabel('$C_L$ [-]',fontsize=14)
# plt.title('Cessna Citation II Lift versus Angle of attack', fontsize=14)
# ax.text(-0.7, 0.84, '$R_e$='+str('{:.3E}'.format(Decimal(Re)))+'\n'
#                     'M: '+str(round(min(M),3))+'-' +str(round(max(M),3))+'\n'
#                     'Clean configuration',style='italic',bbox={'facecolor': 'none', 'edgecolor':'black'})
# plt.savefig('Plots\CLvsalpha.png')

#=======Lift/drag plot======#
# plt.scatter(CD,CL,color='black', marker='o') #CHECK for linear part
# plt.plot(CDtheory,CL2theory,color='#0080FF')
# plt.xlim(0.02,0.08)
# plt.ylim(0,1)
# ax.set_xticks(np.arange(0.02, 0.09, 0.01))
# ax.set_xticks(np.arange(0.02, 0.08, 0.005), minor='True')
# ax.set_yticks(np.arange(0, 1.1, 0.1))
# ax.grid(which='both')
# plt.xlabel('$C_D$ [-]',fontsize=16)
# plt.ylabel('$C_L$ [-]',fontsize=14)
# plt.title('Cessna Citation II Drag polar', fontsize=14)
# ax.text(0.022, 0.84, '$R_e$='+str('{:.3E}'.format(Decimal(Re)))+'\n'
#                     'M: '+str(round(min(M),3))+'-' +str(round(max(M),3))+'\n'
#                     'Clean configuration',style='italic',bbox={'facecolor': 'none', 'edgecolor':'black'})
# plt.savefig('Plots\Dragpolar.png')

plt.show()
