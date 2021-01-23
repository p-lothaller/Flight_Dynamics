import numpy as np
from massbalance import *


#FOR GRAPH
nsteps= 100
x = np.linspace(10,timelst[-1],nsteps)
y = []


for i in range(len(x)):
    if 52*60 < x[i] < 53*60:
        print('hit')
        y.append(cg(x[i], True))
    else:
        y.append(cg(x[i],False))

# y1 = cg(x[x<52*60],False)
# yshift =cg(x[52*60 < x],True)
# y = np.concatenate([y1,yshift])
# y[x > 53*60] = cg(x[53*60 < x],False)

plt.plot(x,y)
plt.grid(True)
plt.xlabel('Time [s]')
plt.ylabel('cg location [-]')
plt.show()


# for i in range(len(x)):
#     y.append(mass(x[i]))
# plt.plot(x,y)
# plt.grid(True)
# plt.xlabel('Time [s]')
# plt.ylabel('Mass [kg]')
# plt.show()

# for i in range(len(x)):
#     y.append(fuelmomentall(x[i]))
# plt.plot(x,y)
# plt.grid(True)
# plt.xlabel('Time [s]')
# plt.ylabel('Fuel Moment Arm ')
# plt.show()
#
# for i in range(len(x)):
#     y.append(fuelonboard(x[i]))
# plt.plot(x,y)
# plt.grid(True)
# plt.xlabel('Time [s]')
# plt.ylabel('Fuel [lb]')
# plt.show()
