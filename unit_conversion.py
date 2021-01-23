####### Converting text files from imperial to SI units ####

import numpy as np 
import csv
import os

### converting knots, lbs, g, ft, ft/min, deg C, lbs/hr, psi
### knots = 0.5144444444444444 m/s
### lbs = 0.45359237 kgs
### g = 9.80665 m/s^2
### ft = 0.3048 m
### ft/min = 0.00508 m/s
### deg C = + 273.15 K
### lbs/hr = 0.0001259979 kg/s
### psi = 6894.75729 pas

def conversion(filename):

    
    
    f = open(filename, "r")
    n = open(filename[:-4]+"SI.txt", "w")
    lines = f.readlines()

    unit = lines[1]
    n.write(lines[0])
    if unit == 'g\n':
        n.write('m/s^2\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())*9.80665)+"\n")
    if unit == 'lbs\n':
        n.write('kg\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())*0.45359237)+"\n")
    if unit == 'knots\n':
        n.write('m/s\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())*0.5144444444444444)+"\n")
    if unit == 'ft\n':
        n.write('m\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())*0.3048)+"\n")
    if unit == 'ft/min\n':
        n.write('m/s\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())*0.00508)+"\n")
    if unit == 'deg C\n':
        n.write('deg K\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())+ 273.15)+"\n")
    if unit == 'lbs/hr\n':
        n.write('kg/s\n')
        for i in range(len(lines)-2):
            n.write(str(float(lines[i+2].strip())*0.0001259979)+"\n")
    if unit == 'psi\n':
        n.write('Pa\n')
        for i in range(len(lines) - 2):
            n.write(str(float(lines[i + 2].strip())*6894.75729)+"\n")
    if unit!="knots\n" and unit!="lbs\n" and unit!="g\n" and unit!="ft\n" and unit!="ft/min\n" and unit!="deg C\n" and unit!="lbs/hr\n" and unit!="psi\n":
        for i in range(len(lines)-1):
            n.write(str(lines[i+1]))

    n.close()

for filename in os.listdir('Data_new'):
    name = "Data_new\\" + filename
    if "SI" in filename:
        os.remove(name)

for filename in os.listdir('Data_new'):
    name = "Data_new\\" + filename
    conversion(name)