# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy.io as sio  
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import animation

#generate P0-P5.(size of P: 5.12 * 2.56 cm^2)
P0 = np.zeros([513,256])
P1 = np.zeros([513,256])
P2 = np.zeros([513,256])
P3 = np.zeros([513,256])
P4 = np.zeros([513,256])
P5 = np.zeros([513,256])
array_coor = 77 * []
#array_x_r(unit:m),array_coor(unit:0.1mm)
array_coor = 76 * []       #76 transducers on a semi-phased-array
del array_x_r[0]           #delete coordinate(0,0)(j=0)
for i in range(76):
    array_coor.append(np.array(10000 * np.array(array_x_r[i]),int))
    P5[array_coor[i][0],array_coor[i][1]] = array_p_t[1][1]

#initialization
c = 1445                   #speed of sound
rou = 921                  #tissue density
f = 1e6                    #frequency of input signal
B_A = 10.3                 #parameter:B/A
beta = 1 + 0.5 * B_A
alpha_2 = 7 / 1e6 * f      #attenuation coefficient(unit:Np/m/MHz)
delta = 2 * c**3 * alpha_2 / (2 * np.pi * f)**2              #diffusivity of sound

alpha = -0.3               #dp/dx=dp/dr=alpha(boundary condition)
delta_t = 1e-8             #unit:s
delta_x = 1e-4             #unit:m
delta_r = 1e-4             #unit:m
#coefficients in formula
a = 1 / (c**2 * delta_t**2)
b = 1 / (12 * delta_x**2)
c1 = -1 / (12 * delta_r**2)
d1 = 4 / (3 * delta_r**2)
e1 = 4 / (3 * delta_r**2)
f1 = -1 / (12 * delta_r**2)
g = -5 / (2 * delta_r**2) - 5 / (2 * delta_x**2) + 2 / (c**2 * delta_t**2) + \
    3 * delta / (4 * c**4 * delta_t**3)
h = -a - 23 * delta / (8 * c**4 * delta_t**3)
k = delta / (8 * c**4 * delta_t**3)
l = beta / (2 * rou * c**4 * delta_t**2)

def c2(r):
        return -1 / (12 * r * delta_r**2)

def d2(r):
        return 2 / (3 * r * delta_r**2)

def e2(r):
        return -2 / (3 * r * delta_r**2)
    
def f2(r):
        return 1 / (12 * r * delta_r**2)

P = 801 * []          #acoustic pressure matrixes are stored in list P(length:8us)
P.append(P5)

#generate matrix A
column0 = [g,d1,c1] + 253 * [0]
column1 = [e1,g,d1,c1] + 252 * [0]
column254 = 252 * [0] + [f1,e1,g,d1]
column255 = 253 * [0] + [f1,e1,g]
column = [column0,column1]
zero_up = 0
zero_bottom = 251
while zero_bottom >= 0:
    temp = zero_up * [0] + [f1,e1,g,d1,c1] + zero_bottom * [0]
    column.append(temp)
    zero_up += 1
    zero_bottom -= 1
column.append(column254)
column.append(column255)
A = np.vstack((column[0],column[1]))
for i in range(2,256):
    A = np.vstack((A,column[i]))
A = A.transpose((1,0))

#generate matrix B
column0 = [0,d2(1),c2(1)] + 253 * [0]
column1 = [e2(2),0,d2(2),c2(2)] + 252 * [0]
column254 = 252 * [0] + [f2(255),e2(255),0,d2(255)]
column155 = 253 * [0] + [f2(256),e2(256),0]
column = [column0,column1]
zero_up = 0
zero_bottom = 251
while zero_bottom >= 0:
    temp = zero_up * [0] + [f2(254-zero_bottom),e2(254-zero_bottom),0,d2(254-zero_bottom),c2(254-zero_bottom)] + zero_bottom * [0]
    column.append(temp)
    zero_up += 1
    zero_bottom -= 1  
column.append(column254)
column.append(column255)
B = np.vstack((column[0],column[1]))
for i in range(2,256):
    B = np.vstack((B,column[i]))
B = B.transpose((1,0))

#generate matrix C
row0 = [0,16*b,-b] + 510 * [0]
row1 = [16*b,0,16*b,-b] + 509 * [0]
row511 = 509 * [0] + [-b,16*b,0,16*b]
row512 = 510 * [0] + [-b,16*b,0]
row = [row0,row1]
zero_left = 0
zero_right = 508
while zero_right >= 0:
    temp = zero_left * [0] + [-b,16*b,0,16*b,-b] + zero_right * [0]
    row.append(temp)
    zero_left += 1
    zero_right -= 1
row.append(row511)
row.append(row512)
C = np.vstack((row[0],row[1]))
for i in range(2,513):
    C = np.vstack((C,row[i]))

#generate matrix D
def calculateD():
    Pim_0 = P5[:,0] + alpha * delta_r
    Pim_1 = P5[:,0] + 2 * alpha * delta_r
    Pip_1 = P5[:,255] + alpha * delta_r
    Pip_2 = P5[:,255] + 2 * alpha * delta_r
    Pjm_1 = P5[0,:] + alpha * delta_x
    Pjm_2 = P5[0,:] + 2 * alpha * delta_x
    Pjp_1 = P5[512,:] + alpha * delta_x
    Pjp_2 = P5[512,:] + 2 * alpha * delta_x

    D1 = (e1 + e2(1)) * np.column_stack((Pim_0,np.zeros((513,255))))
    D2 = (f1 + f2(1)) * np.column_stack((Pim_1,np.zeros((513,255))))
    D3 = (f1 + f2(2)) * np.column_stack((np.zeros((513,1)),Pim_0,np.zeros((513,254))))
    D4 = (c1 + c2(256)) * np.column_stack((np.zeros((513,255)),Pip_2))
    D5 = (d1 + d2(256)) * np.column_stack((np.zeros((513,255)),Pip_1))
    D6 = (c1 + c2(255)) * np.column_stack((np.zeros((513,254)),Pip_1,np.zeros((513,1))))
    D7 = 16 * b * np.vstack((Pjm_1,np.zeros((512,256))))
    D8 = -b * np.vstack((Pjm_2,np.zeros((512,256))))
    D9 = -b * np.vstack((np.zeros((1,256)),Pjm_1,np.zeros((511,256))))
    D10 = -b * np.vstack((np.zeros((511,256)),Pjp_1,np.zeros((1,256))))
    D11 = 16 * b * np.vstack((np.zeros((512,256)),Pjp_1))
    D12 = -b * np.vstack((np.zeros((512,256)),Pjp_2))

    D = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12
    return D

#calculation itenerations
#with input signal(within the first 8us)
#calculate P6 and after
for i in range(1,801):
    D = calculateD()
    P6 = (np.dot(P5,A+B) + np.dot(C,P5) + D + h * P4 + \
         k * (-P0 + 8 * P1 -24 * P2 + 34 * P3) + \
         l * (4*P5 * (2*P5 - 5*P4 + 4*P3 - P2) + (3*P5 - 4*P4 + P3)**2)) / a
    for j in range(76):
        P6[array_coor[j][0],array_coor[j][1]] = array_p_t[i][1]
    P.append(P6)
    P0 = P1
    P1 = P2
    P2 = P3
    P3 = P4
    P4 = P5
    P5 = P6

#without input signal(after 8us)
"""
i = 800
while (i > 0):
    D = calculateD()
    P6 = (np.dot(P5,A+B) + np.dot(C,P5) + D + h * P4 + \
         k * (-P0 + 8 * P1 -24 * P2 + 34 * P3) + \
         l * (4*P5 * (2*P5 - 5*P4 + 4*P3 - P2) + (3*P5 - 4*P4 + P3)**2)) / a
    P.append(P6)
    P0 = P1
    P1 = P2
    P2 = P3
    P3 = P4
    P4 = P5
    P5 = P6
    i -= 1
del P[0]
"""
#save the output data in ".mat" format
sio.savemat('D:/FDTD/nonlinear_non-repeating_1MPa_1MHz_fat_1.mat',{'P':P})

#display animation
#If animation cannot be displayed in Spyder, try running it in Ipython.
"""
data = sio.loadmat('D:/FDTD/output data.mat')
P = data['P']
fig = plt.figure()

def init_animation():
    fig = plt.imshow(P[0])
    return fig,
    
def animate(i): 
    fig = plt.imshow(P[i])
    return fig,

anim = matplotlib.animation.FuncAnimation(fig,animate,init_func=init_animation,frames = 801,interval=0.01)
plt.show()
"""