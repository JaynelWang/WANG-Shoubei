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

#generate P0 and P1
P0 = np.zeros([513,257])
P1 = np.zeros([513,257])
array_coor = 77 * []
#array_x_r(unit:m),array_coor(unit:0.1mm)
for i in range(77):
    array_coor.append(np.array(10000 * np.array(array_x_r[i]),int))
    P1[array_coor[i][0],array_coor[i][1]] = array_p_t[1][1]

#initialization
c = 1445                   #speed of sound
rou = 921                  #tissue density
f = 1e6                    #frequency of input signal
B_A = 10.3                 #parameter:B/A
beta = 1 + 0.5 * B_A
alpha_2 = 7 / 1e6 * f      #attenuation coefficient(unit:Np/m/MHz)
delta = 2 * c**2 * alpha_2 / (2 * np.pi * f)**2

alpha = -0.3               #dp/dx=dp/dr=alpha(boundary condition)
delta_t = 1e-8             #unit:s
delta_x = 1e-4             #unit:m
delta_r = 1e-4             #unit:m

a = 1 / (c**2 * delta_t**2) - beta / (rou * c**4 * delta_t**2)
b = 1 / (12 * delta_x**2)
c1 = -1 / (12 * delta_r**2)
d1 = 4 / (3 * delta_r**2)
e1 = 4 / (3 * delta_r**2)
f1 = -1 / (12 * delta_r**2)
g = -5 / (2 * delta_r**2) - 5 / (2 * delta_x**2) + 2 / (c**2 * delta_t**2) - \
    2 * beta / (rou * c**4 * delta_t**2)
h = -a - 23 * delta / (8 * c**4 * delta_t**3)
k = delta / (8 * c**4 * delta_t**3)

def c2(r):
        return -1 / (12 * r * delta_r**2)

def d2(r):
        return 2 / (3 * r * delta_r**2)

def e2(r):
        return -2 / (3 * r * delta_r**2)
    
def f2(r):
        return 1 / (12 * r * delta_r**2)

P = 801 * []          #acoustic pressure matrixes are stored in list P
P.append(P0)
P.append(P1)

#generate matrix A
column0 = [g,d1,c1] + 254 * [0]
column1 = [e1,g,d1,c1] + 253 * [0]
column255 = 253 * [0] + [f1,e1,g,d1]
column256 = 254 * [0] + [f1,e1,g]
column = [column0,column1]
zero_up = 0
zero_bottom = 252
while zero_bottom >= 0:
    temp = zero_up * [0] + [f1,e1,g,d1,c1] + zero_bottom * [0]
    column.append(temp)
    zero_up += 1
    zero_bottom -= 1
column.append(column255)
column.append(column256)
A = np.vstack((column[0],column[1]))
for i in range(2,257):
    A = np.vstack((A,column[i]))
A = A.transpose((1,0))

#generate matrix B
column0 = 257 * [0]
column1 = [e2(1),0,d2(1),c2(1)] + 253 * [0]
column255 = 253 * [0] + [f2(255),e2(255),0,d2(255)]
column256 = 254 * [0] + [f2(256),e2(256),0]
column = [column0,column1]
zero_up = 0
zero_bottom = 252
while zero_bottom >= 0:
    temp = zero_up * [0] + [f2(zero_up+2),e2(zero_up+2),0,d2(zero_up+2),c2(zero_up+2)] + zero_bottom * [0]
    column.append(temp)
    zero_up += 1
    zero_bottom -= 1
column.append(column255)
column.append(column256)
B = np.vstack((column[0],column[1]))
for i in range(2,257):
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
Pim_1 = P1[:,0] + 2 * alpha * delta_r
Pim_2 = P1[:,0] + 4 * alpha * delta_r
Pip_1 = P1[:,256] + 2 * alpha * delta_r
Pip_2 = P1[:,256] + 4 * alpha * delta_r
Pjm_1 = P1[0,:] + 2 * alpha * delta_x
Pjm_2 = P1[0,:] + 4 * alpha * delta_x
Pjp_1 = P1[512,:] + 2 * alpha * delta_x
Pjp_2 = P1[512,:] + 4 * alpha * delta_x

D1 = e1 * np.column_stack((Pim_1,np.zeros((513,256))))
D2 = f1 * np.column_stack((Pim_2,np.zeros((513,256))))
D3 = (f1 + f2(1)) * np.column_stack((np.zeros((513,1)),Pim_1,np.zeros((513,255))))
D4 = (c1 + c2(256)) * np.column_stack((np.zeros((513,256)),Pip_2))
D5 = (d1 + d2(256)) * np.column_stack((np.zeros((513,256)),Pip_1))
D6 = (c1 + c2(255)) * np.column_stack((np.zeros((513,255)),Pip_1,np.zeros((513,1))))
D7 = 16 * b * np.vstack((Pjm_1,np.zeros((512,257))))
D8 = -b * np.vstack((Pjm_2,np.zeros((512,257))))
D9 = -b * np.vstack((np.zeros((1,257)),Pjm_1,np.zeros((511,257))))
D10 = -b * np.vstack((np.zeros((511,257)),Pjp_1,np.zeros((1,257))))
D11 = 16 * b * np.vstack((np.zeros((512,257)),Pjp_1))
D12 = -b * np.vstack((np.zeros((512,257)),Pjp_2))

D = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12

#calculation itenerations
#with input signal
#calculate P2 ~ P5 (n has to be >=5)
for i in range(2,6):
    P2 = (np.dot(P1,A) + np.dot(P1,B) + np.dot(C,P1) - a * P0 + D) / a
    for j in range(77):
        P2[array_coor[j][0],array_coor[j][1]] = array_p_t[i][1]
    P.append(P2)
    P0 = P1
    P1 = P2
    
#recalculate matrix A
g = -5 / (2 * delta_r**2) - 5 / (2 * delta_x**2) + 2 / (c**2 * delta_t**2) + \
    -2 * beta / (rou * c**4 * delta_t**2) + 3 * delta / (4 * c**4 * delta_t**3)
column0 = [g,d1,c1] + 254 * [0]
column1 = [e1,g,d1,c1] + 253 * [0]
column255 = 253 * [0] + [f1,e1,g,d1]
column256 = 254 * [0] + [f1,e1,g]
column = [column0,column1]
zero_up = 0
zero_bottom = 252
while zero_bottom >= 0:
    temp = zero_up * [0] + [f1,e1,g,d1,c1] + zero_bottom * [0]
    column.append(temp)
    zero_up += 1
    zero_bottom -= 1
column.append(column255)
column.append(column256)
A = np.vstack((column[0],column[1]))
for i in range(2,257):
    A = np.vstack((A,column[i]))
A = A.transpose((1,0))

#calculate P6 and after
P0 = P[0]
P1 = P[1]
P2 = P[2]
P3 = P[3]
P4 = P[4]
P5 = P[5]

for i in range(6,801):
    P6 = (np.dot(P5,A) + np.dot(P5,B) + np.dot(C,P5) + h * P4 + D + \
         k * (-P0 + 8 * P1 -24 * P2 + 34 * P3)) / a
    for j in range(77):
        P6[array_coor[j][0],array_coor[j][1]] = array_p_t[i][1]
    P.append(P6)
    P0 = P1
    P1 = P2
    P2 = P3
    P3 = P4
    P4 = P5
    P5 = P6

#without input signal
#i = 800
#while (i > 0):
#    P6 = (np.dot(P5,A) + np.dot(P5,B) + np.dot(C,P5) + h * P4 + D + \
#         k * (-P0 + 8 * P1 -24 * P2 + 34 * P3)) / a
#    P.append(P6)
#    P0 = P1
#    P1 = P2
#    P2 = P3
#    P3 = P4
#    P4 = P5
#    P5 = P6
#    i -= 1
#
#del P[0:6]
sio.savemat('D:/lab_Jaynel/Lab/FDTD/nonlinear_non-repeating_1MPa_1MHz_fat_1.mat',{'P':P})


#display animation
#data = sio.loadmat('D:/lab_Jaynel/Lab/FDTD/Output data_1.mat')
#P = data['P']
#fig = plt.figure()
#
#def init_animation():
#    fig = plt.imshow(P[0])
#    return fig,
#    
#def animate(i): 
#    fig = plt.imshow(P[i])
#    return fig,
#
#anim = matplotlib.animation.FuncAnimation(fig,animate,init_func=init_animation,frames = 801,interval=0.01)
#plt.show()