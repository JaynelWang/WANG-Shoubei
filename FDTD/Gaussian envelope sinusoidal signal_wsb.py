# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:37:39 2016

@author: Louie
"""

import math

# Definition of constants for coordinates
RADIUS = 3 * pow(10, -2)            # m, Radius
RADIUS_ARRAY = 2 * pow(10, -2)      # m, Radius of phase array
DELTA_X = 0.1 * pow(10, -3)         # m, Delta_x
DELTA_R = 0.1 * pow(10, -3)         # m, Delta_r

DEPTH = RADIUS - math.sqrt(RADIUS ** 2 - RADIUS_ARRAY ** 2)         # m, The depth of phase array
INDEX_X = math.floor(DEPTH / DELTA_X) + 1                           # The index number of x array


# Definition of constants for input signal
PI = math.pi
DELTA_T = 10 * pow(10, -9)          # s, Delta_t

F_S = pow(10, 8)                    # Sample frequency
SAMPLES = 801                       # Sample number
n = 1
AMPL = n * pow(10, 6)               # Pa, Maximum amplitude
FREQ = 1 * pow(10, 6)               # Hz, Frequency
A_SIN = n * 2.6

MU_T = 4 * pow(10, -6)
SIGMA_T = A_SIN / (math.sqrt(2 * PI) * AMPL)


# Generate the array of (x,t) pairs
def func_gen_xr():
    array_x_r_l = [[0, 0]] * INDEX_X
    for i in range(INDEX_X):
        array_temp = [0, 0]
        r_raw = math.sqrt(RADIUS ** 2 - (RADIUS - i * DELTA_X) ** 2)
        r_pre = round(r_raw / DELTA_R) * DELTA_R
        array_temp[0] = DELTA_X * i
        array_temp[1] = r_pre
        # array_temp[0] = DELTA_X * i * pow(10, 3)          # If want data be stored as mm unit
        # array_temp[1] = r_pre * pow(10, 3)                # If want data be stored as mm unit
        array_x_r_l[i] = array_temp

    return array_x_r_l


# Generate the array of input P
def func_gen_input():
    array_p_t_l = [[0, 0]] * SAMPLES
    for j in range(SAMPLES):
        array_temp = [0, 0]
        t = DELTA_T * j
        fsin = A_SIN * math.sin(2 * PI * FREQ * (t - 0.75 * 1e-6))
        fg = (1 / (SIGMA_T * math.sqrt(2*PI))) * math.exp(-((t - MU_T)**2) / (2 * SIGMA_T**2))
        f = fsin * fg
        array_temp[0] = t
        array_temp[1] = f
        array_p_t_l[j] = array_temp

    return array_p_t_l


# Main function
array_x_r = func_gen_xr()           # The array of the (x,r) pairs
# print(array_x_r)
array_p_t = func_gen_input()        # The array of input signal
