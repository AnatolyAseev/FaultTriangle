#!/usr/bin/env python3
""" 
data.py is a powerhorse for the Fault Triangle tool.
It takes two .las files and converts them to numpy arrays using lasio Python library. 

Then it:
1. Assigns the shorter .las file as las1 and the longer .las file as las2
2. Extends las1 by connecting the missing part from las2, so as las1.depth == las2.depth
3. Converts GR to Vshale for both las1 and las2. It converts GR to Vshale using 5 most common equation: a.Linear, b.Larionov_Young,
c. Larionov_old, d.Steiber and e.Clavier
4. Estimates Shale Gouge Ratio (SGR) using chosen GR->Vshale conversion method and list of hypothetical
fault throw values. For the sake of proper plotting, the default fault throw is a thickness of the longest .las file, that means
throw = the last depth value - the first depth value.
5. Estimates a maximum burial depth of the fault - Zmax. Which is estimated as current depth + current throw.
6. Estimates critical across fault difference pressures using Yielding (2012) and Bretan(2003)
7. Estimates mercury-air permeability (Kfma) -> mercury-air threshold pressure (Pfma) -> threshold pressure for any pair of fluids (Pfhw).
By default light oil and water constants are used for converting Pfma to Pfhw (threshold pressure for hydrocarbons). All equations are
from Sperrevik et al(2002)
8. Having threshold pressures from all three methods we are converting them to the hydrocarbon column height (HCCH)
9. Then we are combining all the data in one dataframe for further plotting

"""

import numpy as np
from numpy import arange
import pandas as pd
import random
import lasio as las
import math
import itertools


def loader(path1, path2):
    """
    loads two .las files. Assigns shorter well as las1

    """
    las1_load=las.read(path1)
    las2_load=las.read(path2)
    las1=las1_load.data
    las2=las2_load.data
    if las1.shape[0] > las2.shape[0]:
        las1, las2=las2, las1
    else:
        las1, las2
    return las1, las2

def las1_extension(las1, las2):
    """ 
    extend the shorter .las which is las1 
    by attaching the lower part of the longer .las (las2) to the shorter .las (las1)
    """
    #how many more indexes in the longer .las
    diff=las2.shape[0] - las1.shape[0] 

    #step between indices in meters
    step=las1[1][0] - las1[0][0]

    # the total depth in meters for the shorter .las (las1)
    total=las1[-1,0] + (diff * step) 

    #making additional depth 
    las1depth_add=np.asarray([i + step for i in arange(las1[-1,0], total, step)])

    #making additional GR values
    las1GR_add=np.asarray([i for i in las2[las1.shape[0]:, 1]])

    #connecting additional depth and additional GR to the bottom of depth and GR arrays of the las1
    las1depth_true=np.append(las1[:,0], las1depth_add)
    las1GR_true=np.append(las1[:,1], las1GR_add)

    #stacking two arrays in 2D array: depth and GR
    las1=np.stack((las1depth_true, las1GR_true), axis = 1)
    return las1

def GR_to_Vshale(las1, las2):
    """

    add Vshale for every GR->Vshale convertion method to the numpy array

    """
    #Linear
    LGR1=(las1[:,1] - np.amin(las1[:,1])) / (np.amax(las1[:,1] - np.amin(las1[:,1])))
    LGR2 = (las2[:,1] - np.amin(las2[:,1])) / (np.amax(las2[:,1] - np.amin(las2[:,1])))
    
    #Larionov_Young Vshale
    LAR_YOUNG1=0.083 * (np.power(2, 3.7 * LGR1) - 1)
    LAR_YOUNG2=0.083 * (np.power(2, 3.7 * LGR2) - 1)
    
    #Larionov_Old Vshale
    LAR_OLD1=0.33 * (np.power(2, 2 * LGR1) - 1)
    LAR_OLD2 = 0.33 * (np.power(2, 2 * LGR2) - 1)
    
    #Steiber Vshale
    STEIBER1=LGR1 / (3 - 2 * LGR1)
    STEIBER2=LGR2 / (3 - 2 * LGR2)
    
    #Clavier Vshale
    CLAVIER1=1.7 - ((3.38 - (LGR1 + 0.7)**2))**0.5
    CLAVIER2=1.7 - ((3.38 - (LGR2 + 0.7)**2))**0.5
    
    #stack Vshale to las1 and las2
    las1=np.column_stack((las1, LGR1, LAR_YOUNG1, LAR_OLD1, STEIBER1, CLAVIER1))
    las2=np.column_stack((las2, LGR2, LAR_YOUNG2, LAR_OLD2, STEIBER2, CLAVIER2))
    
    return las1, las2

def throw(las2):
    """
    making throw list which is the thickness of the .las, i.e. last depth - first depth, every 1 m.

    """

    return [i for i in range(0, int(las2[-1,0] - las2[0,0]), 1)]

def SGR_estimator(las1, las2, throw, method):
    """
    Here we calculate Shale Gauge Ratio (SGR), using Vshale from one of the five possible GR->Vshale conversion methods
    """

    #step between two nearest indexes  in meters
    step=las1[1][0] - las1[0][0] 
    if method == 'linear' or method == '':
        las_mean=(las1[:,2] + las2[:,2]) / 2 
    elif method == 'larionov_young':
        las_mean=(las1[:,3] + las2[:,3]) / 2
    elif method == 'larionov_old':
        las_mean=(las1[:,4] + las2[:,4]) / 2
    elif method == 'steiber':
        las_mean=(las1[:,5] + las2[:,5]) / 2
    elif method == 'clavier':
        las_mean=(las1[:,6] + las2[:,6]) / 2
    SGRlist=[]
    for i in throw:
        X=int(round(i / step))
        depth=las2[X:,0].tolist()
        inner=[]
        j = 0
        while j < len(depth):
            inner.append(((((np.mean(las_mean[j:X + j])) * (las2[X + j, 0] - las2[j,0])) / i) * 100))
            j+=1
        SGRlist.append(inner)
    return SGRlist

def Zmax_estimator(las2, throw):
    step=las2[1][0] - las2[0][0]
    Zmaxlist=[]
    for i in throw:
        X=int(round(i / step))
        depth=las2[X:,0].tolist()
        inner=[]
        j=0
        while j < len(depth):
            inner.append((depth[j] + i))
            j+=1
        Zmaxlist.append(inner)
    return Zmaxlist

def BP_critical_yielding(SGR, Zmax):
    """ 
    Yielding Buoyancy Pressure (BP) in psi

    """
    BP_Y=[]
    for s, z in zip(SGR, Zmax):
        inner=[]
        for i in range(0, len(s)):
            if z[i] < 3000:
                inner.append((s[i] * 0.175 - 3.5) * 14.5038)
            elif z[i] > 3500:
                inner.append((s[i] * 0.15 + 1.9) * 14.5038)
            else:
                inner.append((s[i] * 0.17 + 0.92) * 14.5038)
        BP_Y.append(inner)
    return BP_Y

def BP_critical_bretan(SGR, Zmax):
    """ 
    Bretan Buoyancy Pressure (BP) in psi

    """
    BP_B=[]
    for s, z in zip(SGR, Zmax):
        inner=[]
        for i in range(0, len(s)):
            if z[i] < 3000:
                inner.append((10 ** (s[i] / 27 - 0.5)) * 14.5038)
            elif z[i] > 3500:
                inner.append((10 ** (s[i] / 27)) * 14.5038)
            else:
                inner.append((10 ** (s[i] / 27 - 0.25)) * 14.5038)
        BP_B.append(inner)
    return BP_B

def Kf_ma(SGR, Zmax, Zf_0):
    """
    Sperrevik mercucy-air permeability from SGR in mD.
    Zf by default is 100m
    """

    #Sperrevik's Constants
    A1=80000
    A2=19.4
    A3=0.00403
    A4=0.0055
    A5=12.5

    Kf_ma=[]
    for s, z in zip(SGR, Zmax):
        inner=[]
        for i in range(0, len(s)):
            inner.append(A1 * math.exp((-(A2 * (s[i] / 100) + A3 * z[i] + ((z[i] - z[0] + Zf_0) * A4 - A5) * ((1 - s[i] / 100) ** 7)))))
        Kf_ma.append(inner)
    return Kf_ma


def Pf_ma(Kf_ma):
    """
    mercury-Air fault threshold pressure from Kf in psi
    """
    Pf_ma=[]
    for throw in Kf_ma:
        inner=[]
        for i in range(len(throw)):
            inner.append(31.838 * (throw[i] ** -0.3848))
        Pf_ma.append(inner)
    return Pf_ma

def Pf_hw(Pf_ma, Y_HW, Y_MA, THETA_HW, THETA_MA):
    """
    Sperrevik HC-water capillary threshold pressure (psi)
    default parameters are:
    Y_HW = 30 - fluid tension for HC, dynes/cm
    Y_MA = 480 - fluid tension for mercury and air, dynes/cm
    THETA_HW = 30 - contact angle for HC, degrees
    THETA_MA = 40 - contact angle for for mercury and air, degrees
    """
    Pf_hw=[]
    for throw in Pf_ma:
        inner=[]
        for i in range(len(throw)):
            inner.append((Y_HW * math.cos(np.deg2rad(THETA_HW)) * throw[i]) / (Y_MA * math.cos(np.deg2rad(THETA_MA))))
        Pf_hw.append(inner)
    return Pf_hw

def column_height(Pf_hw, BP_Y, BP_B, DEN_WATER, DEN_HW):

    """
    HCCH estimation usin three methods

    """
    ## HCCH with Sperrevik, m
    HCCH_S=[]
    for throw in Pf_hw:
        inner=[]
        for i in range(len(throw)):
            inner.append((throw[i] / (0.433 * (DEN_WATER - DEN_HW))) * 0.3048)
        HCCH_S .append(inner)
        
    ## HCCH with Yielding, m
    HCCH_Y=[]
    for throw in BP_Y:
        inner=[]
        for i in range(len(throw)):
            inner.append((throw[i] / (0.433 * (DEN_WATER - DEN_HW))) * 0.3048)
        HCCH_Y .append(inner)

    ## HCCH with Bretan, m    
    HCCH_B=[]
    for throw in BP_B:
        inner=[]
        for i in range(len(throw)):
            inner.append((throw[i] / (0.433 * (DEN_WATER - DEN_HW))) * 0.3048)
        HCCH_B.append(inner)
        
    return HCCH_S, HCCH_Y, HCCH_B

def convert_to_array(list_of_lists):
    """ 
    we need this function to fill our triangle with NONE values
    """
    length=max(map(len, list_of_lists))
    return  np.array([[None] * (length - len(throw)) + throw for throw in list_of_lists])

def dataframe(throw, las, SGR, Kf_ma, HCCH_S, HCCH_Y, HCCH_B):
    """
    At last we are making dataframe with coordinate X - throw, Y - depth and 5 properties
    """
    #x coordinate - throw
    x_throw=[]
    i = 0
    while i < len(las[:,0]):
        x_throw.append(throw)
        i += 1
        
    #y coordinate - depth (the longest las)
    y_depth=[]
    for i in las[:,0]:
        inner=[]
        j=0
        while j < len(throw):
            inner.append(i)
            j+=1
        y_depth.append(inner)
        
    #combining in pandas dataframe
    return pd.DataFrame(data = {'throw, m': np.asarray(x_throw).flatten(), 'depth, m': np.asarray(y_depth).flatten(),
                               'SGR, %': convert_to_array(SGR).T.flatten(), 
                                'Kfma, mD': convert_to_array(Kf_ma).T.flatten(),
                               'HCCH Sperrevik(2002), m': convert_to_array(HCCH_S).T.flatten(),
                               'HCCH Yielding(2012), m': convert_to_array(HCCH_Y).T.flatten(),
                               'HCCH Bretan(2003), m': convert_to_array(HCCH_B).T.flatten()})