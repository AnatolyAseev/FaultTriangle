#!/usr/bin/env python3

"""
we are running file and compiling everything
"""

import numpy as np
from numpy import arange
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import random
import lasio as las
import math
import itertools

import data
import plot

def get_with_default(type_, message, default):
    user_input = input(message)
    try:
        r=type_(user_input)
    except ValueError:
        r=default
    return r

def main():
    input("### Press *Enter* to load .las files with Depth and Gamma Ray tracks ### ...   ")
    path1=input("### Please copy-paste ot type the path to the *fisrt* .las file ###   ")
    path2=input("### Please copy-paste ot type the path to the *second* .las file ###   ")

    #load data
    las1, las2=data.loader(path1, path2)

    #extend the shorter las
    las1=data.las1_extension(las1, las2)

    #estimate Vshale
    las1, las2=data.GR_to_Vshale(las1, las2)

    #estimate SGR
    method="### Enter GR to Vshale conversion method. The options are: 'linear', 'larionov_young', 'larionov_old', 'steiber', 'clavier'. The default method is 'linear' ###.....   "
    SGR=data.SGR_estimator(las1, las2, data.throw(las2), method=get_with_default(str, method, 'linear'))

    #estimate buoyancy pressure or AFPD for Yielding and Bretan, in psi
    BP_Y=data.BP_critical_yielding(SGR, data.Zmax_estimator(las2, data.throw(las2)))
    BP_B=data.BP_critical_bretan(SGR, data.Zmax_estimator(las2, data.throw(las2)))

    #estimate mercury-air fault permeability using Sperrevik(2002) equation.
    Zf_0="### Enter integer value representing the depth of faulting or Zf (from Sperrevik, 2002) in meters for the first depth point. The default value is *100* meters ###......   "
    Kf_ma=data.Kf_ma(SGR, data.Zmax_estimator(las2, data.throw(las2)), Zf_0=get_with_default(int, Zf_0, 100))

    #estimate Pf_hw - hydrocarbon-water threshold capillary pressure. The liquid pair can be different. However, user should know input parameters 
    #of fluid tension, contact angle and density
    Y_HW="### Enter integer value representing a fluid tension for HC and water in dynes/cm. The default value for the light oil is *30* dynes/cm ###......   "
    Y_MA = "### Enter integer value representing a fluid tension for mercury and air in dynes/cm. The default value is *480* dynes/cm ###...... "
    THETA_HW="### Enter integer value representing a contact angle for HC and water, degrees. The default value for the light oil is *30* degrees ###......   "
    THETA_MA="### Enter integer value representing a contact angle for mercury and air, degrees. The default value is *40* degrees ###......   "
    Pf_hw=data.Pf_hw(data.Pf_ma(Kf_ma), Y_HW=get_with_default(int, Y_HW, 30), Y_MA=get_with_default(int, Y_MA, 480), THETA_HW=get_with_default(int, THETA_HW, 30), THETA_MA=get_with_default(int, THETA_MA, 40))

    #estimate HC column height in meters using three methods: Sperrevik(2002), Yielding(2012) and Bretan(2003). 
    DEN_WATER="### Enter float value representing a water density, g/cm^3. The default value is *1.030* g/cm^3 ###......   "
    DEN_HW="### Enter float value representing a HC density, g/cm^3. The default value is *0.700* g/cm^3 ###......   "
    HCCH_S, HCCH_Y, HCCH_B=data.column_height(Pf_hw, BP_Y, BP_B, DEN_WATER = get_with_default(float, DEN_WATER, 1.030), DEN_HW = get_with_default(float, DEN_HW, 0.700))

    # making a dataframe
    input ("Now, let's create a dataframe with X coordinate - Fault throw, m; Y coordinate - Depth, m; and five properties: SGR, Kfma, HCCH_Sperrevikm HCCH_Yielding and HCCH_Bretan... ### Press *Enter* to continue ###")
    df=data.dataframe(data.throw(las2), las2, SGR, Kf_ma, HCCH_S, HCCH_Y, HCCH_B)
    print (df.tail(10))

    #plot
    input ("Now, let's plot and save five triangle plots in any folder. X - Fault throw, m; Y - Depth, m; and five properties: SGR, Kfma, HCCH_Sperrevikm HCCH_Yielding and HCCH_Bretan... ### Press *Enter* to continue ###")
    output_path="### Please copy-paste or type a *PATH* to an output folder. Default path is a current folder *./* ###   "
    plot.plotHCCH(df, output_path=get_with_default(str, output_path, './'))
    plot.plotKf(df, output_path=get_with_default(str, output_path, './'))
    plot.plotSGR(df, output_path=get_with_default(str, output_path, './'))

    print("####################### Thank you for using the Fault Triangle script ##########################")

if __name__ == '__main__':
    main()




