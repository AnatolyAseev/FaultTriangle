#!/usr/bin/env python3
""" 
plot.py is a plotting tool.

we are plotting 1)SGR triangle, 2)Kfma, 3)HCCH Sperrevik 2002, 4)HCCH Yielding 2012, 5)HCCH Bretan 2003

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


def plotHCCH(df, output_path):
    """ 
    plotting three HCCH plots
    """
    property= list(df)[4:]
    x=df['throw, m']
    y=df['depth, m']
    for i in property:
        c=df[i]
        fig, ax=plt.subplots(figsize=(20, 20))
        ax.scatter(x, y, c=c, cmap=plt.get_cmap('jet', 20), norm=matplotlib.colors.LogNorm(vmin=10, vmax=1000))
        ax.legend()
        ax.grid(True)
        cax=ax.scatter(x, y, c=c, cmap=plt.get_cmap('jet', 20), norm=matplotlib.colors.LogNorm(vmin=10, vmax=1000))
        cbar=fig.colorbar(cax)
        cbar.ax.set_ylabel(str(i), rotation=90)
        ax.set_xlabel('Throw, m')
        ax.set_ylabel('Depth, m')
        ax.set_xlim(x.iloc[0], x.iloc[-1])
        ax.set_ylim(y.iloc[0], y.iloc[-1])
        ax.invert_yaxis()
        plt.rcParams.update({'font.size': 20})
        plt.savefig(output_path + i + '.png')
        plt.close('all')

def plotKf(df, output_path):
    x=df['throw, m']
    y=df['depth, m']
    c=df['Kfma, mD']
    fig, ax=plt.subplots(figsize=(20, 20))
    ax.scatter(x, y, c=c, cmap=plt.get_cmap('jet', 10), norm=matplotlib.colors.LogNorm())
    ax.legend()
    ax.grid(True)
    cax=ax.scatter(x, y, c=c, cmap=plt.get_cmap('jet', 10), norm=matplotlib.colors.LogNorm())
    cbar=fig.colorbar(cax)
    cbar.ax.set_ylabel('Kfma, mD', rotation=90)
    ax.set_title('Kf mercury-air using Sperrevik et al (2002)')
    ax.set_xlabel('Throw, m')
    ax.set_ylabel('Depth, m')
    ax.set_xlim(x.iloc[0], x.iloc[-1])
    ax.set_ylim(y.iloc[0], y.iloc[-1])
    ax.invert_yaxis()
    plt.savefig(output_path + 'Kf.png')
    plt.close('all')

def plotSGR(df, output_path):
    x=df['throw, m']
    y=df['depth, m']
    c=df['SGR, %']
    fig, ax=plt.subplots(figsize=(20, 20))
    ax.scatter(x, y, c=c, cmap=plt.get_cmap('jet', 50), norm=matplotlib.colors.DivergingNorm(vcenter=18, vmin=0, vmax=100))
    ax.legend()
    ax.grid(True)
    cax=ax.scatter(x, y, c=c, cmap=plt.get_cmap('jet', 50), norm=matplotlib.colors.DivergingNorm(vcenter=18,vmin=0, vmax=100))
    cbar=fig.colorbar(cax)
    cbar.ax.set_ylabel('SGR,%', rotation=90)
    ax.set_title('SGR simulations, %')
    ax.set_xlabel('Throw, m')
    ax.set_ylabel('Depth, m')
    ax.set_xlim(x.iloc[0], x.iloc[-1])
    ax.set_ylim(y.iloc[0], y.iloc[-1])
    ax.invert_yaxis()
    plt.savefig(output_path + 'SGR.png')
    plt.close('all')