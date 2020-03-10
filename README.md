# Fault Triangle
A Python script for subsurface geoscientists. The script estimates Shale Gauge Ratio (SGR) and hydrocarbon column height (HCCH) of the subsurface fault between two wells using .las files with Gamma Ray. The script plots SGR and HCCH as a triangle diagram, where X-coordinate is a fault throw (vertical displacement along the fault), Y-coordinate is a depth and color - SGR or HCCH.

## Short technical background

This tool may seem too specific for everyone but geoscientists. Here I'll try to explain technical background.

It is very important to understand if subsurface fault holds or transmits fluids. Why? Imagine CO2 sequestration project. In order to get rid of CO2 from burning fuels we can inject this gas inside the Earth. But if there are conductive faults, CO2 can leak out. Similar with oil or gas. To predict hydrocarbon accumulation in faulted area it is always neccessary to understand if the fault leaks or seals. That knowledge can be significant for risk reduction during petroleum exploration, where one well can cost up to hundreds of millions of $.


### How to look inside the Earth?

Unfortunately, we cannot have a look inside, but we can estimate the geophysical properties of the rock by running some measurements inside well bores. This is called logging. These measurements are proxies to real rocks. One of the paramount logging methods is a Gamma ray logging - "a method of measuring naturally occurring gamma radiation to characterize the rock or sediment in a borehole or drill hole." Let's omit all the complicated Earth physics, and here I can tell you that geoscientists use Gamma Ray (GR) for estimation of the fraction of shale in the rock. The more GR value, the more shale content in the rock. As you may guess, the more shale the less permeable our rock is. There are five well-known equations of conversion of GR to Vshale (amount of shale in the rock volume): linear, Larionov_young, Larionov_old, Steiber and Clavier. I'm not going to bother you with formulas. We need to know just one thing: every conversion method gives slightly different results. It's up to user which method to use.


### Estimation of Fault Seal Capacity

 Let's talk about faults now. If we assume that we have ONLY specific rocks (shale, sand and silt) under subsurface, then we can estimate SGR along the fault plane - "At any point on a fault surface, the SGR is equal to the net shale/clay content of the rocks that have slipped past that point." I.e. if we know how much shale in our rocks and what the vertical fault displacement (or *throw*) is, then we can calculate this SGR property. If you are still reading, you may notice, that SGR is an amount of shale at every point of the fault plane.
 Now let's discuss how we can estimate sealing properties of the rock from SGR. The easiest way is to estimate how much fluid fault rock can hold with given SGR. If it can hold 100 m, then, let's say it's impermeable for that fluid. It it can hold < 1m of the fluid, then the fault conducts this fluid and the fluid can leak away. Before going to column height estimation, we have to estimate a maximum fluid pressure with the given SGR. There are three different methods of pressure estimation:

 * Sperrevik et al 2002. In order to estimate the column of the fluid we need to know: SGR, maximum depth of the fault and an initial depth of the fault.
 * For two other methods Yielding 2012, and Bretan 2003, we need to know: SGR and maximum depth of the fault.

After that, we convert our maximum pressures to column heights of any fluid.

It is important to know, that if we want to calculate everything above, we need to know some basic parameters for every fluid: 1. Density of the fluids (usually water and petroleum), 2.Contact angle between two fluids and 3. Fluid tension. 

Hence, to summarize, using this tool we: convert GR to Vshale using one of the five methods, calculate SGR for every throw step (vertical displacement of the fault), calculate critical pressures using three methods, calculate column heights, make a useful plot where we can see on which depth our fault is leaking.

## Code design

My code design is straightforward. I split it in three parts: data.py takes care of data processing and calculations; plot.py takes care of plotting and main.py aggregates everything in one working script. 

Code is written on Python 3.8 

In my code I used NumPy, Pandas and Matplotlib as main working libraries, as well library called *lasio* for taking care of .las files from the wells. All the nessessary libraries are listed in requirements.txt

## Installation instructions

1. Download data.py, plot.py and main.py files, as well as .las files with the required information from the wells.
2. Install requred packages on your virtual environment using requirements.txt
3. Run main.py and follow instructions

## Known Bugs
There are several warning messages, which I wanted to remove but didn't have time.

