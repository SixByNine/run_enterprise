#!/usr/bin/env python
'''Read pulsar glitches solution files and plot the correlation between parameters.
    Written by Y.Liu (yang.liu-50@postgrad.manchester.ac.uk).'''

#import argparse
import glob
#import sys, os
#import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import numpy as np
#import numpy.ma as ma
import pandas as pd
#import pulsar_glitch as pg
#from __future__ import print_function
from uncertainties import ufloat
#from scipy.optimize import curve_fit
#from scipy import linalg
#from astropy import coordinates as coord
#from astropy import units as u
#from matplotlib.ticker import FormatStrFormatter
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#from mpl_toolkits.axes_grid.inset_locator import inset_axes


    
file_name_list = glob.glob('slt_*.csv')
file_name_list.sort()
pulsar_name = [file_name.split('_')[-1].split('.')[0] for file_name in file_name_list]

summary = 'summary.csv'
sum_table = pd.DataFrame()

for i, solution in enumerate(file_name_list):
    table = pd.read_csv(solution, index_col=[0,1])
    sum_table = pd.concat([sum_table, table])

sum_table.to_csv(summary)
print('Summary solution', sum_table)
sum_table.plot(x='F0', y='GLF0(instant)')
