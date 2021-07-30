# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:06:36 2020

@author: Wilhelm Hodder
"""

import numpy as np
import iris
from iris.analysis import trajectory
import iris.quickplot as qplt
import iris.plot as iplt
import iris.coord_categorisation
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import mlab,colors,cm
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as PathEffects
from scipy import interpolate,stats
#import scipy
import os
import sys
from copy import deepcopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import sin, cos, sqrt, atan2, log10
import datetime




