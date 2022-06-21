# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 13:45:56 2021

@author: kse18nru
"""

import numpy as np
import pandas as pd
import iris


def potential_temperature(T,p, p0 = 100000, R = 287, cp = 1004):
    return T*(p0/p)**(R/cp)

def specific_humidity(w):
    return w/(1+w)