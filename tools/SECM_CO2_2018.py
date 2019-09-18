import numpy as np
from matplotlib import pyplot as plt
from jdutil import *
from netCDF4 import date2num, num2date

# This is an implementation of Max Reuter's IDL code
'''
def ufunc1_2018(lats, dates, a):
    # LATS are in radians!

    jdates = np.array([date_to_jd(d.year, d.month, d.day)
                       for d in dates])
    t = (jdates - date_to_jd(2003, 1, 1)) / 365.25
    dumt1 = np.tanh(a[3] * lats + a[4])
    dumt2 = np.tanh(a[5] * lats + a[6]) + a[7]
    dumt3 = np.tanh(a[8] * lats + a[9])

    season = (
        a[11] * dumt2 * np.sin(t * 2.0 * np.pi + a[10] * dumt3) +
        a[12] * dumt2 * np.cos(t * 2.0 * np.pi + a[10] * dumt3) +
        a[13] * dumt2 * np.sin(t * 4.0 * np.pi + a[10] * dumt3) +
        a[14] * dumt2 * np.cos(t * 4.0 * np.pi + a[10] * dumt3)
        )

    f = a[0] + a[1] * t + a[2] * dumt1 + season
    return f, season

def secm_co2_2018(lats, dates):

    c = np.array([38.779246,   2.0831194])
'''

def ufunc1_2018(lats, dates, a):

    dumt2 = np.tanh(a[5] * lats + a[6]) + a[7]
    dumt3 = np.tanh(a[8] * lats + a[9])

    season = (
        a[11] * dumt2 * np.sin(t * 2.0 * np.pi + a[10] * dumt3) +
        a[12] * dumt2 * np.cos(t * 2.0 * np.pi + a[10] * dumt3) +
        a[13] * dumt2 * np.sin(t * 4.0 * np.pi + a[10] * dumt3) +
        a[14] * dumt2 * np.cos(t * 4.0 * np.pi + a[10] * dumt3)
        )

    return season


def ufunc2_2018(lats, dates, a):

    season = ufunc1_2018(lats, dates, a)

    dumt1 = np.tanh(a[3] * lats + a[4])
    return a[0] + a[1]*t + a[2] * dumt1 + season


def ufunc3_2018(lats, dates, a, c, p, pt=0.2):

    result = np.zeros((len(dates), len(p)))

    result += c[0] * (0.5 * pt * pt - pt)
    result += c[1] * (pt - 0.5 - 0.5 * pt * pt)

    print(result.shape)
    result = (result.T * ufunc1_2018(lats, dates, a)).T

    p_low = p <= pt
    p_high = p > pt

    result[:, p_low] += c[0] * p[p_low]
    result[:, p_high] += c[0] * pt + np.outer(ufunc1_2018(lats, dates, a),
                                              c[1] * (p[p_high] - pt))

    return (result.T + ufunc2_2018(lats, dates, a)).T
