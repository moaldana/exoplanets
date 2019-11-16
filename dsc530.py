"""This file contains code for use with DSC530 exercises 
by Mario Aldana.

"""

import thinkplot
import thinkstats2

import survival

import pyvo
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

 
    
def SpearmanCorr(xs, ys):
    xranks = pd.Series(xs).rank()
    yranks = pd.Series(ys).rank()
    return Corr(xranks, yranks)


def Corr(xs, ys):
    xs = np.asarray(xs)
    ys = np.asarray(ys)

    meanx, varx = thinkstats2.MeanVar(xs)
    meany, vary = thinkstats2.MeanVar(ys)

    corr = Cov(xs, ys, meanx, meany) / np.sqrt(varx * vary)
    return corr


def Cov(xs, ys, meanx=None, meany=None):
    xs = np.asarray(xs)
    ys = np.asarray(ys)

    if meanx is None:
        meanx = np.mean(xs)
    if meany is None:
        meany = np.mean(ys)

    cov = np.dot(xs-meanx, ys-meany) / len(xs)
    return cov


def Jitter(values, jitter=0.5):
    n = len(values)
    return np.random.normal(0, jitter, n) + values


def ConvertDataframe(sql):
    """Convert a pyvo.dal.tap.TAPResults object by a pandas dataframe.

    sql: query of exoplanets database

    returns: result into a pandas dataframe
    """

    engine = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")
    
    # Get the records
    sqlResult = engine.search(sql)

    # Convert initial TAPResults object to astropy table,
    tblResult = sqlResult.to_table()

    # Convert astropy table object to numpy array,
    npResult = tblResult.as_array()

    # Convert the numpy array to Pandas dataframe.
    pdResult = pd.DataFrame(np.array(npResult))
    
    return pdResult



# From Week 13,

def EstimateSurvival(resp):
    """Estimates the survival curve.

    resp: DataFrame of respondents

    returns: pair of HazardFunction, SurvivalFunction
    """
    complete = resp[resp.notdivorced == 0].duration.dropna()
    ongoing = resp[resp.notdivorced == 1].durationsofar.dropna()

    hf = survival.EstimateHazardFunction(complete, ongoing)
    sf = hf.MakeSurvival()

    return hf, sf


def EstimateSurvivalByDecade(groups, **options):
    """Groups respondents by decade and plots survival curves.

    groups: GroupBy object
    """
    thinkplot.PrePlot(len(groups))
    for name, group in groups:
        _, sf = EstimateSurvival(group)
        thinkplot.Plot(sf, **options)


def main():
    pass
    

if __name__ == '__main__':
    main()
