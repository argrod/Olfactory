import scipy as sp
import numpy as np
import pandas as pd
import utm
from statistics import mode
import torch
import os, re, glob, pyproj, math, datetime
from sys import platform
from scipy.optimize import minimize
import matplotlib.pyplot as plt

timeWindow = 51
cutlength = 45
cutv = 4.1667
constv = 34.7/3.6

def Likelihoodww(data1: np.ndarray,data2: np.ndarray,cv: np.ndarray): # calculate log-likelihood of the model
    def f(par):
        b = cv/sp.special.gamma(1+1/par[0])
        L = 0
        for i,j in zip(data1,data2):
            r1 = np.sqrt(pow((i*np.cos(j) - par[3]),2) + pow((i*np.sin(j) - par[4]),2))
            rx = (i*np.cos(j)-par[3])/r1
            ry = (i*np.sin(j)-par[4])/r1
            lp = (par[0]-2) * math.log(r1) - (r1/b)**par[0] + par[1]*rx + par[2]*ry + math.log(par[0]) - math.log(b) + (1-par[0])*math.log(b) - math.log(sp.special.iv(0,np.sqrt(pow(par[1],2) + pow(par[2],2)),))
            L = L+lp
        return -1/L
    return f
for i,j in zip(df.track_speed[windows[x]],df.track_direction[windows[x]]):
    r1 = np.sqrt(pow((i*np.cos(j) - tester[3]),2) + pow((i*np.sin(j) - tester[4]),2))
    rx = (i*np.cos(j)-tester[3])/r1
    ry = (i*np.sin(j)-tester[4])/r1
    lp = (tester[0]-2) * math.log(r1) - (r1/constv/sp.special.gamma(1+1/tester[0]))**tester[0] + tester[1]*rx + tester[2]*ry + math.log(tester[0]) - math.log(constv/sp.special.gamma(1+1/tester[0])) + (1-tester[0])*math.log(constv/sp.special.gamma(1+1/tester[0])) - math.log(sp.special.iv(0,np.sqrt(pow(tester[1],2) + pow(tester[2],2)),))


def Weibull_sd(a,b): # standard deviation of Weibull distribution
    return b*np.sqrt(sp.special.gamma(1+2/a) - sp.special.gamma(1+1/a)*sp.special.gamma(1+1/a))

def Weibull_mean(a,b): # mean of Weibull distribution
    return b*sp.special.gamma(1+1/a)

def Von_Mises_sd(kappa): # standard deviation of von Mises distribution
    return 1/np.sqrt(kappa)

def readAxyGPS(filename): # read in AxyTrek GPS data (txt files)
    df = pd.read_csv(filename, sep = "\t", header = None, usecols = [0,1,2,3],
    names = ['Date','Time','lat','lon'])
    df['DT'] = pd.to_datetime(df['Date'] + " " + df['Time'],format="%d/%m/%Y %H:%M:%S")
    return df

def readBIPAxy(filename): # read in AxyTrek data as formatted by BIP system
    df = pd.read_csv(filename, sep = ",", header = 0, usecols = [0,1,2], names = ['DT','lat','lon']).dropna().reset_index()
    df['DT'] = pd.to_datetime(df['DT'].str[0:-6], format = "%Y-%m-%d %H:%M:%S")
    return df

def nearest(items, pivot): # find the nearest time position
    return min(items, key=lambda x: abs(x - pivot))

def timeRescale(dat,tdiff): # calculated indeces for rescaling time (DT) for regular sampling every tdiff mins
    return dat.iloc[np.arange(0,len(dat),step=np.timedelta64(tdiff,'m')/np.timedelta64(mode(np.diff(dat['DT'])),'s')).astype(int),:].reset_index()
np.timedelta64(np.diff(df.DT),'m')
def distLatLon(lat,lon):
    X,Y,_,_ = utm.from_latlon(np.array(lat),np.array(lon))
    return np.append(np.nan,np.sqrt((np.diff(X))**2 + (np.diff(Y))**2))

def spTrav(DT,X,Y,threshold=0): # distance and speed from time (DT), lat, and lon with option for a threshold speed
    dist = distLatLon(X,Y)
    speed = (dist)/np.array(np.append(np.nan,np.diff(DT)/np.timedelta64(1,'s')))
    if threshold != 0:
        while np.nanmax(speed) > threshold:
            lat = lat[speed < threshold]
            lon = lon[speed < threshold]
            dist = distLatLon(X,Y)
            speed = (dist)/np.array([np.nan,np.diff(DT)/np.timedelta64(1,'s')])
    return dist, speed

def prePare(filename, convertToMin: bool = True): # prepare BIP data as per required for Goto original method. Adds columns 'dt' (elapsed time from fix to previous time point in seconds), 'dist' (distance travelled from previous point in m), 'track_speed' (in m/sec), 'track_direction' in rad
    df = readBIPAxy(filename)
    if convertToMin:
        df = timeRescale(df,1)
    df['dt'] = np.append(np.nan,(np.diff(df['DT']) / np.timedelta64(1,'s')).astype(float))
    X,Y,_,_ = utm.from_latlon(np.array(df['lat']),np.array(df['lon']))
    vg_x_obs = np.diff(X)
    vg_y_obs = np.diff(Y)
    df['dist'], df['track_speed'] = spTrav(df['DT'],df['lat'],df['lon'])
    df['track_direction'] = np.append(np.nan,[math.atan2(vg_y_obs[x],vg_x_obs[x]) for x in range(len(vg_x_obs))])
    return df

def A1inv(x):
    if ((0 <= x) * (x < 0.53)):
        return 2 * x + pow(x,3) + (5 * pow(x,5))/6
    else:
        if (x < 0.85):
            return -0.4 + 1.39 * x + 0.43/(1-x)
        else:
            return 1/(pow(x,3) - 4 * pow(x,2) + 3 * x)

def windowFit(DF,tp,cutT,cutV,cutLength,windowSize):
    nxtPoint = min(range(len(DF.DT)), key=lambda i: abs(DF.DT[i]-DF.DT[tp]))
    if (nxtPoint - tp) > cutLength:
        if sum((DF.track_speed[start:end] > cutV) & (DF.dt[start:end] < cutT) & (DF.track_direction != 100)) >= cutLength:
            return np.array(range(tp,nxtPoint))[[i for i, x in enumerate((DF.track_speed[tp:nxtPoint] > cutV) & (DF.dt[tp:nxtPoint] > cutT) & (DF.track_direction[tp:nxtPoint] != 100)) if x]]

# generate a function to create windows capable of running the estimation method. Requirements are 51 minutes of data, with 75% of expected samples
def windowFit(DF,start,end,end2,cutT,cutV,cutLength):
    if (DF.DT[end] - DF.DT[start]) > np.timedelta64(cutT,'m'):
        if sum((DF.track_speed[start:end] > cutV) & (DF.dt[start:end] < cutT) & (DF.track_direction[start:end] != 100)) >= cutLength:
            if (DF.DT[end2 - 1] - DF.DT[start]) < np.timedelta64(cutT,'m'):
                while (DF.DT[end2] - DF.DT[start]) < np.timedelta64(cutT,'m'):
                    end2 = end2 + 1
            return range(start,end2)
        else:
            return range(start,end)

def findWindows(DF,windowlength=51): # calculate appropriate windows from datetimes using minimum 
    # start from minimum possible point
    fs = 1/np.abs(np.timedelta64(mode(np.diff(DF.DT)),'m')).astype(int) # in Hz
    expSamp = round(51 * fs) # expected number of samples
    cutlength = round(45/51 * fs * 51)
    error_of_sampling_interval = 5 * fs # give 5 seconds of leeway for samples to be found in
    cutt = ((60 * fs) + error_of_sampling_interval).astype(int)
    return list(filter(None,[windowFit(DF,x,x+expSamp,x+cutlength,cutt,cutv,cutlength) for x in range(len(DF) - expSamp)]))

def initPars(head,spd,hed,cv):
    inita = 0
    while (inita < 5):
        inita = np.abs(np.random.normal(12.5,5))
    meangd = np.arctan2(np.mean(np.sin(hed)),np.mean(np.cos(hed)))
    inithd = meangd + (head/3 * np.pi/2)
    initkappa = A1inv(np.mean(np.cos(hed - meangd)))
    initmux = initkappa * np.cos(inithd)
    initmuy = initkappa * np.sin(inithd)
    initwx = np.mean(spd) * np.cos(meangd) - cv * np.cos(inithd)
    initwy = np.mean(spd) * np.sin(meangd) - cv * np.sin(inithd)
    return meangd,initkappa,[inita,initmux,initmuy,initwx,initwy]

def windOptims(spd,hed,cv,pars):
    return minimize(Likelihoodww(spd,hed,cv),pars,method="L-BFGS-B")

def yokoTate(optimAns,cv):
    yoko = Von_Mises_sd(np.sqrt(optimAns[2]*optimAns[2] + optimAns[1]*optimAns[1])) * Weibull_mean(optimAns[0],cv/sp.special.gamma(1+1/optimAns[0]))
    tate = Weibull_sd(optimAns[0], cv/sp.special.gamma(1 + 1/optimAns[0]))
    return yoko, tate

def ensureOptimConv(optimAns,spd,hed,cv):
    newAns = optimAns
    for try_no in range(100):
        pars = [newAns.x[0],newAns.x[1],newAns.x[2],newAns.x[3],newAns.x[4]]
        newAns = windOptims(spd,hed,cv,pars)
        if newAns.success:
            break
    return newAns, yokoTate(newAns.x,cv)

def windOptim(initHead,spd,head,cv):
    # define initial parameters
    try:
        _,_,par = initPars(initHead,spd,head,cv)
    except:
        _,_,par = initPars(initHead,spd,head,cv)
    # run first optimisation
    answ = windOptims(spd,head,cv,par)
    # SD of heading vector perpendicular (yoko) and along (tate) mean direction
    yoko,tate=yokoTate(answ.x,cv)
    # repeat MLE to ensure convergence
    if not np.isnan(tate):
        answ,[yoko,tate] = ensureOptimConv(answ,spd,head,cv)
        return answ,yoko,tate
    else:
        return np.nan,np.nan,np.nan

def headSpdDir(spd,hed,answ):
    nvx = spd * np.cos(hed) - answ.x[3]
    nvy = spd * np.sin(hed) - answ.x[4]
    nr = np.sqrt(nvx**2 + nvy**2)
    nd = np.arctan2(nvy,nvx)
    return nr,nd

def GOFtests(hed,nr,nd,answ):
    # mean track direction
    meangd = np.arctan2(np.mean(np.sin(hed)),np.mean(np.cos(hed)))
    nrp = sp.stats.kstest(nr,'weibull_min',args=[answ.x[0],constv/sp.special.gamma(1+1/answ.x[0])])
    # direction of heading vector
    mu = math.atan2(answ.x[2],answ.x[1])
    kappa = np.sqrt(pow(answ.x[1],2) + answ.x[2]*answ.x[2])
    ndp = sp.stats.kstest(nd,'vonmises',args=[kappa,mu])
    # correlation test between direction and speed of heading vector
    cnrnd = sp.stats.pearsonr(nr,nd)
    # describe conditions required
    cond1 = (yoko/tate) > 1
    cond2 = (np.cos(meangd)*np.cos(mu) + np.sin(meangd)*np.sin(mu)) > 0
    cond3 = (nrp.pvalue > 0.05) * (ndp.pvalue > 0.05) * (cnrnd.pvalue > 0.05)
    return [cond1,cond2,cond3]

outfile = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/PyOut.csv"
with open(outfile,'a') as fd:
    fd.write(['Time','Lat','Lon','X','Y','WindHead','BirdHead'])



# Process
# read in sample from BIP system
if platform == "darwin":
    filename = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/SampleAxyTrek.csv"
else:
    filename = "I:/My Drive/PD/Data/TestingData/SampleAxyTrek.csv"
df = prePare(filename)

cutv = 4.1667
constv = 34.7/3.6
# first generate windows over which the estimation method will be run
windows = findWindows(df,51)
# set initial value for max likelihood and best answer
max_like = np.nan
answ_best = np.nan
# then perform MLE for initial headings ranging from -3 to 3
hd_try = 3
for x in range(len(windows)):
    for id_hd in range(-hd_try,hd_try):
        # perform optimisation routine
        answ,yoko,tate = windOptim(id_hd,df.track_speed[windows[x]],df.track_direction[windows[x]],constv)
        if (not np.isnan(tate)) & (answ.status==0):
            # calculate speeds and headings
            nr,nd = headSpdDir(df.track_speed[windows[x]],df.track_direction[windows[x]],answ)
            # GOF tests
            if np.isnan(max_like) & (np.prod(GOFtests(df.track_direction[windows[x]],nr,nd,answ)) == 1):
                max_like = answ.fun
                answ_best = answ
                print(x)
            if (answ.fun > max_like) & (np.prod(GOFtests(df.track_direction[windows[x]],nr,nd,answ)) == 1):
                max_like = answ.fun
                answ_best = answ
            else:
                answ_best = np.nan

#[range(750, 801), range(751, 802), range(752, 803), range(753, 804), range(754, 805), range(755, 806), range(756, 807), range(757, 808), range(758, 809), range(759, 810)]


nr,nd = headSpdDir(df.track_speed[windows[0]],df.track_direction[windows[0]],testAns)
GOFtests(df.track_direction[windows[0]],nr,nd,testAns)

_,_,testPars = initPars(-3,df.track_speed[windows[0]],df.track_direction[windows[0]],constv)
testAns = windOptims(df.track_speed[windows[0]],df.track_direction[windows[0]],constv,testPars)
testYoko,testTate = yokoTate(testAns.x,constv)



testmgd,testkappa,testpars = initPars(-3,df.track_speed[windows[0]],df.track_direction[windows[0]],constv)
testAns=windOptims(df.track_speed[windows[0]],df.track_direction[windows[0]],constv,testpars)
ensureOptimConv(testAns,df.track_speed[windows[0]],df.track_direction[windows[0]],constv)
testAns.x
Von_Mises_sd(np.sqrt(testAns.x[2]*testAns.x[2] + testAns.x[1]*testAns.x[1])) * Weibull_mean(testAns.x[0],constv/sp.special.gamma(1+1/testAns.x[0]))
yokoTate(testAns.x,constv)
testAns,[yoko,tate]=ensureOptimConv(testAns,df.track_speed[windows[0]],df.track_direction[windows[0]],constv)

testAns[2]*testAns[2] + testAns[1]**testAns[1]

# # set initial value for max likelihood and best answer
# max_like = np.nan
# answ_best = np.nan
# # loop through each hd_try value
# # for heading in range(-3,3): 
# heading = -3   
# # create initial parameters to maximise likelihood function
# initPars = initHead(heading,rrow[windows[0]],drow[windows[0]],constv)
# # run first optimisation
# answ = windOptims(rrow[windows[0]],drow[windows[0]],constv,initPars)
# # generate sd of heading vector perp to mean direction and along mean direction
# yoko,tate = yokoTate(answ,constv)
# # repeat MLE to ensure optimisation convergence
# if not np.isnan(tate):
#     answ = ensureOptimConv(answ,rrow[windows[0]],drow[windows[0]],constv)
#     yoko,tate = yokoTate(answ,constv)
# # if successfully converged
# if answ.success:
#     nr, nd = headSpDir(rrow[windows[0]],drow[windows[0]],answ.x[3],answ.x[4])
#     # define mean heading
    # meangd = np.arctan2(np.mean(np.sin(drow[windows[0]])),np.mean(np.cos(drow[windows[0]])))
    # # test goodness of fit
    # if np.isnan(max_like) & (corTests(nr,nd,answ,constv,meangd) == 1):
    #     max_like = answ.fun
    #     answ_best = answ
    # else:
    #     answ_best = np.nan
    # if (answ.fun > max_like) & (corTests(nr,nd,answ,constv,meangd) == 1):
    #     max_like = answ.fun
    #     answ_best = answ
    # else:
    #     answ_best = np.nan
# *** goodness of fit tests ***
# speed of heading vector
nrp = sp.stats.kstest(nr,'weibull_min',args=[answ.x[0],constv/sp.special.gamma(1+1/answ.x[0])])
# direction of heading vector
mu = math.atan2(answ.x[2],answ.x[1])
kappa = np.sqrt(pow(answ.x[1],2) + answ.x[2]*answ.x[2])
ndp = sp.stats.kstest(nd,'vonmises',args=[kappa,mu])
# correlation test between direction and speed of heading vector
cnrnd = sp.stats.pearsonr(nr,nd)
# describe conditions required
cond1 = (yoko/tate) > 1
cond2 = (np.cos(meangd)*np.cos(mu) + np.sin(meangd)*np.sin(mu)) > 0
cond3 = (nrp.pvalue > 0.05) * (ndp.pvalue > 0.05) * (cnrnd.pvalue > 0.05)
[cond1,cond2,cond3]