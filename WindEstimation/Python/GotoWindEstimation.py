import scipy as sp
import numpy as np
import pandas as pd
import csv
import utm
from statistics import mode
import pyproj, math
from sys import platform
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from functools import partial
import cmath
import glob, re

def Likelihoodww(data1: np.ndarray,data2: np.ndarray,cv: np.ndarray): # calculate log-likelihood of the model
    def f(par):
        b = cv/sp.special.gamma(1+1/par[0])
        L = 0
        for i,j in zip(data1,data2):
            r1 = np.sqrt(pow((i*np.cos(j) - par[3]),2) + pow((i*np.sin(j) - par[4]),2))
            rx = (i*np.cos(j)-par[3])/r1
            ry = (i*np.sin(j)-par[4])/r1
            lp = (par[0]-2) * math.log(r1) - (r1/b)**par[0] + par[1]*rx + par[2]*ry + math.log(par[0]) - math.log(b) + (1-par[0])*math.log(b) - math.log(sp.special.iv(0,np.sqrt(par[1]**2 + par[2]**2),))
            L = L+lp
        return L/(-1)
    return f

# logs used on r1 (always +ve), par[0], b (cv/gamma(1+1/par[0]), always +ve), and besselI(+ve value)

def Weibull_sd(a,b): # standard deviation of Weibull distribution
    return b*np.sqrt(sp.special.gamma(1+2/a) - sp.special.gamma(1+1/a)*sp.special.gamma(1+1/a))

def Weibull_mean(a,b): # mean of Weibull distribution
    return b*sp.special.gamma(1+1/a)

def Von_Mises_sd(kappa): # standard deviation of von Mises distribution
    return 1/np.sqrt(kappa)

def readAxyGPS(filename): # read in AxyTrek GPS data (txt files)
    df = pd.read_csv(filename, sep = "\t", usecols = [0,1,2,3],
    names = ['Date','Time','lat','lon'])
    df['DT'] = pd.to_datetime(df['Date'] + " " + df['Time'],format="%d/%m/%Y %H:%M:%S")
    return df

def readBIPAxy(filename): # read in AxyTrek data as formatted by BIP system
    df = pd.read_csv(filename, sep = ",", header = 0, usecols = [0,1,2], names = ['DT','lat','lon']).dropna().reset_index()
    df['DT'] = pd.to_datetime(df['DT'].str[0:-6], format = "%Y-%m-%d %H:%M:%S")
    return df

def nearest(items, pivot): # find the nearest time position
    return min(items, key=lambda x: abs(x - pivot))

def nearestInd(items, pivot): # find the nearest index
    return min(range(len(items)), key=lambda i: abs(items[i] - pivot))

def timeRescale(dat,tdiff): # calculated indeces for rescaling time (DT) for regular sampling every tdiff mins
    return dat.iloc[np.arange(0,len(dat),step=np.timedelta64(tdiff,'m')/np.timedelta64(mode(np.diff(dat['DT'])),'s')).astype(int),:]

def angles(longitudes,latitudes):
    lon1 = longitudes.values[:-1]                                                                                       
    lat1 = latitudes.values[:-1]                                                                                        
    lon2 = longitudes.values[1:]                                                                                        
    lat2 = latitudes.values[1:]
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])   

    X = lon2 - lon1
    Y = lat2 - lat1
    return np.append(np.nan,np.arctan2(Y,X))

def gps_speed(longitudes, latitudes, timestamps):                                                                       
    """                                                                                                                 
    Calculates the instantaneous speed from the GPS positions and timestamps. The distances between the points          
    are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
    the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: pandas series of longitudes                                                                         
        latitudes:  pandas series of latitudes                                                                          
        timestamps: pandas series of timestamps                                                                         
 
    Returns:                                                                                                            
        Speed is returned an array in m/s.                                                                             
 
    Example:                                                                                                            
        >>> df['gpsSpeed'] = gps_speed(df.longitude, df.latitude, df.recordedAt)
    """                                                                                                                 
 
    assert longitudes.shape[0] > 1                                                                                      
    assert latitudes.shape[0] > 1                                                                                       
    assert timestamps.shape[0] > 1                                                                                      
 
    lon1 = longitudes.values[:-1]                                                                                       
    lat1 = latitudes.values[:-1]                                                                                        
    lon2 = longitudes.values[1:]                                                                                        
    lat2 = latitudes.values[1:]                                                                                         
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    dist = 6371 * 2 * np.arcsin(np.sqrt(a)) * 1000
    time_array = (timestamps.diff().dt.seconds ).values[1:]                                                       
 
    # Calculate the speed                                                                                               
    time_array[time_array == 0] = np.nan  # To avoid division by zero                                                   
    speed = dist / time_array                                                                                       
 
    # Make the arrays as long as the input arrays                                                                        
    speed = np.insert(speed, 0, np.nan, axis=0)
    dist = np.insert(dist, 0, np.nan, axis=0)                                                                         
 
    return dist,speed

def prePare(filename, convertToMin: bool = True, isBip: bool = True): # prepare BIP data as per required for Goto original method. Adds columns 'dt' (elapsed time from fix to previous time point in seconds), 'dist' (distance travelled from previous point in m), 'track_speed' (in m/sec), 'track_direction' in rad
    if isBip:
        df = readBIPAxy(filename)
    else:
        df = readAxyGPS(filename)
    if convertToMin:
        df = timeRescale(df,1)
    df['dt'] = np.append(np.nan,(np.diff(df['DT']) / np.timedelta64(1,'s')).astype(float))
    df['dist'], df['track_speed'] = gps_speed(df['lon'],df['lat'],df['DT'])
    df['track_direction'] = angles(df['lon'],df['lat'])
    return df.dropna().reset_index()[['DT','lon','lat','dt','dist','track_speed','track_direction']]

def A1inv(x):
    # copy of A1inv function from circular package in R
    if ((0 <= x) & (x < 0.53)):
        return 2 * x + pow(x,3) + (5 * pow(x,5))/6
    else:
        if (x < 0.85):
            return -0.4 + 1.39 * x + 0.43/(1-x)
        else:
            return 1/(pow(x,3) - 4 * pow(x,2) + 3 * x)

# def windowFit(dataF,start,cutT,cutV,cutLength,windowlength):
#     # generate windows capable of running the estimation method. Requirements are 51 minutes of data, with 75% of expected samples
#     end = start + next(filter(lambda i: i[1] >= windowlength*60, enumerate(np.cumsum(dataF.dt[0::]))))[0]
#     if sum((dataF.track_speed[start:end] > cutV) & (dataF.dt[start:end] < cutT) & (dataF.track_direction[start:end] != 100)) >= cutLength:
#         return np.array(range(start,end))[(dataF.track_speed[start:end] > cutV) & (dataF.dt[start:end] < cutT) & (dataF.track_direction[start:end] != 100)]

# def findWindows(dataF,cutv,windowlength=51): # calculate appropriate windows from datetimes using minimum 
#     # start from minimum possible point
#     fs = 60/(np.abs(np.timedelta64(mode(np.diff(pDat.DT)),'s'))).astype(int)  # in Hz
#     expSamp = round(windowlength * fs) # expected number of samples
#     cutlength = round(45/51 * fs * expSamp)
#     error_of_sampling_interval = 5 * fs # give 5 seconds of leeway for samples to be found in
#     cutt = ((60 * fs) + error_of_sampling_interval).astype(int)
#     return list(filter(lambda item: item is not None,[windowFit(dataF,x,cutt,cutv,cutlength,windowlength) for x in range(0,len(dataF)-next(filter(lambda i: i[1] >= windowlength*60, enumerate(np.cumsum(dataF.dt[::-1]))))[0]-1)]))

# starting from a center point, search reverse and forward until half the window width is available.
# as the center point *should* have n seconds of interval in it, we can then subtract fs * 1min from the window width and search for half this new value
# searchGap = (windowWidth - (fs*1)/2) where windowWidth in minutes


def windowFit(DF,start,end,cutT,cutV,cutLength,windowlength): # generate windows capable of running the estimation method. Requirements are 51 minutes of data, with 75% of expected samples
    if sum(((DF.DT[end] - DF.DT[start]) <= np.timedelta64(cutT*windowlength,'s')) & (DF.track_speed[start:end] > cutV) & (DF.dt[start:end] < cutT) & (DF.track_direction[start:end] != 100)) >= cutLength:
        return np.array(range(start,end))[(DF.track_speed[start:end] > cutV) & (DF.dt[start:end] < cutT) & (DF.track_direction[start:end] != 100)]

def rangeGen(DF,cen,wws,cutV,cutT,cutLength): # generate windows of appropriate time durations (wws seconds on either side)
    st = DF.dt[:cen][::-1].cumsum().gt(wws).idxmax() + 1
    end = DF.dt[(cen+1):].cumsum().gt(wws).idxmax() + 1
    if sum((DF.track_speed[st:end] > cutV) & (DF.dt[st:end] < cutT) & (DF.track_direction[st:end] != 100)) > cutLength:
        return range(st,end),cen

def findWindows(DF,cutv,windowlength=51): # calculate appropriate windows from datetimes using minimum 
    # start from minimum possible point
    fs = (1/np.abs(np.timedelta64(mode(np.diff(DF.DT)),'m')).astype(int)).astype(int) # in fixes per minute
    expSamp = round(51 * fs) # expected number of samples
    cutlength = round(45/51 * fs * 51)
    error_of_sampling_interval = 5 * fs # give 5 seconds of leeway for samples to be found in
    cutt = ((60 * fs) + error_of_sampling_interval).astype(int)
    windwidthsec = ((windowlength * fs)/2) * 60 + (60 / (4*fs))
    # go through each possible center value
    centr = DF.dt.cumsum().gt(windwidthsec).idxmax() + 1
    windows,centers = zip(*filter(lambda item: item is not None,[rangeGen(DF,center,windwidthsec,cutv,cutt,cutlength) for center in range(centr,DF.dt[::-1].cumsum().gt(windwidthsec).idxmax())]))
    return windows,centers

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
    return minimize(Likelihoodww(spd,hed,cv),pars,bounds=([0.01,None],[None,None],[None,None],[None,None],[None,None]),method="L-BFGS-B")

def yokoTate(optimAns,cv):
    yoko = Von_Mises_sd(np.sqrt(optimAns[2]*optimAns[2] + optimAns[1]*optimAns[1])) * Weibull_mean(optimAns[0],cv/sp.special.gamma(1+1/optimAns[0]))
    tate = Weibull_sd(optimAns[0], cv/sp.special.gamma(1 + 1/optimAns[0]))
    return yoko, tate

def ensureOptimConv(optimAns,spd,hed,cv):
    newAns = optimAns
    pars = [newAns.x[0],newAns.x[1],newAns.x[2],newAns.x[3],newAns.x[4]]
    # substitute first parameter if estimated < 0
    if newAns.x[0] >= 0.01:
        newAns = windOptims(spd,hed,cv,pars)
        return newAns, yokoTate(newAns.x,cv)
    else:
        return optimAns, yokoTate(optimAns.x,cv)

def windOptim(initHead,spd,head,cv):
    # run first optimisation
    _,_,par = initPars(initHead,spd,head,cv)
    answ = windOptims(spd,head,cv,par)
    # SD of heading vector perpendicular (yoko) and along (tate) mean direction
    yoko,tate=yokoTate(answ.x,cv)
    # repeat MLE to ensure convergence
    if 'tate' in locals():
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

def GOFtests(hed,nr,nd,answ,yoko,tate,cv):
    # mean track direction
    meangd = np.arctan2(np.mean(np.sin(hed)),np.mean(np.cos(hed)))
    nrp = sp.stats.kstest(nr,partial(sp.stats.weibull_min.cdf,c=answ.x[0],scale=cv/sp.special.gamma(1+1/answ.x[0])))
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
    
def stringify(vars):
    sep = ","
    return [str(x) for x in vars]

def maxLikeWind(r,d,cv):
    hd_try = 3
    answ_best = np.nan
    max_like = np.nan
    for id_hd in range(-hd_try,hd_try):
    
        # perform optimisation routine
        answ,yoko,tate = windOptim(id_hd,r,d,cv)

        if hasattr(answ,'success'):
            if (not np.isnan(tate)) & (answ.status==0):
        
                # calculate speeds and headings
                nr,nd = headSpdDir(r,d,answ)
        
                # GOF tests
                if np.isnan(max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
                    max_like = answ.fun
                    answ_best = answ

                if (answ.fun > max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
                    max_like = answ.fun
                    answ_best = answ

    return answ_best

def windEstimation(file,outfile,cutv: float=4.1667,cv=34.7/3.6,windowLength: int =51,rescaleTime: bool =True,isBp: bool = True):
    # wind estimation from bird GPS track
    # filename      - location of BiP formatted file
    # cutv          - minimum ground speed in m/s (default: 4.1667) 
    # cv            - mean air speed in m/s (default: 34.7 kph)
    # windowLength  - window size for estimation in mins (default: 51)
    # rescaleTime   - should original data be rescaled to 1 fix/min (default: True)

    # write the header of the outfile
    with open(outfile,'w',newline='') as f:
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(['Time','Lat','Lon','MeanHead','X','Y'])

    dat = prePare(file, rescaleTime,isBp)

    # generate windows over which estimation method will be run
    windows,centers = findWindows(dat,cutv,windowLength)


    # max likelihood calculations for wind estimation
    for win in range(len(windows)):
        
        r = dat.track_speed[windows[win]]
        d = dat.track_direction[windows[win]]
        answ_best = maxLikeWind(r,d,cv)

        if type(answ_best) != float:
            with open(outfile,'a',newline='') as f:
                writer = csv.writer(f, delimiter = ',')
                writer.writerow(stringify([(dat.DT[centers[win]]),dat.lat[centers[win]],dat.lon[centers[win]],np.arctan2(answ_best.x[2],answ_best.x[1]),answ_best.x[3],answ_best.x[4]]))

def roundTo(x,rn):
    return np.format_float_positional(x, precision=rn, unique=False, fractional=False, trim='k')

def makeForGraphText(num,rn):
    if (num < 0.1) & (num != 0):
        afterE = round(math.log10(num) - (0.5 if math.log10(num)<0 else 0))
        precede = roundTo(num * 10**abs(afterE),3)
        return str(precede) + "e" + str(afterE)
    else:
        return roundTo(num,rn)

# taken from 'https://gist.github.com/kn1cht/89dc4f877a90ab3de4ddef84ad91124e', circular python package by kn1cht
def Circcorrcoef(x, y, deg=True, test=False):
    '''Circular correlation coefficient of two angle data(default to degree)
    Set `test=True` to perform a significance test.
    '''
    convert = np.pi / 180.0 if deg else 1
    sx = np.frompyfunc(np.sin, 1, 1)((x - Circmean(x, deg)) * convert)
    sy = np.frompyfunc(np.sin, 1, 1)((y - Circmean(y, deg)) * convert)
    r = (sx * sy).sum() / np.sqrt((sx ** 2).sum() * (sy ** 2).sum())

    if test:
        l20, l02, l22 = (sx ** 2).sum(),(sy ** 2).sum(), ((sx ** 2) * (sy ** 2)).sum()
        test_stat = r * np.sqrt(l20 * l02 / l22)
        p_value = 2 * (1 - sp.stats.norm.cdf(abs(test_stat)))
        return tuple(round(v, 7) for v in (r, test_stat, p_value))
    return round(r, 7)

def Circmean(angles, deg=True):
    '''Circular mean of angle data(default to degree)
    '''
    a = np.deg2rad(angles) if deg else np.array(angles)
    angles_complex = np.frompyfunc(cmath.exp, 1, 1)(a * 1j)
    mean = cmath.phase(angles_complex.sum()) % (2 * np.pi)
    return round(np.rad2deg(mean) if deg else mean, 7)

# location of sample BiP data and resultant estimates CSV
if platform == "darwin":
    filename = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/SampleAxyTrek.csv"
    outfile = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/PyNewOut.txt"
else:
    filename = "I:/My Drive/PD/Data/TestingData/SampleAxyTrek.csv"
    outfile = "I:/My Drive/PD/Data/TestingData/PyNewOut.txt"

rTestFile = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/Rdataframe.txt"


dat = prePare(filename)


ogDat = pd.read_csv(filename)
ogDat['DT'] = pd.to_datetime(ogDat.time,format="%Y-%m-%d %H:%M:%S")

pyOut = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/PyOut.txt")
pyOut.Time = pd.to_datetime(pyOut.Time,format="%Y-%m-%d %H:%M:%S")
rOut = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/ROut.txt",on_bad_lines='skip',header=None)
rOut.rename(columns={0:'row',1:'DT',2:'lat',3:'lon',4:'meanHead',5:'X',6:'Y'},inplace=True)
rOut.DT = pd.to_datetime(rOut.DT,format="%Y-%m-%d %H:%M:%S")

nearestInd(RDat.DT,pd.to_datetime("2021-08-29 19:16:05",format="%Y-%m-%d %H:%M:%S"))

combDat = []
for x in range(len(rOut)):
    ind = nearestInd(pyOut.Time,rOut.DT[x])
    combDat.append([rOut.DT[x],rOut.X[x],pyOut.X[ind],rOut.Y[x],pyOut.Y[ind]])
combDat = pd.DataFrame(combDat)
combDat.rename(columns={0:'DT',1:'rX',2:'pyX',3:'rY',4:'pyY'},inplace=True)
plt.scatter(combDat.rX,combDat.pyX)
plt.plot([np.min(combDat.rX),np.max(combDat.rX)],[np.min(combDat.pyX),np.max(combDat.pyX)])
plt.show()
plt.scatter(combDat.rY,combDat.pyY)
plt.plot([np.min(combDat.rY),np.max(combDat.rY)],[np.min(combDat.pyY),np.max(combDat.pyY)])
plt.show()

plt.scatter(np.sqrt(combDat.rX**2 + combDat.rY**2),np.sqrt(combDat.pyX**2 + combDat.pyY**2))
plt.show()
windows=findWindows(pDat,34.7/3.6,51)

plt.scatter(np.arctan2(combDat.rY,combDat.rX),np.arctan2(combDat.pyY,combDat.pyX))
plt.show()
RDat.DT[windows[0]]
windEstimation(filename,outfile)

pyGen = pd.read_csv(outfile)
pyGen.Time = pd.to_datetime(pyGen.Time,format="%Y-%m-%d %H:%M:%S")
rGen = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/ROut.txt",header=None)
rGen.rename(columns={0:'row',1:'DT',2:'lat',3:'lon',4:'meanHead',5:'X',6:'Y'},inplace=True)
rGen.DT = pd.to_datetime(rGen.DT,format="%Y-%m-%d %H:%M:%S")

combDat = []
for x in range(len(rGen)):
    ind = nearestInd(pyGen.Time,rGen.DT[x])
    combDat.append([rGen.DT[x],rGen.X[x],pyGen.X[ind],rGen.Y[x],pyGen.Y[ind]])
combDat = pd.DataFrame(combDat)
combDat.rename(columns={0:'DT',1:'rX',2:'pyX',3:'rY',4:'pyY'},inplace=True)

plt.scatter(combDat.DT,combDat.rX,color='orange',label="R")
plt.scatter(combDat.DT,combDat.pyX,label="Py")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Estimated wind X component")
plt.show()

# compare X
plt.scatter(combDat.rY,combDat.pyY)
plt.xlabel("Time")
plt.ylabel("Estimated wind Y component")
plt.show()

props = dict(boxstyle='round', facecolor='lightgrey')

rAngCorr = Circcorrcoef(np.arctan2(combDat.rY,combDat.rX),np.arctan2(combDat.pyY,combDat.pyX),deg=False,test=True)
fig,ax=plt.subplots()
ax.scatter(np.arctan2(combDat.rY,combDat.rX),np.arctan2(combDat.pyY,combDat.pyX))
ax.plot(np.arange(-np.pi,np.pi),np.arange(-np.pi,np.pi),color='red',linestyle='dashed')
ax.set_title("Estimated wind heading comparison")
ax.set_xlabel("R generated wind heading (rad)")
ax.set_ylabel("Python generated wind heading (rad)")
ax.text(0.1, 0.9,
     "corr = " + makeForGraphText(rAngCorr[0],3) + "\n" + "p < 0.05",
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes,
     bbox = props)
plt.show()

Ypr = sp.stats.pearsonr(np.sqrt(combDat.rY**2+combDat.rX**2),np.sqrt(combDat.pyY**2 + combDat.pyX**2))
fig,ax=plt.subplots()
ax.scatter(np.sqrt(combDat.rY**2+combDat.rX**2),np.sqrt(combDat.pyY**2 + combDat.pyX**2))
ax.plot(np.arange(0,12),np.arange(0,12),color='red',linestyle='dashed')
ax.set_title("Estimated wind speed comparison")
ax.set_xlabel("R generated wind speed (m/s$^2$)")
ax.set_ylabel("Python generated wind speed (m/s$^2$)")
ax.text(0.1, 0.9,
     "corr = " + makeForGraphText(Ypr.statistic,3) + "\n" + "p = " + makeForGraphText(Ypr.pvalue,3),
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes,
     bbox = props)
plt.show()

fig, ax = plt.subplots()
ax.scatter(RDat.DT,RDat.speed,color='orange',label="R")
ax.plot(pDat.DT,pDat.track_speed,label="Py")
ax.legend()
ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
ax.set_ylabel("Bird ground speed (m/s$^2$)")
plt.show()

fig, ax = plt.subplots()
ax.scatter(RDat.DT,RDat['head'],color='orange',label="R")
ax.scatter(pDat.DT,pDat.track_direction,label="Py",s=.7)
ax.legend()
ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
ax.set_ylabel("Bird heading (rad)")
plt.show()

pyTest = pd.read_csv(outfile)

plt.scatter(RDat.DT,RDat.speed,color='orange')
plt.scatter(pDat.DT,pDat.track_speed,marker=2)
plt.show()


tdiffs = np.abs(RDat.DT - (pd.to_datetime("2021-08-24 19:43:01",format="%Y-%m-%d %H:%M:%S")))

RDat.DT[1499]
np.max(windows[0]) - np.min(windows[0])
plt.plot(RDat.track_speed[windows[0]])

windows = findWindows(pDat,4.1667,51)
out = [windMLE(df.track_speed[win],df.track_direction[win],34.7/3.6) for win in windows]

for win in windows:
    windMLE(df.track_speed[win],df.track_direction[win],34.7/3.6)

windMLE(df.track_speed[windows[0]],df.track_direction[windows[0]],34.7/3.6)

testr = [9.573340,6.870717,6.888984,9.775826,8.889476,9.087295,6.502933,7.621181,8.563742,6.304837,9.630127,9.503375,5.990518,10.615128,8.977472,11.682448,9.937076,8.986239,11.529867,12.037547,11.027331,11.958726,11.674352,9.236950,11.482676,11.099358,10.216385,10.680228,11.922089,11.138566,12.012096,10.430966,10.898318,11.373726,10.612603,12.276527,12.064162,9.506190,10.835805,10.837890,11.434897,11.542841,11.880499,11.840336,11.322806]
testd = [1.1928020,1.1245011,1.1341103,1.2624509,1.1234189,1.1661455,1.3566603,1.0684093,1.0587566,0.9347287,1.2419793,1.2388479,1.0015028,1.4442361,1.8459093,1.9182309,1.6511230,1.5017697,1.5017876,1.5411053,1.5031676,1.5000700,1.3827510,1.1922057,1.3152792,1.2855452,1.1636894,1.3092694,1.2440532,1.2413781,1.2151118,1.0770955,1.2063922,1.1382877,1.0852326,1.0732983,1.1655605,1.2568569,1.0391123,0.9735457,0.9969128,1.1527833,1.2271071,1.1667219,1.0243778]
np.sum(testd)
windMLE(testr,testd,34.7/3.6)



x = 1
for id_hd in range(-3,3):
    try:
        answ,yoko,tate = windOptim(id_hd,df.track_speed[windows[x]],df.track_direction[windows[x]],34.7/3.6)    
    except ValueError:
        pass
    

_,_,par = initPars(1,df.track_speed[windows[x]],df.track_direction[windows[x]],34.7/3.6)
    answ = windOptims(spd,head,cv,par)


np.array(df.track_speed[windows[x]])

hd_try = 3
# for id_hd in range(-hd_try,hd_try):

# perform optimisation routine        
answ,yoko,tate = windOptim(id_hd,testr,testd,cv)




if (not np.isnan(tate)) & (answ.status==0):

    # calculate speeds and headings
    nr,nd = headSpdDir(testr,testd,answ)

    # GOF tests
    if np.isnan(max_like) & (np.prod(GOFtests(testd,nr,nd,answ,cv)) == 1):
        max_like = answ.fun
        answ_best = answ

    elif (answ.fun > max_like) & (np.prod(GOFtests(d,nr,nd,answ,cv)) == 1):
        max_like = answ.fun
        answ_best = answ

    else:
        answ_best = np.nan

constv = 34.7/3.6
_,_,par = initPars(-3,testr,testd,constv)
answ = windOptims(testr,testd,constv,par)

    # then perform MLE for initial headings ranging from -3 to 3
hd_try = 3

# 976:1025
cutv = 4.1667
constv = 34.7/3.6
# first generate windows over which the estimation method will be run
windows = findWindows(df,51)
# set initial value for max likelihood and best answer


#5.692397 31.284772 -5.470312 -6.131341 11.298614
#[range(750, 801), range(751, 802), range(752, 803), range(753, 804), range(754, 805), range(755, 806), range(756, 807), range(757, 808), range(758, 809), range(759, 810)]
pars=[5.698634,31.316981,-5.485192,-6.125966,11.301431]

nr,nd = headSpdDir(df.track_speed[windows[0]],df.track_direction[windows[0]],testAns)
GOFtests(df.track_direction[windows[0]],nr,nd,testAns)

_,_,testPars = initPars(-3,df.track_speed[windows[0]],df.track_direction[windows[0]],constv)
testAns = windOptims(df.track_speed[windows[0]],df.track_direction[windows[0]],constv,testPars)
testYoko,testTate = yokoTate(testAns.x,constv)


# start from minimum possible point
fs = 1/np.abs(np.timedelta64(mode(np.diff(RDat.DT)),'m')).astype(int) # in Hz
expSamp = round(51 * fs) # expected number of samples
cutlength = round(45/51 * fs * 51)
error_of_sampling_interval = 5 * fs # give 5 seconds of leeway for samples to be found in
cutt = ((60 * fs) + error_of_sampling_interval).astype(int)

x=965
end2 = x+cutlength
(RDat.DT[x+expSamp] - RDat.DT[x]) <= np.timedelta64(cutt,'m'):
    if sum((RDat.track_speed[x:(x+expSamp)] > 4.1667) & (RDat.dt[x:(x+expSamp)] < cutt) & (RDat.track_direction[x:(x+expSamp)] != 100)) >= cutlength:
        if (RDat.DT[end2 - 1] - RDat.DT[x]) < np.timedelta64(cutt,'m'):
            while (RDat.DT[end2] - RDat.DT[x]) < np.timedelta64(cutt,'m'):
                end2 = end2 + 1
        return range(x,end2)
        else:
            return range(start,end)



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

RDat = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/Rdata.txt")
RDat.rename(columns = {'Unnamed: 0':'Index','Unnamed: 1':'DT','Unnamed: 2':'lat','Unnamed: 3':'lon','X':'X','Y':'Y','Unnamed: 6':'speed','Unnamed: 7':'head','Unnamed: 8':'distance'},inplace=True)
RDat.DT = pd.to_datetime(RDat.DT,format="%Y-%m-%d %H:%M:%S")

dat = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/SampleAxyTrek.csv")
dat = dat.loc[~np.isnan(dat.latitude),]
dat['DT'] = pd.to_datetime(dat['time'].str[0:-6], format = "%Y-%m-%d %H:%M:%S")
dat = timeRescale(dat,1)
dat = dat.iloc[1:].reset_index()


X,Y,zone,extra = utm.from_latlon(np.array(dat.latitude),np.array(dat.longitude))

combDat = []
for x in range(len(RDat)):
    ind = nearestInd(dat.DT,RDat.DT[x])
# inds = [nearestInd(dat.DT,RDat.DT[x]) for x in range(len(RDat))]
    combDat.append([RDat.DT[x],RDat.lat[x],RDat.lon[x],RDat.X[x],RDat.Y[x],RDat.speed[x],RDat['head'][x],dat.latitude[ind],dat.longitude[ind],X[ind],Y[ind],pDat.track_speed[ind],pDat.track_direction[ind]])
combDat = pd.DataFrame(combDat)
combDat.rename(columns={0:'DT',1:'rlat',2:'rlon',3:'rX',4:'rY',5:'rSpeed',6:'rDir',7:'pylat',8:'pylon',9:'pyX',10:'pyY',11:'pySpeed',12:'pyDir'},inplace=True)

plt.plot(combDat.rlat,combDat.pylat)
plt.show()

plt.scatter(combDat.DT,combDat.rX,color='orange')
plt.plot(combDat.DT,combDat.pyX)
plt.show()


plt.scatter(combDat.DT,combDat.rSpeed,color='orange')
plt.plot(combDat.DT,combDat.pySpeed)
plt.show()

spFromR = np.sqrt(np.diff(combDat.rX)**2 + np.diff(combDat.rY)**2)
spFromPy = np.sqrt(np.diff(combDat.pyX)**2 + np.diff(combDat.pyY)**2)
plt.scatter(combDat.DT.iloc[1:],spFromR,color='orange')
plt.plot(combDat.DT.iloc[1:],spFromPy)
plt.show()

plt.plot(combDat.rX-combDat.pyX)
plt.show()

pyResults = pd.read_csv(outfile)


combDat.rX-combDat.pyX

next(filter(lambda i: i[1] > 1, enumerate(combDat.rX-combDat.pyX)))

combDat.iloc[350:355]

pstart=prproj.Proj()
p=pyproj.Proj(proj='utm',zone=54,ellps='WGS84')
newx,newy = p(combDat.pylat[0],combDat.pylon[0])
plt.plot(combDat.rX - newx)
plt.show()

p = pyproj.Proj(proj='utm',zone=54,ellps='WGS84')
testx,testy = p(combDat.pylon, combDat.pylat)

plt.plot(combDat.rX - testx)
plt.show()

plt.hist(pDat.track_speed,1000)
plt.show()

def haversine(longitudes, latitudes):                                                                                  
    ''' Returns the distance between the points in m. Vectorized and will work with arrays and return an array of       
    distances, but will also work with scalars and return a scalar.'''
    lon1 = longitudes.values[:-1]                                                                                       
    lat1 = latitudes.values[:-1]                                                                                        
    lon2 = longitudes.values[1:]                                                                                        
    lat2 = latitudes.values[1:]                                                    
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    distance = 6371 * 2 * np.arcsin(np.sqrt(a)) * 1000
    return distance

dist = haversine(pDat.lon,pDat.lat)



rAngles = angles(RDat.lon,RDat.lat)

angles([-94.581213,-90.200203],[39.099912,38.627089])

lon1 = [-94.581213,-90.200203][0]                                                                                       
lat1 = [39.099912,38.627089][0]                                                                                        
lon2 = [-94.581213,-90.200203][1]                                                                                        
lat2 = [39.099912,38.627089][1]
lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])   

X = np.cos(lat2) * np.sin(abs(lon1 - lon2))
Y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(abs(lon1 - lon2))
np.arctan2(X,Y)


plt.plot(rAngles-RDat['head'][1:])
plt.show()

plt.scatter(pDat.dist.iloc[1:],dist)
plt.show()

rDist = haversine(RDat.lon,RDat.lat)

plt.scatter(rDist,RDat.distance[1:])
plt.show()

plt.scatter(gps_speed(RDat.lon,RDat.lat,RDat.DT),RDat.speed)
plt.show()




602765.105332 - 602732.104824
4.357874e+06 - 4.357927e+06 


plt.plot(combDat.rX - combDat.pyX)
plt.show()

import geopandas as gpd
from shapely.geometry import Point

s = gpd.GeoSeries([Point(x,y) for x, y in zip(dat.longitude, dat.latitude)])

geo_df = gpd.GeoDataFrame(df[['id', 'population']], geometry=s)
# Define crs for our geodataframe:
geo_df.crs = {'init': 'epsg:4326'} 


plt.plot(RDat.DT[0:100],df.DT[0:100])
plt.plot(df.DT[0:100],df.track_speed[0:100])
plt.show()


RDat = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/Rdataframe.txt")
# RDat.rename(columns = {'Unnamed: 0':'Index','Unnamed: 1':'DT','Unnamed: 2':'lat','Unnamed: 3':'lon','X':'X','Y':'Y','Unnamed: 6':'speed','Unnamed: 7':'head','Unnamed: 8':'distance'},inplace=True)
RDat.DT = pd.to_datetime(RDat.DT,format="%Y-%m-%d %H:%M:%S")

windows = findWindows(RDat,4.1667,51)
out = list()
for win in windows:
    out.append(windMLE(RDat.track_speed[win],RDat.track_direction[win],34.7/3.6))

x=0
r = RDat.track_speed[windows[x]]
d = RDat.track_direction[windows[x]]
cv = 34.7/3.6
# max likelihood calculations for wind estimation
# set the initial values for max likelihood
max_like = np.nan
answ_best = np.nan

hd_try = 3
# for id_hd in range(-hd_try,hd_try):
id_hd=-3

    # perform optimisation routine        
    answ,yoko,tate = windOptim(id_hd,r,d,cv)
    if (not np.isnan(tate)) & (answ.status==0):

        # calculate speeds and headings
        nr,nd = headSpdDir(r,d,answ)

        # GOF tests
        if np.isnan(max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
            max_like = answ.fun
            answ_best = answ

        elif (answ.fun > max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
            max_like = answ.fun
            answ_best = answ

        else:
            answ_best = np.nan

plt.plot(RDat.track_speed[win])
plt.show()


Xpr = sp.stats.pearsonr(combDat.rX,combDat.pyX)

def roundTo(x,rn):
    return np.format_float_positional(x, precision=rn, unique=False, fractional=False, trim='k')

def makeForGraphText(num,rn):
    if num < 0.1:
        afterE = round(math.log10(num) - (0.5 if math.log10(num)<0 else 0))
        precede = roundTo(num * 10**abs(afterE),3)
        return str(precede) + "e" + str(afterE)
    else:
        return roundTo(num,rn)

roundTo(Xpr.statistic,3)
makeForGraphText(Xpr.pvalue,3)


testX = [616346.3,616508.1,616525.7,616662.2,616896.6,617050.6,617262.6,617437.6,617612.4,617790.4,618021.1,618239.4,618275.7,618360.0,618576.5,618828.3,618885.4,619110.2,619296.7,619482.6,619610.9,619581.8,619518.4,619712.1,619792.5,619646.2,619407.5,619359.7,619396.9,619444.6,619466.0,619510.7,619561.4,619692.4,619897.2,620071.4,620258.8,620505.5,620668.4,620898.0,621117.8,621368.8,621660.5,621893.5,622179.5,622476.7,622828.2,623118.4,623291.6,623621.2,623993.0]
testY = [4374798,4375076,4375051,4375289,4375823,4376271,4376805,4377170,4377545,4378104,4378585,4379094,4379300,4379688,4380082,4380530,4380613,4380917,4381464,4382003,4382205,4382087,4382131,4382434,4383066,4383584,4384243,4384838,4385375,4386066,4386788,4387448,4388163,4388851,4389366,4390033,4390672,4391244,4391853,4392531,4393173,4393849,4394391,4395002,4395622,4396185,4396832,4397508,4398042,4398602,4399149]

#################################################################

# Run over all AxyTrek files

#################################################################

fileloc = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/2018Shearwater/AxyTrek/"
outf = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/PythonWinds/"
# list all files
files = glob.glob(fileloc + "**/*.txt")

for file in files:
    tag = re.search(re.escape(fileloc) + r".*/(.*?).txt", file).group(1)
    outfile = outf + tag + ".csv"

    windEstimation(files[0],outfile,isBp=False)