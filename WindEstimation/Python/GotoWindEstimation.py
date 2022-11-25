import scipy as sp
import numpy as np
import pandas as pd
import csv
import utm
from statistics import mode
import os, pyproj, math, datetime
from sys import platform
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from functools import partial

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
        print(L)
    return f

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

def distLatLon(lat,lon):
    X,Y,_,_ = utm.from_latlon(np.array(lat),np.array(lon))
    return np.append(np.nan,np.sqrt(np.diff(X)**2,np.diff(Y)**2))

def spTrav(DT,lat,lon,threshold=0): # distance and speed from time (DT), lat, and lon with option for a threshold speed
    dist = distLatLon(lat,lon)
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
    # copy of A1inv function from circular package in R
    if ((0 <= x) & (x < 0.53)):
        return 2 * x + pow(x,3) + (5 * pow(x,5))/6
    else:
        if (x < 0.85):
            return -0.4 + 1.39 * x + 0.43/(1-x)
        else:
            return 1/(pow(x,3) - 4 * pow(x,2) + 3 * x)

def windowFit(DF,start,end,cutT,cutV,cutLength):
    # generate windows capable of running the estimation method. Requirements are 51 minutes of data, with 75% of expected samples
    if (DF.DT[end] - DF.DT[start]) <= np.timedelta64(cutT,'m'):
        if sum((DF.track_speed[start:end] > cutV) & (DF.dt[start:end] < cutT) & (DF.track_direction[start:end] != 100)) >= cutLength:
            return np.array(range(start,end))[(DF.track_speed[start:end] > cutV) & (DF.dt[start:end] < cutT) & (DF.track_direction[start:end] != 100)]

def findWindows(DF,cutv,windowlength=51): # calculate appropriate windows from datetimes using minimum 
    # start from minimum possible point
    fs = 1/np.abs(np.timedelta64(mode(np.diff(DF.DT)),'m')).astype(int) # in Hz
    expSamp = round(51 * fs) # expected number of samples
    cutlength = round(45/51 * fs * 51)
    error_of_sampling_interval = 5 * fs # give 5 seconds of leeway for samples to be found in
    cutt = ((60 * fs) + error_of_sampling_interval).astype(int)
    return list(filter(lambda item: item is not None,[windowFit(DF,x,x+expSamp,cutt,cutv,cutlength) for x in range(len(DF) - expSamp)]))

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
    return minimize(Likelihoodww(spd,hed,cv),pars,bounds=((0,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf)),method="L-BFGS-B")

def yokoTate(optimAns,cv):
    yoko = Von_Mises_sd(np.sqrt(optimAns[2]*optimAns[2] + optimAns[1]*optimAns[1])) * Weibull_mean(optimAns[0],cv/sp.special.gamma(1+1/optimAns[0]))
    tate = Weibull_sd(optimAns[0], cv/sp.special.gamma(1 + 1/optimAns[0]))
    return yoko, tate

def ensureOptimConv(optimAns,spd,hed,cv):
    newAns = optimAns
    # for try_no in range(100):
    pars = [newAns.x[0],newAns.x[1],newAns.x[2],newAns.x[3],newAns.x[4]]
        # try:
    newAns = windOptims(spd,hed,cv,pars)
        #     if newAns.success:
        #         break
        # except:
        #     continue
    return newAns, yokoTate(newAns.x,cv)

def windOptim(initHead,spd,head,cv):
    # run first optimisation
    _,_,par = initPars(initHead,spd,head,cv)
    try:
        answ = windOptims(spd,head,cv,par)
        # SD of heading vector perpendicular (yoko) and along (tate) mean direction
        yoko,tate=yokoTate(answ.x,cv)
        # repeat MLE to ensure convergence
    except ValueError:
        pass
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

def windEstimation(file: str,outfile: str,cutv: float=4.1667,cv=34.7/3.6,windowLength: int =51,rescaleTime: bool =True):
    # wind estimation from bird GPS track
    # filename      - location of BiP formatted file
    # cutv          - minimum ground speed in m/s (default: 4.1667) 
    # cv            - mean air speed in m/s (default: 34.7 kph)
    # windowLength  - window size for estimation in mins (default: 51)
    # rescaleTime   - should original data be rescaled to 1 fix/min (default: True)

    dat = prePare(file, rescaleTime)

    # generate windows over which estimation method will be run
    windows = findWindows(dat,cutv,windowLength)

    # write the header of the outfile
    with open(outfile,'w',newline='') as f:
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(['Time','Lat','Lon','MeanHead','X','Y'])

    # max likelihood calculations for wind estimation
    for win in windows:

        # set the initial values for max likelihood
        max_like = np.nan
        answ_best = np.nan
        
        r = dat.track_speed[win]
        d = dat.track_direction[win]
        hd_try = 3
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

                    elif (answ.fun > max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
                        max_like = answ.fun
                        answ_best = answ

                    else:
                        answ_best = np.nan

        if type(answ_best) != float:
            with open(outfile,'a',newline='') as f:
                writer = csv.writer(f, delimiter = ',')
                writer.writerow(stringify([(RDat.DT[round((np.max(win)+np.min(win))/2)]),RDat.lat[round((np.max(win)+np.min(win))/2)],RDat.lon[round((np.max(win)+np.min(win))/2)],np.arctan2(answ_best.x[2],answ_best.x[1]),answ_best.x[3],answ_best.x[4]]))

# location of sample BiP data and resultant estimates CSV
if platform == "darwin":
    filename = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/SampleAxyTrek.csv"
    outfile = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/PyOut.txt"
else:
    filename = "I:/My Drive/PD/Data/TestingData/SampleAxyTrek.csv"
    outfile = "I:/My Drive/PD/Data/TestingData/PyOut.txt"

rTestFile = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/Rdataframe.txt"

windEstimation(filename,outfile)

pyTest = pd.read_csv(outfile)

tdiffs = np.abs(RDat.DT - (pd.to_datetime("2021-08-24 19:43:01",format="%Y-%m-%d %H:%M:%S")))

RDat.DT[1499]

windows = findWindows(df,4.1667,51)
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

plt.plot(RDat.DT[0:100],df.DT[0:100])
plt.plot(df.DT[0:100],df.track_speed[0:100])
plt.show()


RDat = pd.read_csv("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/Rdataframe.txt")
RDat.rename(columns = {'Unnamed: 0':'Index','Unnamed: 1':'DT','Unnamed: 2':'lat','Unnamed: 3':'lon','X':'X','Y':'Y','Unnamed: 6':'speed','Unnamed: 7':'head','Unnamed: 8':'distance'},inplace=True)
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