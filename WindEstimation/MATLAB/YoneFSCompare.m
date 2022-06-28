%% ANALYSE WIND AS USUAL
for b = 1:length(dat)
    [time, lat, lon] = gettimelatlon(dat, b);
%     forage = dat{b,6};
    [x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM
    DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
    tdiff = diff(time); % time difference between GPS points
    spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
    dir = atan2(diff(y),diff(x));
    [flight,fs,fe] = flightmask(spd,4,5);
    [ss,se] = getsection(1,300,60,fs,fe);
%     [vw,wd,va,resn] = windestimates(spd,dir,ss,se);
    [vw,wd,va,resnorm,bh,rwh,wInd] = windestimates5(spd,dir,ss,se);
    % create vector for wind values same length as lat lon
    wSpd = NaN(length(lat),1);
    wDir = NaN(length(lat),1);
    vA = NaN(length(lat),1);
    wSpd(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    vA(wInd(~isnan(wInd))) = va(~isnan(wInd));
    Resnorm = NaN(length(lat),1);
    Resnorm(wInd(~isnan(wInd))) = resnorm(~isnan(wInd)); 

    outW = table(time,lat,lon,wSpd,wDir,vA,Resnorm);
    % output the data
    writetable(outW, strcat(outloc,"1sFix/",tags(b),"WindYone.txt"));
end