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
    [vw,wd,va,resnorm] = windestimates(spd,dir,ss,se);
    % create vector for wind values same length as lat lon
    wSpd = NaN(length(lat),1);
    wDir = NaN(length(lat),1);
    vA = NaN(length(lat),1);
    wSpd(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    vA(wInd(~isnan(wInd))) = va(~isnan(wInd));
    Resnorm = NaN(length(lat),1);
    Resnorm(wInd(~isnan(wInd))) = resnorm(~isnan(wInd)); 

    orig = table(time,lat,lon,wSpd,wDir,vA,Resnorm);
end

% subsets (1s, 5s, 10s, 30s, 1min)
subs = [5, 10, 30, 60];

% create one with the same time window

% create another with expanding time windows

wy([3 0]).*sin(gd) + wx([3 0]).*cos(gd) + ...
sqrt((wy([3 0]).*sin(gd)+wx([3 0]).*cos(gd)).^2 ...
    - wy([3 0]).^2 - wx([3 0]).^2 + 9^2) - vg