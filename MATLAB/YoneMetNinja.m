if ismac()
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/GypsyTranslocation/";
else
    fileloc = "F:/UTokyoDrive/PhD/Data/2018Shearwater/GypsyTranslocation/";
end
files = dir2(strcat(fileloc,"*.csv"));
%files = files(3:end);
% extract file names
filesNs = string();
for b = 1:length(files)
    filesNs(b) = files(b).name;
end
tags = extractBetween(filesNs,10,15);

dat = cell(length(tags),5);
for tg = 1:length(tags)
    fileID = fopen(strcat(fileloc,filesNs(tg)));
    GPSDat = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter',',','HeaderLines',0);
    fclose(fileID);
    GPSDat{16} = datetime(GPSDat{1},GPSDat{2},GPSDat{3},GPSDat{4},GPSDat{5},GPSDat{6});
    if tg < 6
        ind = find(GPSDat{16} > datetime(2018,09,02,04,0,0),1);
    else
        ind = find(GPSDat{16} > datetime(2018,09,01,04,0,0),1);
    end
    dat{tg,3} = GPSDat{16}(ind:end);
    dat{tg,4} = GPSDat{7}(ind:end);
    dat{tg,5} = GPSDat{8}(ind:end);
end
save('/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/GypsyTranslocation/AllGPS.mat','dat')

%% READ IN STRAIGHT

if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/GypsyTranslocation/AllGPS.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/GypsyTranslocation/WindEst/YoneMet/";
else
    load('F:/UTokyoDrive/PhD/Data/2018Shearwater/GypsyTranslocation/AllGPS.mat');
    outloc = "F:/UTokyoDrive/PhD/Data/2018Shearwater/GypsyTranslocation/WindEst/YoneMet/";
end

%% RUN WIND EST

for b = 1:length(dat)
    [time, lat, lon] = gettimelatlon(dat, b);
    [x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM
    DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
    tdiff = diff(time); % time difference between GPS points
    spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
    dir = atan2(diff(y),diff(x));
    [flight,fs,fe] = flightmask(spd,4,5);
    [ss,se] = getsection(5,300,60,fs,fe);
%     windestimates(spd,dir,ss,se)
    [vw,wd,va,resn,bh,rwh,wInd] = windestimates5(spd,dir,ss,se);
    aveDir = zeros(length(lat),1);
    aveDir(wInd(~isnan(wInd))) = bh(~isnan(wInd));
    wDir = zeros(length(lat),1);
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    wSp = zeros(length(lat),1);
    wSp(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    Resnorm = NaN(length(lat),1);
    Resnorm(wInd(~isnan(wInd))) = resn(~isnan(wInd));     
    outW = table(time,lat,lon,aveDir,wDir,wSp,Resnorm);
    % output the data
    writetable(outW, strcat(outloc,tags(b),"WindYone.txt"));
end
