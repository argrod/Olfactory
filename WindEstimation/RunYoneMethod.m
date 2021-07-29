if ismac()
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/TxtDat/AxyTrek/BehaviourDetection/PredictedForage/";
else
    fileloc = "F:/UTokyoDrive/PhD/Data/2019Shearwater/TxtDat/AxyTrek/BehaviourDetection/PredictedForage/";
end
files = dir2(strcat(fileloc,"*ForageGPS.txt"));
%files = files(3:end);
% extract file names
filesNs = string(split([files.name],'.txt'));
filesNs = filesNs(contains(filesNs,"2019"));
% get file names and dates
tags = strings(length(filesNs),1);
for b = 1:length(filesNs)
    file = char(filesNs(b));
    tags(b) = file(1:strfind(file,'-2019')-1);
end
clear file i
% select only unique tag names and find the index of last occurrence
[tags,te] = unique(tags,'last');

dat = cell(length(tags),5);
for tg = 1:length(tags)
    tagfiles = filesNs(startsWith(filesNs,tags(tg)));
    for b = 1:length(tagfiles)
        fileID = fopen(strcat(fileloc,tagfiles(b),".txt"));
        GPSDat = textscan(fileID, '%{dd-MM-yyyy HH:mm:ss}D %f %f %f','Delimiter',',','HeaderLines',1);
        fclose(fileID);
        if b > 1
            dat{tg,3} = vertcat(dat{tg,3},GPSDat{1});
            dat{tg,4} = vertcat(dat{tg,4},GPSDat{2});
            dat{tg,5} = vertcat(dat{tg,5},GPSDat{3});
            dat{tg,6} = vertcat(dat{tg,6},GPSDat{4});
        else
            dat{tg,3} = GPSDat{1};
            dat{tg,4} = GPSDat{2};
            dat{tg,5} = GPSDat{3};
            dat{tg,6} = GPSDat{4};
        end
    end
end

%% READ IN STRAIGHT

if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/MatlabDat/AxyTrek/ReadIn/AllGPSForage.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/";
else
    load('F:/UTokyoDrive/PhD/Data/2019Shearwater/MatlabDat/AxyTrek/ReadIn/AllGPSForage.mat');
    outloc = "F:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/YoneMet/";
end

%% RUN WIND EST

for b = 1:length(dat)
    [time, lat, lon] = gettimelatlon(dat, b);
    forage = dat{b,6};
    [x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM
    DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
    tdiff = diff(time); % time difference between GPS points
    spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
    dir = atan2(diff(y),diff(x));
    [flight,fs,fe] = flightmask(spd,4,5);
    [ss,se] = getsection(.2,1500,60,fs,fe);
    [vw,wd,va,resn,bh,rwh,wInd] = windestimates5(spd,dir,ss,se);
    aveDir = zeros(length(forage),1);
    aveDir(wInd(~isnan(wInd))) = bh(~isnan(wInd));
    wDir = zeros(length(forage),1);
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    wSp = zeros(length(forage),1);
    wSp(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    % find the foraging points
    forSt = find(diff(forage) == 1) + 1;
    forEd = find(diff(forage) == -1);
    if forSt(1) > forEd(1)
        forSt = [1; forSt];
    end
    if forEd(end) < forSt(end)
        forEd = [forEd; length(forage)];
    end
    distTo = zeros(length(forage),1);
    for nxt = 1:length(forSt)
        if nxt == 1 && forSt(nxt) ~= 1
            distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
        else
            distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
        end
    end
    outW = table(time,lat,lon,forage,distTo,aveDir,wDir,wSp);
    % output the data
    writetable(outW, strcat(outloc,tags(b),"WindYone.txt"));
end

%% RUN THE 2016 DATA 
if ismac()
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/AxyTrek/";
else
    fileloc = "F:/UTokyoDrive/PhD/Data/2016Shearwater/AxyTrek/";
end
files = dir2(strcat(fileloc,"**/*.txt"));
%files = files(3:end);
% extract file names
filesNs = string(split([files.name],'.txt'));
filesNs = filesNs(contains(filesNs,"2016"));
% get file names and dates
tags = strings(length(filesNs),1);
for b = 1:length(filesNs)
    file = char(filesNs(b));
    tags(b) = file(1:strfind(file,'2016')-1);
end
clear file i
% select only unique tag names and find the index of last occurrence
[tags,te] = unique(tags,'last');

dat = cell(length(tags),5);
for tg = 1:length(tags)
    tagfiles = filesNs(startsWith(filesNs,tags(tg)));
    for b = 1:length(tagfiles)
        fileID = fopen(strcat(fileloc,tags(tg),"/",tagfiles(b),".txt"));
        GPSDat = textscan(fileID, '%{dd/MM/yyyy HH:mm:ss}D %f %f %f %f %f %f %f %f %f');
        fclose(fileID);
        if b > 1
            dat{tg,3} = vertcat(dat{tg,3},GPSDat{1});
            dat{tg,4} = vertcat(dat{tg,4},GPSDat{2});
            dat{tg,5} = vertcat(dat{tg,5},GPSDat{3});
        else
            dat{tg,3} = GPSDat{1};
            dat{tg,4} = GPSDat{2};
            dat{tg,5} = GPSDat{3};
        end
    end
end

%% LOAD IN (2016)
if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/";
else
    load('F:/UTokyoDrive/PhD//Data/2016Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "F:/UTokyoDrive/PhD//Data/2016Shearwater/WindEst/YoneMet/";
end

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
    [vw,wd,va,resn] = windestimates(spd,dir,ss,se);
    aveDir = zeros(length(forage),1);
    aveDir(wInd(~isnan(wInd))) = bh(~isnan(wInd));
    wDir = zeros(length(forage),1);
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    wSp = zeros(length(forage),1);
    wSp(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    % find the foraging points
    forSt = find(diff(forage) == 1) + 1;
    forEd = find(diff(forage) == -1);
    if forSt(1) > forEd(1)
        forSt = [1; forSt];
    end
    if forEd(end) < forSt(end)
        forEd = [forEd; length(forage)];
    end
    distTo = zeros(length(forage),1);
    for nxt = 1:length(forSt)
        if nxt == 1 && forSt(nxt) ~= 1
            distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
        else
            distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
        end
    end
    outW = table(time,lat,lon,forage,distTo,aveDir,wDir,wSp);
    % output the data
    writetable(outW, strcat(outloc,tags(b),"WindYone.txt"));
end


%% SUBSAMPLE

%% RUN 2017 DATA
if ismac()
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/TxtDat/BehaviourDetection/PredictedForage/";
else
    fileloc = "F:/UTokyoDrive/PhD/Data/2016Shearwater/TxtDat/BehaviourDetection/PredictedForage/";
end
files = dir2(strcat(fileloc,"*ForageGPS.txt"));
%files = files(3:end);
% extract file names
filesNs = string(split([files.name],'.txt'));
filesNs = filesNs(contains(filesNs,"2017"));
% get file names and dates
tags = strings(length(filesNs),1);
for b = 1:length(filesNs)
    file = char(filesNs(b));
    tags(b) = file(1:strfind(file,'-2017')-1);
end
clear file i
% select only unique tag names and find the index of last occurrence
[tags,te] = unique(tags,'last');

dat = cell(length(tags),5);
for tg = 1:length(tags)
    tagfiles = filesNs(startsWith(filesNs,tags(tg)));
    for b = 1:length(tagfiles)
        fileID = fopen(strcat(fileloc,tagfiles(b),".txt"));
        GPSDat = textscan(fileID, '%{dd-MM-yyyy HH:mm:ss}D %f %f %f','Delimiter',',','HeaderLines',1);
        fclose(fileID);
        if b > 1
            dat{tg,3} = vertcat(dat{tg,3},GPSDat{1});
            dat{tg,4} = vertcat(dat{tg,4},GPSDat{2});
            dat{tg,5} = vertcat(dat{tg,5},GPSDat{3});
            dat{tg,6} = vertcat(dat{tg,6},GPSDat{4});
        else
            dat{tg,3} = GPSDat{1};
            dat{tg,4} = GPSDat{2};
            dat{tg,5} = GPSDat{3};
            dat{tg,6} = GPSDat{4};
        end
    end
end

%% LOAD IN (2017)
if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/";
else
    load('F:/UTokyoDrive/PhD//Data/2017Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "F:/UTokyoDrive/PhD//Data/2017Shearwater/WindEst/YoneMet/";
end

%% RERUN 2016/17 WHILE SUBSAMPLED

