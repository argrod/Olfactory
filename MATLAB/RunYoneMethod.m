if ismac()
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/AxyTrek/";
else
    fileloc = "F:/UTokyoDrive/PhD/Data/2019Shearwater/AxyTrek/";
end
files = dir2(strcat(fileloc,"**/*.txt"));
%files = files(3:end);
% extract file names
filesNs = string();
for b = 1:length(files)
    filesNs(b) = files(b).name;
end
% filesNs = string(split([files.name],'.txt'));
% filesNs = filesNs(contains(filesNs,"2019"));
% get file names and dates
% tags = strings(length(filesNs),1);
% for b = 1:length(filesNs)
%     file = char(filesNs(b));
%     tags(b) = file(1:strfind(file,'2019')-1);
% end
% clear file i
tgs = dir2(fileloc);
tagsFolds = string();
for b = 1:length(tgs)
    tagsFolds(b) = tgs(b).name;
end
% select only unique tag names and find the index of last occurrence
tags = extractBefore(filesNs,".txt");

dat = cell(length(tags),5);
for tg = 1:length(tags)
    tagfiles = filesNs(startsWith(filesNs,tags(tg)));
    for b = 1:length(tagfiles)
        fileID = fopen(strcat(fileloc,tagsFolds(tg),"/",tagfiles(b)));
        GPSDat = textscan(fileID, '%{yyyy/MM/dd,HH:mm:ss}D %f %f %f %f %f %f %s','Delimiter','\t','HeaderLines',0);
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
    [flight,fs,fe] = flightmask(spd,5,1);
    [ss,se] = getsection(.2,300,60,fs,fe);
%     windestimates(spd,dir,ss,se)
    [vw,wd,va,resn,bh,rwh,wInd] = windestimates5(spd,dir,ss,se);
    aveDir = zeros(length(forage),1);
    aveDir(wInd(~isnan(wInd))) = bh(~isnan(wInd));
    wDir = zeros(length(forage),1);
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    wSp = zeros(length(forage),1);
    wSp(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    Resnorm = NaN(length(lat),1);
    Resnorm(wInd(~isnan(wInd))) = resn(~isnan(wInd)); 
    
    % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(time,lat,lon,aveDir,wDir,wSp,Resnorm);
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
        GPSDat = textscan(fileID, '%q %f %f %f %f %f %f %f %s','Delimiter','\t');
%         GPSDat = textscan(fileID, '%{dd/MM/yyyy HH:mm:ss}D %f %f %f %f %f %f %f %s');
%         textscan(fileID, '%{dd/MM/yyyy HH:mm:ss}D %f %f %f','Delimiter','\t','HeaderLines',1);
        fclose(fileID);
        if b > 1
            dat{tg,3} = vertcat(dat{tg,3},datetime(GPSDat{1},"InputFormat","dd/MM/yyyy HH:mm:ss"));
            dat{tg,4} = vertcat(dat{tg,4},GPSDat{2});
            dat{tg,5} = vertcat(dat{tg,5},GPSDat{3});
        else
            dat{tg,3} = datetime(GPSDat{1},"InputFormat","dd/MM/yyyy HH:mm:ss");
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
    outloc = "F:/UTokyoDrive/PhD/Data/2016Shearwater/WindEst/YoneMet/";
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
    
    % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(time,lat,lon,wSpd,wDir,vA,Resnorm);
    % output the data
    writetable(outW, strcat(outloc,"1sFix/",tags(b),"WindYone.txt"));
end


%% SUBSAMPLE

datSub = dat;
for b = 1:length(dat)
    indeces = 1:5:length(dat{b,3});
    datSub{b,3} = datSub{b,3}(indeces);
    datSub{b,4} = datSub{b,4}(indeces);
    datSub{b,5} = datSub{b,5}(indeces);

    [timeSub, latSub, lonSub] = gettimelatlon(datSub, b);
%     forage = datSub{b,6};
    [x,y,zone] = deg2utm(latSub,lonSub); % convert from dec degs to UTM
    DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
    tdiff = diff(timeSub); % time difference between GPS points
    spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
    dir = atan2(diff(y),diff(x));
    [flight,fs,fe] = flightmask(spd,5,1);
    [ssSub,seSub] = getsection(.2,300,60,fs,fe);
    [vwSub,wdSub,vaSub,resnSub,bh,rwh,wInd] = windestimates5(spd,dir,ssSub,seSub);
    wSpd = NaN(length(latSub),1);
    wDir = NaN(length(latSub),1);
    vA = NaN(length(latSub),1);
    wSpd(wInd(~isnan(wInd))) = vwSub(~isnan(wInd));
    wDir(wInd(~isnan(wInd))) = wdSub(~isnan(wInd));
    vA(wInd(~isnan(wInd))) = vaSub(~isnan(wInd));
    ResnSub = NaN(length(latSub),1);
    ResnSub(wInd(~isnan(wInd))) = resnSub(~isnan(wInd));

    % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(timeSub,latSub,lonSub,wSpd,wDir,vA,ResnSub);
    writetable(outW, strcat(outloc,"5sFix/",tags(b),"WindYone.txt"));
end

%% RUN 2017 DATA
if ismac()
    fileloc = '/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/AxyTrek/AxyTrek/';
else
    fileloc = 'F:/UTokyoDrive/PhD/Data/2017Shearwater/AxyTrek/AxyTrek/';
end
files = dir2(strcat(fileloc,"**/*.txt"));
%files = files(3:end);
% extract file names
filesNs = string(split([files.name],'.txt'));
filesNs = filesNs(1:end-1);
% get file names and dates
tagstruct = dir2(fileloc(1:end-1));
tags = strings(length(filesNs),1);
for b = 1:length(filesNs)
    tags(b) = string(tagstruct(b).name);
end
clear file i

dat = cell(length(tags),5);
for tg = 1:length(tags)
    fileID = fopen(strcat(fileloc,tags(tg),'/',filesNs(tg),".txt"));
    GPSDat = textscan(fileID, '%{yyyy/MM/dd,HH:mm:ss}D %f %f %f %f %f %f %f','Delimiter','\t','HeaderLines',1);
    fclose(fileID);
    GPSDat{1}.Format = "dd-MMM-yyyy HH:mm:ss";
    dat{tg,3} = GPSDat{1};
    dat{tg,4} = GPSDat{2};
    dat{tg,5} = GPSDat{3};
    dat{tg,6} = GPSDat{4};
end

%% LOAD IN (2017)
if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/";
else
    load('F:/UTokyoDrive/PhD//Data/2017Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "F:/UTokyoDrive/PhD//Data/2017Shearwater/WindEst/YoneMet/";
end

%% ANALYSE WIND AS USUAL
for b = 1:size(dat,1)
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
    
    % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(time,lat,lon,wSpd,wDir,vA,Resnorm);
    % output the data
    writetable(outW, strcat(outloc,"1sFix/",tags(b),"WindYone.txt"));
end


%% SUBSAMPLE

datSub = dat;
for b = 1:size(dat,1)
    indeces = 1:5:length(dat{b,3});
    datSub{b,3} = datSub{b,3}(indeces);
    datSub{b,4} = datSub{b,4}(indeces);
    datSub{b,5} = datSub{b,5}(indeces);

    [timeSub, latSub, lonSub] = gettimelatlon(datSub, b);
%     forage = datSub{b,6};
    [x,y,zone] = deg2utm(latSub,lonSub); % convert from dec degs to UTM
    DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
    tdiff = diff(timeSub); % time difference between GPS points
    spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
    dir = atan2(diff(y),diff(x));
    [flight,fs,fe] = flightmask(spd,5,1);
    [ssSub,seSub] = getsection(.2,300,60,fs,fe);
    [vwSub,wdSub,vaSub,resnSub,bh,rwh,wInd] = windestimates5(spd,dir,ssSub,seSub);
    wSpd = NaN(length(latSub),1);
    wDir = NaN(length(latSub),1);
    vA = NaN(length(latSub),1);
    wSpd(wInd(~isnan(wInd))) = vwSub(~isnan(wInd));
    wDir(wInd(~isnan(wInd))) = wdSub(~isnan(wInd));
    vA(wInd(~isnan(wInd))) = vaSub(~isnan(wInd));
    ResnSub = NaN(length(latSub),1);
    ResnSub(wInd(~isnan(wInd))) = resnSub(~isnan(wInd));

    % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(timeSub,latSub,lonSub,wSpd,wDir,vA,ResnSub);
    writetable(outW, strcat(outloc,"5sFix/",tags(b),"WindYone.txt"));
end

%% RUN 2014 DATA
if ismac()
    loc14 = '/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/AxyTrek/';
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/WindEst/YoneMet/";
else
    loc14 = 'F:/UTokyoDrive/PhD//Data/2014Shearwater/AxyTrek/';
    outloc = "F:/UTokyoDrive/PhD//Data/2014Shearwater/WindEst/YoneMet/";
end
files14 = dir2(strcat(loc14,'*.txt'));
files14 = string({files14.name})';
files14 = files14(~contains(files14,'yobi.txt'));
% remove tags with strange looking data
files14 = files14(~contains(files14,'st14n140f'));
files14 = files14(~contains(files14,'st14n64m'));
files14 = files14(~contains(files14,'st14n97f'));
dat = cell(length(files14),4);
for b = 1:length(files14)
    fileID = fopen(strcat(loc14,files14(b)));
    GPSDat = textscan(fileID, '%s %s %f %f %f %f %f %f %f %f',...
    'Delimiter',' ','HeaderLines',0);
    fclose(fileID);
    if b > 1
        dat{b,3} = vertcat(dat{b,3},datetime(cell2mat(strcat(GPSDat{1}," ",GPSDat{2})),'InputFormat',"HH:mm:ss.SS dd/MM/yyyy"));
        dat{b,4} = vertcat(dat{b,4},GPSDat{3});
        dat{b,5} = vertcat(dat{b,5},GPSDat{4});
    else
        dat{b,3} = datetime(cell2mat(strcat(GPSDat{1}," ",GPSDat{2})),'InputFormat',"HH:mm:ss.SS dd/MM/yyyy");
        dat{b,4} = GPSDat{3};
        dat{b,5} = GPSDat{4};
    end
end

save(strcat(outloc,"Dat14.mat"),"dat")


tags = erase(files14,'.txt');

%% LOAD IN 2014
if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD//Data/2014Shearwater/WindEst/YoneMet/Dat14.mat');
else
    load("F:/UTokyoDrive/PhD/Data/2014Shearwater/WindEst/YoneMet/Dat14.mat");
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
%     [vw,wd,va,resn] = windestimates(spd,dir,ss,se);
    [vw,wd,va,resnorm,bh,rwh,wInd] = windestimates5(spd,dir,ss,se);
    % create vector for wind values same length as lat lon
    wSpd = NaN(length(lat),1);
    wDir = NaN(length(lat),1);
    vA = NaN(length(lat),1);
    Resnorm = NaN(length(lat),1);
    wSpd(wInd(~isnan(wInd))) = vw(~isnan(wInd));
    wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
    vA(wInd(~isnan(wInd))) = va(~isnan(wInd));
    Resnorm(wInd(~isnan(wInd))) = resnorm(~isnan(wInd));
    
%     aveDir = zeros(length(forage),1);
%     aveDir(wInd(~isnan(wInd))) = bh(~isnan(wInd));
%     wDir = zeros(length(forage),1);
%     wDir(wInd(~isnan(wInd))) = wd(~isnan(wInd));
%     wSp = zeros(length(forage),1);
%     wSp(wInd(~isnan(wInd))) = vw(~isnan(wInd));
%     % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(time,lat,lon,wSpd,wDir,vA,Resnorm);
    % output the data
    writetable(outW, strcat(outloc,"1sFix/",tags(b),"WindYone.txt"));
end

%% SUBSAMPLE

datSub = dat;
for b = 1:length(dat)
    indeces = 1:5:length(dat{b,3});
    datSub{b,3} = datSub{b,3}(indeces);
    datSub{b,4} = datSub{b,4}(indeces);
    datSub{b,5} = datSub{b,5}(indeces);

    [timeSub, latSub, lonSub] = gettimelatlon(datSub, b);
%     forage = datSub{b,6};
    [x,y,zone] = deg2utm(latSub,lonSub); % convert from dec degs to UTM
    DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
    tdiff = diff(timeSub); % time difference between GPS points
    spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
    dir = atan2(diff(y),diff(x));
    [flight,fs,fe] = flightmask(spd,5,1);
    [ssSub,seSub] = getsection(.2,300,60,fs,fe);
    [vwSub,wdSub,vaSub,resnSub,bh,rwh,wInd] = windestimates5(spd,dir,ssSub,seSub);
    wSpd = NaN(length(latSub),1);
    wDir = NaN(length(latSub),1);
    vA = NaN(length(latSub),1);
    ResnSub = NaN(length(latSub),1);
    wSpd(wInd(~isnan(wInd))) = vwSub(~isnan(wInd));
    wDir(wInd(~isnan(wInd))) = wdSub(~isnan(wInd));
    vA(wInd(~isnan(wInd))) = vaSub(~isnan(wInd));
    ResnSub(wInd(~isnan(wInd))) = resnSub(~isnan(wInd));

    % find the foraging points
%     forSt = find(diff(forage) == 1) + 1;
%     forEd = find(diff(forage) == -1);
%     if forSt(1) > forEd(1)
%         forSt = [1; forSt];
%     end
%     if forEd(end) < forSt(end)
%         forEd = [forEd; length(forage)];
%     end
%     distTo = zeros(length(forage),1);
%     for nxt = 1:length(forSt)
%         if nxt == 1 && forSt(nxt) ~= 1
%             distTo(1:(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(1:(forSt(nxt) - 1))).^2 + (y(forSt(nxt)) - y(1:(forSt(nxt) - 1))).^2);
%         else
%             distTo((forEd(nxt-1) + 1):(forSt(nxt) - 1)) = sqrt((x(forSt(nxt)) - x(forEd(nxt-1)+1:(forSt(nxt)-1))).^2 + (y(forSt(nxt)) - y(forEd(nxt-1)+1:(forSt(nxt)-1))).^2);
%         end
%     end
    outW = table(timeSub,latSub,lonSub,wSpd,wDir,vA,ResnSub);
    writetable(outW, strcat(outloc,"5sFix/",tags(b),"WindYone.txt"));
end


%% NETCDF Files
if ismac()
    netloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/netCDF/ASCAT/2016/";
else
    netloc = "F:/UTokyoDrive/PhD/Data/netCDF/ASCAT/2016/";
end

netfiles = dir2(netloc);
fileTimes = cat(1,netfiles.name);
filesTimes = strcat(fileTimes(:,7:10),"/",fileTimes(:,11:12),"/",fileTimes(:,13:15),fileTimes(:,16:17),":",fileTimes(:,18:19),":",fileTimes(:,20:21));
filesTimes = datetime(filesTimes,'InputFormat','yyyy/MM/dd_HH:mm:ss');

lat = ncread(strcat(netloc,netfiles(12).name),"lat");
lon = ncread(strcat(netloc,netfiles(12).name),"lon");
lats = reshape(lat,[size(lat,1)*size(lat,2),1]);
lons = reshape(lon,[size(lon,2)*size(lon,2),1]);
times = ncread(strcat(netloc,netfiles(12).name),"time");
times = datetime(1990,01,01,0,0,0) + seconds(times);
times = reshape(times,[size(times,1)*size(times,2),1]);
% find unique times
times = times(1,:)';
range(times)
[x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM

timesall = dat{1,3};
latsall = dat{1,4};
lonsall = dat{1,5};
for b = 1:length(timesall)
    tpoint = timesall(b);
    if abs(filesTimes(find(filesTimes > tpoint,1)) - tpoint) < minutes(30)
    % find the nearest timepoint file
    ind = find(filesTimes > tpoint,1);
    ts = ncread(strcat(netloc,netfiles(ind).name),"time");
    ts = datetime(1990,01,01,0,0,0) + seconds(ts) + hours(9); % presumably in GMT, so add 9 hours
    lts = ncread(strcat(netloc,netfiles(ind).name),"lat");
    lns = ncread(strcat(netloc,netfiles(ind).name),"lon");
    wSp = ncread(strcat(netloc,netfiles(ind).name),"model_speed");
    wDir = ncread(strcat(netloc,netfiles(ind).name),"model_dir");
    % find the nearest time (column)
    timeind = find(ts(1,:)'>tpoint, 1);
    ts(1,:)
    [x,y,zone] = deg2utm(lts(:,timeind),lns(:,timeind)); % convert from dec degs to UTM
    [locx,locy,loczone] = deg2utm(latsall(b),lonsall(b));
    distFrom = sqrt((x - locx).^2 + (y - locy).^2)
    end
end
    
    
    ts(1,1)
    filesTimes(ind-1)
tst=cat(1,convertCharsToStrings(netfiles.name))
tst(,1)
[netfiles(1:end).name]
extractfield(netfiles,'name')

find(min(abs(latsall(b) - lts)) == min(min(abs(latsall(b) - lts))))
find(abs(latsall(b)-lts(:,339)) == min(abs(latsall(b)-lts(:,339))))

lts(29,337)
lts(23,339)
lns(29,337)
lns(23,339)
lonsall(b)
diff(lns(1,:))

latsall(b)

% presumably time is in GMT, so add 9 hours
ts = ts + hours(9);


% first select time that is appropriate


index = zeros(length(latsall),2);
for b = 1:length(latsall)
    index(b,:) = [find(abs(latsall(b)-lts(:,(min(abs(latsall(b) - lts)) == min(min(abs(latsall(b) - lts)))))) == min(abs(latsall(b)-lts(:,(min(abs(latsall(b) - lts)) == min(min(abs(latsall(b) - lts)))))))) find(min(abs(latsall(b) - lts)) == min(min(abs(latsall(b) - lts))))];
end
lts(29,337)
lts(105,339)

year = 2019;
start_jd = 182;
end_jd = 212;
for jd = start_jd:end_jd
    filename=['weather_data_' num2str(year) '_' num2str(jd) '.txt'];
    if isfile(filename)
        fprintf('already have the file |%s|\n',filename);
    else
        url=['http://www.weather.unh.edu/data/' num2str(year) '/' num2str(jd) '.txt'];
        outname=websave(filename,url);
        fprintf('got weather data file |%s|\n',outname);
    end
    
end

