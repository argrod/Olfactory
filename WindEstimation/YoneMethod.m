set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

if ismac()
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/AxyTrek/";
else
    fileloc = "E:/My Drive/PhD/Data/2016Shearwater/AxyTrek/";
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
    tagfiles = files(startsWith(filesNs,tags(tg))).name;
    fileID = fopen(strcat(fileloc,tags(tg),'/',tagfiles));
    GPSDat = textscan(fileID, '%{dd/MM/yyyy HH:mm:ss}D %f %f %f %f %f %f %f %s','Delimiter','\t','HeaderLines',1);
    fclose(fileID);
    dat{tg,3} = GPSDat{1};
    dat{tg,4} = GPSDat{2};
    dat{tg,5} = GPSDat{3};
end

[time, lat, lon] = gettimelatlon(dat, 1);
[x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM
DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
tdiff = diff(time); % time difference between GPS points
spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
dir = atan2(diff(y),diff(x));
% [spd,dir] = getspddir(lat,lon);

[flight,fs,fe] = flightmask(spd,4,5);
[ss,se] = getsection(1,300,60,fs,fe);

[vw,wd,va,re sn,rwh] = windestimates(spd,dir,ss(1:end-1),se(1:end-1));

plot(rwh/(180/pi))

for b = 405:length(se)
[vw,wd,va] = windestimates(spd,dir,ss(b),se(b));
end

ss(402:404)

[vw,wd,va] = windestimates(spd,dir,ss,se);
[wlat,wlon] = estimateposition(lat,lon,ss,se);

