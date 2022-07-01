% LOAD DATA AND TAGNAMES
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
% LOAD IN (2016)
if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/";
else
    load('E:/My Drive/PhD/Data/2016Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
    outloc = "F:/UTokyoDrive/PhD/Data/2016Shearwater/WindEst/YoneMet/";
end

% ANALYSE WIND AS USUAL
orig = cell(length(dat),2);
% subsets (1s, 5s, 10s, 30s, 1min)
subs = [1, 5, 10, 30, 60];

if ismac()
    outloc = "/Volumes/GoogleDrive/My Drive/PD/BiP/YoneMethodComparison/";
else
    outloc = "I:/My Drive/PD/BiP/YoneMethodComparison/";
end
% create one with the same time window
for c = 1:length(tags)
    for b = 1:length(subs)
        datSub = dat;
        indeces = 1:subs(b):length(dat{c,3});
        datSub{c,3} = datSub{c,3}(indeces);
        datSub{c,4} = datSub{c,4}(indeces);
        datSub{c,5} = datSub{c,5}(indeces);
        
        [time, lat, lon] = gettimelatlon(datSub, c);
        %     forage = datSub{b,6};
        [x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM
        DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
        tdiff = diff(time); % time difference between GPS points
        spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
        dir = atan2(diff(y),diff(x));
        [flight,fs,fe] = flightmask(spd,4,5);
        [ss,se] = getsection(1/subs(b),300,60,fs,fe);
        %     [vw,wd,va,resn] = windestimates(spd,dir,ss,se);
        [vw,wd,va,resnorm] = windestimates(spd,dir,ss,se);
        % create vector for wind values same length as lat lon
        orig{b,1} = table(time(ss + round((se-ss)/2)),lat(ss + round((se-ss)/2)),...
            lon(ss + round((se-ss)/2)),vw,wd,va,resnorm,...
            repmat("5min",length(vw),1),repmat(string(subs(b)),length(vw),1),...
            'VariableNames', {'Time', 'lat','lon','wSpeed','wDir','aSpeed','resnorm','treatment','fs'});
        
        % try again but expanding the window size to always fit 300 points
        [ss,se] = getsection(1/subs(b),300*subs(b),60,fs,fe);
        
        %     [vw,wd,va,resn] = windestimates(spd,dir,ss,se);
        [vw,wd,va,resnorm] = windestimates(spd,dir,ss,se);
        % create vector for wind values same length as lat lon
        orig{b,2} = table(time(ss + round((se-ss)/2)),lat(ss + round((se-ss)/2)),...
            lon(ss + round((se-ss)/2)),vw,wd,va,resnorm,...
            repmat("300points",length(vw),1),repmat(string(subs(b)),length(vw),1),...
            'VariableNames', {'Time', 'lat','lon','wSpeed','wDir','aSpeed','resnorm','treatment','fs'});
    end
    writetable(vertcat(orig{:,1},orig{:,2}),strcat(outloc,tags(c),"ComparisonTable.txt"))
end

figure;
plot(table2array(orig{1,1}(:,1)),table2array(orig{1,1}(:,5)),'.')
hold on
for b = 2:length(subs)
    plot(table2array(orig{b,1}(:,1)),table2array(orig{b,1}(:,5)),'.')
end
legend;
hold off

for b = 1:length(subs)
    figure;
    plot(table2array(orig{b,1}(:,1)),table2array(orig{b,1}(:,5)),'.')
    hold on
    plot(table2array(orig{b,2}(:,1)),table2array(orig{b,2}(:,5)),'.')
    hold off
    legend;
end