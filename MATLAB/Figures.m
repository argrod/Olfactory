if ismac()
    outLoc = '/Volumes/GoogleDrive/My Drive/PhD/Figures/Olfactory/';
else
    outLoc = 'F:/UTokyoDrive/PhD/Figures/Olfactory/';
end
%% YONE METHOD EXAMPLE

%% READ IN
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

if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/MatlabDat/AxyTrek/ReadIn/AllGPSForage.mat');
    outloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/";
else
    load('F:/UTokyoDrive/PhD/Data/2019Shearwater/MatlabDat/AxyTrek/ReadIn/AllGPSForage.mat');
    outloc = "F:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/YoneMet/";
end

%% RUN WIND EST

b = 5;
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

%% EXTRACT A SECTION
% get variables
t = 50; % section number
vg = spd(ss(t):se(t)); % ground speed
gd = dir(ss(t):se(t)); % track direction
lat1 = lat(ss(t):se(t)); % latitude
lon1 = lon(ss(t):se(t)); % longitude

% plot track
plot(lon1,lat1)
axis equal
ylabel("Latitude")
xlabel("Longitude")
set(gcf,'color','white','units',"inches","position",[0,0, 5.2, 5.2])
set(gca,"fontsize",12,"FontName","Arial")

% create a colour gradient
len = 2*pi*10000;
red = [1, 0, 0];
blue = [0, 0, 1];
colors_grad = [linspace(red(1), blue(1), len)', linspace(red(2), blue(2), len)', linspace(red(3), blue(3), len)'];
dirCols = -pi:0.0001:pi;
gdCols = find(dirCols > (gd-0.00005) & dirCols < (gd+0.00005));

scatter(lon1,lat1,50,gd*(180/pi),"filled")
hcbar = colorbar;
colormap(jet(10));
caxis([-pi*(180/pi) pi*(180/pi)])


% plot track direction and ground speed on polar coordinate
gx = vg .* cos(gd);
gy = vg .* sin(gd);
plot(gx,gy,'.')
hold on
plot(0,0,'ok','MarkerFaceColor','k')
grid on
ax1 = gca;
ax1.GridColor = [0 0 0];
ax1.GridAlpha = 1;
hold off
axis equal

x = -pi:0.01:pi; % heading
x = reshape(x,[],1);

[c,resnorm]=wind2dveclsq(vg,gd,[3 0 9]);
cvec = c;
if cvec(1)<0
    cvec(1) = -cvec(1);
    cvec(2) = cvec(2) -pi;
end
cvec

% get the sum of the squared residuals
eavec = sqrt((gx-cvec(1)*cos(cvec(2))).^2 + (gy-cvec(1)*sin(cvec(2))).^2); % estimated air speed
resnormvec = sum((eavec - cvec(3)).^2) % set residuals as the deviance from the estimated air speed
% [xcos,ycos] = vecfit(ccos,x);
[xvec,yvec] = vecfit(cvec,x);

plot(gx,gy,'.k')
hold on
plot(xvec,yvec,'b','LineWidth',2)
% plot(xcos,ycos,'r--','LineWidth',2)
plot(0,0,'ok','MarkerFaceColor','k')
plot(cvec(1)*cos(cvec(2)),cvec(1)*sin(cvec(2)),'k')
quiver(0,0,cvec(1)*cos(cvec(2)),cvec(1)*sin(cvec(2)),'LineWidth',3,'MaxHeadSize',7)
hold off
grid on
ax1 = gca;
ax1.GridColor = [0 0 0];
ax1.GridAlpha = 1;
axis equal


%% try with 2016
if ismac()
    load('/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
else
    load('F:/UTokyoDrive/PhD//Data/2016Shearwater/MatlabDat/ReadIn/AllGPSForage.mat');
end

[time, lat, lon] = gettimelatlon(dat, 1);
[x,y,zone] = deg2utm(lat,lon); % convert from dec degs to UTM
DistTrav = sqrt(diff(x).^2+diff(y).^2); % calculate distance between GPS points
tdiff = diff(time); % time difference between GPS points
spd = DistTrav./seconds(tdiff); % speed travelled, m^-2
dir = atan2(diff(y),diff(x));
[flight,fs,fe] = flightmask(spd,4,5);
[ss,se] = getsection(1,300,60,fs,fe);

% get variables
t = 468; % section number
vg = spd(ss(t):se(t)); % ground speed
gd = dir(ss(t):se(t)); % track direction
lat1 = lat(ss(t):se(t)); % latitude
lon1 = lon(ss(t):se(t)); % longitude

% plot track
f = figure;
u = f.Units;
f.Units='inches';
f.Position = [10 10 3 3];
plot(lon1,lat1,'.')
axis equal
ylabel("Latitude")
xlabel("Longitude")
set(gcf,'color','white','units',"inches")
set(gca,"fontsize",10,"FontName","Arial")
export_fig(strcat(outLoc,"SegPlot.png"),'-r300')

% create a colour gradient
len = 2*pi*10000;
red = [1, 0, 0];
blue = [0, 0, 1];
colors_grad = [linspace(red(1), blue(1), len)', linspace(red(2), blue(2), len)', linspace(red(3), blue(3), len)'];
dirCols = -pi:0.0001:pi;
gdCols = find(dirCols > (gd-0.00005) & dirCols < (gd+0.00005));

f = figure;
u = f.Units;
f.Units='inches';
f.Position = [10 10 3 3];
scatter(lon1,lat1,10,gd*(180/pi),"filled")
hcbar = colorbar;
hcbar.Label.String = ['Heading (' char(176) ')'];
colormap(jet(10));
caxis([-pi*(180/pi) pi*(180/pi)])
ylabel("Latitude")
xlabel("Longitude")
set(gcf,'color','white','units',"inches")
set(gca,"fontsize",10,"FontName","Arial")
export_fig(strcat(outLoc,"HeadPlotCB.png"),'-r300')

% plot speed time series
plot(vg,'.-k');
xlabel('Time series','FontSize',12);
ylabel('Ground speed (m/s)','FontSize',12);

f = figure;
u = f.Units;
f.Units='inches';
f.Position = [10 10 3 3];
plot(gd,vg,'k.','MarkerSize',10)
xlim([-pi pi]);
xlabel('Track direction (radians)','FontSize',12);
ylabel('Ground speed (m/s)','FontSize',12);
set(gcf,'color','white','units',"inches")
set(gca,"fontsize",10,"FontName","Arial")
export_fig(strcat(outLoc,"SinePlot.png"),'-r300')

% plot track direction and ground speed on Cartesian coordinate
plot(gd,vg,'.','MarkerSize',12)
xlim([-pi pi]);
xlabel('Track direction (radians)','FontSize',10);
ylabel('Ground speed (m/s)','FontSize',12);

% plot track direction and ground speed on polar coordinate
gx = vg .* cos(gd);
gy = vg .* sin(gd);
plot(gx,gy,'.')
hold on
plot(0,0,'ok','MarkerFaceColor','k')
grid on
ax1 = gca;
ax1.GridColor = [0 0 0];
ax1.GridAlpha = 1;
hold off
axis equal
x = -pi:0.01:pi; % heading
x = reshape(x,[],1);

[c,resnorm]=wind2dveclsq(vg,gd,[3 0 9]);
cvec = c;
if cvec(1)<0
    cvec(1) = -cvec(1);
    cvec(2) = cvec(2) -pi;
end
cvec
% get the sum of the squared residuals
eavec = sqrt((gx-cvec(1)*cos(cvec(2))).^2 + (gy-cvec(1)*sin(cvec(2))).^2); % estimated air speed
resnormvec = sum((eavec - cvec(3)).^2) % set residuals as the deviance from the estimated air speed
% [xcos,ycos] = vecfit(ccos,x);
[xvec,yvec] = vecfit(cvec,x);

f = figure;
u = f.Units;
f.Units='inches';
f.Position = [10 10 3 3];
plot(gx,gy,'.k','MarkerSize',10);
hold on
plot(xvec,yvec,'b','LineWidth',2)
% plot(xcos,ycos,'r--','LineWidth',2)
plot(0,0,'ok','MarkerFaceColor','k')
plot(cvec(1)*cos(cvec(2)),cvec(1)*sin(cvec(2)),'k')
quiver(0,0,cvec(1)*cos(cvec(2)),cvec(1)*sin(cvec(2)),'LineWidth',3,'MaxHeadSize',7)
hold off
grid on
ax1 = gca;
ax1.GridColor = [0 0 0];
ax1.GridAlpha = 1;
axis equal
ylabel("Y component (m/s)")
xlabel("X component (m/s)")
set(gcf,'color','white','units',"inches")
set(gca,"fontsize",10,"FontName","Arial")
export_fig(strcat(outLoc,"CircFitPlot.png"),'-r300')
