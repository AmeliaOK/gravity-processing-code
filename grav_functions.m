function [East, North] = DMS2ugm(latD,latM,latS,lonD,lonM,lonS)
%finds utm easting and northing in meters from coordinates in UTM zone 12, northern hermisphere 
% in form of 6 columns, lattitude degrees,minutes, seconds, and longitude in degrees, minutes, 
% seconds.  Uses function wgs2utm by Alexandre Schimel

latDMS = [latD latM latS];
lonDMS = [lonD lonM lonS];

lat = dms2degrees(latDMS);
lon = -dms2degrees(lonDMS);

utmzone = 12; %set UTM zone
utmhemi = 'N'; %set hemisphere
[East,North] = wgs2utm(lat,lon,utmzone,utmhemi);
end

function stnmap = map_stns(x_coord,y_coord,stn)

%function that maps location of stations overtop of a map of the university
%of alberta. Uses UalbertaArea by Vadim Kravchinsky

stnlabel = int2str(stn); %convert station labels from intergers to strings to be used as labels on the figure

stnmap = figure(1);

plot(x_coord(1:18),y_coord(1:18),'b.', 'MarkerSize',18);
plot(x_coord(19:end),y_coord(19:end),'r.', 'MarkerSize',18);
xlabel('UTM Easting (m)');ylabel('UTM Northing (m)');title('Gravity Point Locations');
xlim([332550 332670]);
ylim([5934120 5934180]);
UalbertaArea(0.6);
text(x_coord,y_coord,stnlabel)

end

function parammap = map_param(East, North, Parameter, labels,stn)
%function that maps the interpolation of a parameter over a region defined
%by the northing and easting (in meters) of the locations where that
%parameter was recorded. Code adapted from Example2_map.m by Vadim
%Kravchinsky

%title = labels(1)
%unit = labels(2)

Nlin = linspace(min(North),max(North),100);
Elin = linspace(min(East),max(East),100);
% create mesh
method='linear';  %linear interpolation
[N,E]=meshgrid(Nlin,Elin);
% assign values of your Parameter to the values of your created mesh 
% in other words,buld the grid space
%Parameter_interpolated = griddata(North,East,Parameter,N,E);
Parameter_interpolated = griddata(East,North, Parameter,E,N);

index0 = find(stn<1);
parammap = figure(2);
set(gcf, 'Position',  [200, 200, 600, 600])
[c,h] = contourf(E,N,Parameter_interpolated); 
hold on;
plot(East, North,'.','MarkerEdgeColor','k' ,'MarkerSize',18)
plot(East(index0),North(index0),'o','MarkerEdgeColor','k' ,'MarkerSize',18)
xlabel('UTM Easting (m)'), ylabel('UTM Northing (m)'),
title(labels(1))
colorbar
set(get(colorbar,'Title'),'String',labels(2));
xlim([min(East)-7 max(East)+7]) %sets the range of the map, adapting automatically to different area ranges
ylim([min(North)-7 max(North)+7]) %sets the range of the figure, adapting automatically to different area ranges
UalbertaArea(0.6)
grid off
hold off;
end

function profile = makeprofile(Coordin, Parameter,labels)

profile = figure(3);
set(gcf, 'Position',  [50, 50, 500, 300])
plot(Coordin,Parameter,'-*');
grid on;
%axis ([45 46 30 31]) 
axis([min(Coordin), max(Coordin), min(Parameter), max(Parameter)])
    xlabel(labels(1));
    ylabel(labels(2));
    title(labels(3));

end

function tidcor_grav = tide_correct(tidefile,grav,time)

%%calculates the tidally corrected gravity from an initial set of relative gravity
%%readings using  a .dat file with the time in integer hours in the first column and tidal gravity
%in mGal in the second column. Created to be used with tidal correction
%files created by tidalcorrection.p by Vadim Kravchinsky. Requires recorded
%time in hours between 0 and 24
tide = load(tidefile);
tidetime = tide(:,1);
tidetime = tidetime + 6; %correct for difference in time zones between tidal correction data and experimental data
tidegrav = tide(:,2);

tidegravint = interp1(tidetime, tidegrav, time); %interpolate between tidal correction values according to experimental gravity recording times
tidcor_grav= grav - tidegravint; %preform and return the correctioned gravity values

end

function instr_grav = inst_correct(grav,stn,time)
%function that corrects for instrumental drift in gravity measurements,
%using sets of base station recordings and linear interpolation to estimate
%instrumental drift at each gravity recording. Requires time in hours
%between 0 and 24

index = find(stn<1); %find index of base station recordings
gravinst = interp1(time(index), grav(index), time(:)); %interpolate instrumental drift between base station recordings
instr_grav = grav(:)- gravinst; %preform and return the correctioned gravity values

end

function latgrav = latcorrect(North,grav)
%function that preforms lattitude correction for gravity measurements.
%Requires northing in meters
Go = 978032.7;% in mGal, base gravitational acceleration of the earth
a = 0.0053024;% coefficient, part of mathematical description of earth as an ellipsoid
Re = 6371000;% in meters, radius of the earth

bslat = North(1)/Re; % in radians, latitude of the base station from its northing in meters

latconst = Go*a*(sin(2*bslat))/Re; %in mGal/meter, coefficient of lattitude conversation from latitude of the basestation
deltaN =North(1)-North(:); %calculation of relative northing, in meters
latcorr = latconst*deltaN; %in mGal, lattitude correction
latgrav = grav + latcorr; %preform and return the correctioned gravity values

end

function airgrav = air_corr(elev, grav)
%function that calculates the free air gravity correction of the elevation
%in meters of the gravity measurements

deltaH = elev(:)-elev(1); % in meters, calculates relative elevation above base station
aircorrect = 0.3086*deltaH; %free air correction, in mGal
airgrav = grav+aircorrect; %correction adds when site is above base station

end

function rho = findrho(east,north,elev,grav)
%function that find the average density of material using the least squares
%method from easting in meters, northing in meters, elevation in meters and
%gravity in mGal at each location

%%create matrix for least squares
col1 = ones(length(grav),1);
M = [col1, east,north,elev]; 
%invert matrix
G = 6.6743*10^-11*10^5; %gravitational constant, adjusted for gravity values in mGal instead of m/s^2
A = inv(M'*M)*M'*grav; %A is a vector of coeffiecents given by the least squares method
rho = A(4)/(2*pi*G); %calculates least squares density from final coefficient in A vector, related to the bouguer correction
end 

function bouggrav = boug_corr(rho,elev,grav)
%complete bouguer correction for a supplied density, using elevation in
%meters and gravity in mGal
G = 6.6743*10^-11; %gravitational constant
deltaH = elev(:)-elev(1); %calculates relative elevation above base station in meters
bouguer = 2*pi*G*rho*deltaH*10^5; %finds bouguer correction
bouggrav = grav-bouguer; %correction is subtracted where the reading is above the base station

end

%%load inputs
%2022 data
grav22 = xlsread("gravimetery_2022.xlsx");
stn22 = grav22(:,1);
North22 = grav22(:,2);
East22 = grav22(:,3);
gravity22 = grav22(:,4);
elev22 = grav22(:,5);
time22= grav22(:,7);

%2024 data
grav24 = xlsread("gravdata_2024_3.xlsx");
stn = grav24(:,1);
latD = grav24(:,2);
latM = grav24(:,3);
latS = grav24(:,4);
lonD = grav24(:,5);
lonM = grav24(:,6);
lonS = grav24(:,7);
elev = grav24(:,8);
gravity = grav24(:,11);
time = grav24(:,14);


%%run functions
%%for 2022 data
tidcor_grav22 = tide_correct("2022TIDAL.dat",gravity22,time22);
instr_grav22 = inst_correct(tidcor_grav22,stn22,time22);
latgrav22 = latcorrect(North22,instr_grav22);
airgrav22 = air_corr(elev22, latgrav22);
rho22 = findrho(East22,North22,elev22,airgrav22);
bouggrav22 = boug_corr(rho22,elev22,airgrav22);

%for 2024 data
[East,North] = DMS2ugm(latD(:),latM(:),latS(:),lonD(:),lonM(:),lonS(:));
tidcor_grav6 = tide_correct("sept6TIDAL.dat",gravity(1:9),time(1:9));
tidcor_grav13 = tide_correct("sept13TIDAL.dat",gravity(10:18),time(10:18));
tidcor_grav = [tidcor_grav6;tidcor_grav13];
instr_grav = inst_correct(tidcor_grav,stn,time);
latgrav = latcorrect(North,instr_grav);
airgrav = air_corr(elev, latgrav);
rho24 = findrho(East,North,elev,airgrav);
bouggravrho24 = boug_corr(rho24,elev,airgrav);
bouggravrho22 = boug_corr(rho22,elev,airgrav);

%%combine files
totNorth = [North; North22];
totEast = [East;East22];
totlat = [latgrav;latgrav22];
totair = [airgrav;airgrav22];
totbou = [bouggravrho22;bouggrav22];
totstn = [stn;stn22];
reltop = [elev;elev22];

%figures for the actual report


bougravlabels = ["Plot of Bouguer corrected gravity anomaly","mGal"];
airgravlabels = ["Plot of free air corrected gravity anomaly","mGal"];
toplabels = ["Elevation above sea level","m"];

%profile figures
%topography
toplabels22 = ["UTM Northing (m)","Elevation relative to base station (m)", "Topography along North South profile of 2022 gravity data"];
figure(4) = makeprofile(North22(2:6),elev22(2:6)-elev22(1), toplabels22);
exportgraphics(gcf,'22topprofile.pdf','ContentType','vector');

northing = [North(13:15);North(6);North(16:17)];
toponorth = [elev(13:15);elev(6);elev(16:17)];
toplabels24N = ["UTM Northing (m)","Elevation relative to base station (m)", "Topography along North South profile of 2024 gravity data"];
figure(5) = makeprofile(northing,toponorth-toponorth(1), toplabels24N);
exportgraphics(gcf,'24Ntopprofile.pdf','ContentType','vector');

easting = [East(2:8);East(11:12)];
topoeast = [elev(2:8);elev(11:12)];
toplabels24W = ["UTM Easting (m)","Elevation relative to base station (m)", "Topography along East West profile of 2024 gravity data"];
figure(6) = makeprofile(easting,topoeast-topoeast(1), toplabels24W);
exportgraphics(gcf,'24Etopprofile.pdf','ContentType','vector');

%free air gravity
airlabels22 = ["UTM Northing (m)","Gravity Anamoly (mGal)", "Free air corrected gravity anamoly along profile of 2022 data"];
figure(7) = makeprofile(North22(2:6),airgrav22(2:6), airlabels22);
exportgraphics(gcf,'22freeairprofile.pdf','ContentType','vector');

northing = [North(13:15);North(6);North(16:17)];
airgravnorth = [airgrav(13:15);airgrav(6);airgrav(16:17)];
airlabels24N = ["UTM Northing (m)","Gravity Anamoly (mGal)", "Free air corrected gravity anamoly along profile of 2024 data"];
figure(8) = makeprofile(northing,airgravnorth, airlabels24N);
exportgraphics(gcf,'24Nfreeairprofile.pdf','ContentType','vector');

easting = [East(2:8);East(11:12)];
airgraveast = [airgrav(2:8);airgrav(11:12)];
airlabels24W = ["UTM Easting (m)","Gravity Anamoly (mGal)", "Free air corrected gravity anamoly along profile of 2024 data"];
figure(9)= makeprofile(easting,airgraveast, airlabels24W);
exportgraphics(gcf,'24Efreeairprofile.pdf','ContentType','vector');

%Bouguer corrected
bouglabels22 = ["UTM Northing (m)","Gravity Anamoly (mGal)", "Bouguer corrected gravity anamoly along profile of 2022 data"];
figure(10)=makeprofile(North22(2:6),bouggrav22(2:6), bouglabels22);
exportgraphics(gcf,'22bougprofile.pdf','ContentType','vector');

northing = [North(13:15);North(6);North(16:17)];
bouggravnorth = [bouggravrho24(13:15);bouggravrho24(6);bouggravrho24(16:17)];
bouglabels24N = ["UTM Northing (m)","Gravity Anamoly (mGal)", "Bouguer corrected gravity anamoly along profile of 2024 data"];
figure(11) = makeprofile(northing,bouggravnorth, bouglabels24N);
exportgraphics(gcf,'24Nbougprofile.pdf','ContentType','vector');

easting = [East(2:8);East(11:12)];
bouggraveast = [bouggravrho24(2:8);bouggravrho24(11:12)];
bouglabels24W = ["UTM Easting (m)","Gravity Anamoly (mGal)", "Bouguer corrected gravity anamoly along profile of 2024 data"];
figure(12) = makeprofile(easting,bouggraveast, bouglabels24W);
exportgraphics(gcf,'24Ebougprofile.pdf','ContentType','vector');

%2d figures
figure(13) = map_param(totEast, totNorth, totair, airgravlabels,totstn);
exportgraphics(gcf,'totairgrav.pdf','ContentType','vector');
figure(14) = map_param(totEast, totNorth, totbou, bougravlabels,totstn);
exportgraphics(gcf,'totbouggravrho22.pdf','ContentType','vector');
figure(15) = map_param(totEast, totNorth, reltop, toplabels,totstn);
exportgraphics(gcf,'totreltopo.pdf','ContentType','vector');
figure(16) = map_stns(totEast,totNorth,totstn);
exportgraphics(gcf,'stations.pdf','ContentType','vector');