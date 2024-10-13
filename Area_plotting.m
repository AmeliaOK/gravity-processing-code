%% Area plotting

close all; clear all;
format long 

data=xlsread('gravimetery_2022.xlsx');
x_coord = data(:,3);
y_coord = data(:,2);

%% Plot data points
figure(1);
plot(x_coord,y_coord,'g.', 'MarkerSize',40);
xlabel('UTM easting (m)');ylabel('UTM northing (m)');title('gravity point locations')
xlim([332550 332700])
ylim([5934100 5934190])
UalbertaArea(0.6)

%% Requires R2020a or later
exportgraphics(gcf,'vectorfig1.pdf','ContentType','vector')