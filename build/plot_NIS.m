clear all
clc
close all

%%
filename = 'output.txt';
T = readtable(filename,'Delimiter','\t');
T.sensor_type = string(T.sensor_type);
idx_lidar = T.sensor_type == "lidar";
nis_lidar = T.NIS(idx_lidar);

histogram(nis_lidar,'normalization','cdf')
hline(0.95,'b-')
% vline(5.991,'r:')
% plot(nis_lidar)

hline(5.991,'r:')
figure
idx_radar = T.sensor_type == "radar";
nis_radar = T.NIS(idx_radar);
histogram(nis_radar,'normalization','cdf')
hline(0.95,'b-')
vline(7.815,'r:')
% plot(nis_radar)
% hline(7.815,'r:')