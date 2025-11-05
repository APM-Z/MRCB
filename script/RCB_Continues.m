% Simple script to plot RCBs
clear; close all; clc;
% User settings 
filename = 'E:\open_source\result\MRO1_GREJC_2_Continues_RCB.out';
system = 'Galileo';      % Choose: GPS, GLONASS, Galileo, QZSS, BDS-2, BDS-3
numFreq = 2;             % How many frequency bands to plot

% System configuration
bands.GPS = {'P1','P2','P5'}; cols.GPS = 12:14;
bands.GLONASS = {'R1','R2'}; cols.GLONASS = 15:16;
bands.Galileo = {'E1','E5a','E5b','E6','E5'}; cols.Galileo = 17:21;
bands.QZSS = {'P1','P2','P5'}; cols.QZSS = 22:24;
bands.BDS_2 = {'B1I','B3I','B2I'}; cols.BDS_2 = 25:27;
bands.BDS_3 = {'B1I','B3I','B2b','B1C','B2a'}; cols.BDS_3 = 25:29;

% Handle BDS names
if strcmp(system, 'BDS-2'), key = 'BDS_2';
elseif strcmp(system, 'BDS-3'), key = 'BDS_3';
else, key = system;
end

% Read data
data = readmatrix(filename, 'FileType','text', 'CommentStyle','#', 'Delimiter',' ', 'MultipleDelimsAsOne',true, 'NumHeaderLines',0);
hour = data(:,4);

% Get RCB data
RCB_data = data(:, cols.(key)) * 10/3;
band_names = bands.(key);

% Plot
num_plot = min(numFreq, length(band_names));
figure('Position', [100, 100, 1200, 600]);

for i = 1:num_plot
    plot(RCB_data(:, i), 'LineWidth', 2, 'DisplayName', [band_names{i}]);
    hold on;
end

% Set up time axis (0,6,12,18,24 hours)
total_pts = length(hour);

x1 = 1:484:length(hour);
x1 = [x1, total_pts];
xticks(x1);
hours = hour(1:484:end,:);
hours = [hours', 0];
xticklabels(hours);  
% Format plot
legend('Orientation','horizontal');
xlim([0 length(RCB_data(:, 1))]);
ylim([-8 8]);
xlabel('UTC [hour]','fontsize',20,'fontname','Times New Roman')	
ylabel('Variability (ns)','fontsize',20,'fontname','Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',20);
title(system,'fontsize',25,'fontname','Times');
grid on;