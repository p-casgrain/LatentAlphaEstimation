addpath('../../Data/') % Add Data Folder to the Path

files = dir('../../Data/*.mat'); % List of .mat files in the directory sorted by name

time_step = 1000; % Time Step in ms -> 1000 corresponds to 1 sec.

LOB = cell(numel(files),1);
h = waitbar(0,'Starting Import');

for i=1:numel(files)
    load(files(i).name); % File loaded into memory with variable name 'data'
    LOB{i} = ExtractLOBSnapshots(data,time_step,1); % Extract touch price data
    clear data
    
    progress = i/numel(files);
    waitbar(progress,h,sprintf('Extracting Data - %d%% Complete',round(progress*100)));
end
clear data; % Remove 'data' item from memory

save('LOB_Parsed_Data.mat',LOB)

close(h);
delete(h);

clear all;



