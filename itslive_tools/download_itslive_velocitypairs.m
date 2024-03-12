% A script to load in and crop its_live netCDFs
clear all; close all

%% Download image pair velocities

%urlListFile = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_shearmargin_90p_30days/velocity_pairs_TEIS_all_shearmargin_90p_30days.txt';
urlListFile = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/velocity_pairs_TWITandTEIS.txt';

% Read the list of URLs from the text file
urlList = importdata(urlListFile);

% Specify the folder where you want to save the downloaded files
%downloadFolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_shearmargin_90p_30days/';
downloadFolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/';

% Loop through each URL and download the corresponding file
for i = 1:numel(urlList)
    url = urlList{i};
    
    % Extract the filename from the URL
    [~,filename,ext] = fileparts(url);
    
    % Create the full path to save the downloaded file
    fullFilePath = fullfile(downloadFolder, [filename, ext]);
    
    try
        % Use the webread function to download the file from the URL
        data = webread(url);
        
        % Write the downloaded data to the local file
        fid = fopen(fullFilePath, 'wb');
        fwrite(fid, data);
        fclose(fid);
        
        disp(['Downloaded: ' filename]);
    catch
        warning(['Failed to download: ' filename]);
    end
end

disp('Download completed.');