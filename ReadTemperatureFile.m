function ReadTemperatureFileToMat()
% Read the temperatures of the 160 cites in China from 1951 to 2000, and write
% them into a .mat file for later use
% Use the matlab uiimport to generate the file reading script and read in
% the data to a matrix

clear all; close all; 
prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
cd(dir);
filename = [dir, '\Temperature\ch160temp']

%=================================================================
% 以下读取程序由matlab的uiimport的 generate script产生
delimiter = ' ';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
ch160temp = [dataArray{1:end-1}];
ch160temp = ch160temp';
clearvars filename delimiter formatSpec fileID dataArray ans;
%=============读取完毕==============================================

% 分离出站点位置 和 站点的数据 
% 设共有N个站点，共有50年数据，每年12个月。则
% （1）站点信息：序号 经度 纬度 保存为一个 3-by-N 的 矩阵，每列一个站点 
% （2）所有站点数据保存为数据体，每个slice对应一个站点，在每个Slice中，每列表示一年，共有50列
%       每行对应一个月。

nSites = 160;
nYears = 50;


tempSiteData = zeros(13, nYears+1, nSites);
tempSiteData = reshape(ch160temp, 13, nYears+1, nSites);
siteInfo = zeros(3,nSites);
tempData  = zeros(12, nYears, nSites); %

% 分离出站点位置  数据
for i = 1:nSites
   siteInfo(:,i) = tempSiteData(1:3,1,i);
end

% 分离出温度数据
for i = 1:nSites
    tempData(:,:,i) = tempSiteData(2:end,2:end, i)*0.1; %数据以0.1度为单位
end

save('tempData.mat','tempData','siteInfo');

