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
% ���¶�ȡ������matlab��uiimport�� generate script����
delimiter = ' ';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
ch160temp = [dataArray{1:end-1}];
ch160temp = ch160temp';
clearvars filename delimiter formatSpec fileID dataArray ans;
%=============��ȡ���==============================================

% �����վ��λ�� �� վ������� 
% �蹲��N��վ�㣬����50�����ݣ�ÿ��12���¡���
% ��1��վ����Ϣ����� ���� γ�� ����Ϊһ�� 3-by-N �� ����ÿ��һ��վ�� 
% ��2������վ�����ݱ���Ϊ�����壬ÿ��slice��Ӧһ��վ�㣬��ÿ��Slice�У�ÿ�б�ʾһ�꣬����50��
%       ÿ�ж�Ӧһ���¡�

nSites = 160;
nYears = 50;


tempSiteData = zeros(13, nYears+1, nSites);
tempSiteData = reshape(ch160temp, 13, nYears+1, nSites);
siteInfo = zeros(3,nSites);
tempData  = zeros(12, nYears, nSites); %

% �����վ��λ��  ����
for i = 1:nSites
   siteInfo(:,i) = tempSiteData(1:3,1,i);
end

% ������¶�����
for i = 1:nSites
    tempData(:,:,i) = tempSiteData(2:end,2:end, i)*0.1; %������0.1��Ϊ��λ
end

save('tempData.mat','tempData','siteInfo');

