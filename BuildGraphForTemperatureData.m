function BuildGraphForTemperatureData()
% Build graph for temperature data
% 1. geodesic distance based approach
% 2. correlation based approach
close all; 
clear all;

nNeighbours = 6; % number of neighbors for each vertex
fracCoeff = 0.8;    % fraction of graph Fourier transform coefficients used in 
                    % restruction
nSites = 160;
nYears = 50;

load tempData.mat   % load variable 'tempData','siteInfo', 
                    % see ReadTemperatureFile.m for definition


% % % 在matlab地图上画出位置
% h = worldmap('China');
% landareas = shaperead('landareas.shp','UseGeoCoords',true);
% geoshow(landareas, 'FaceColor',[1 1 0.5]);

% 在google地图上画出位置
lat = siteInfo(2,:);
lon = siteInfo(3,:);
maph = figure(10); 
plot(lon,lat,'.k','MarkerSize',20)
plot_google_map;
title('Location of 160 cities');

%=====================================
%基于距离构建图
%=====================================
% 根据各个站点的经纬度计算相互之间的距离，得到160乘以160的距离矩阵 
distMt = zeros(nSites, nSites);
for i = 1:nSites
    for j = 1:nSites
        lat1 = siteInfo(2,i);
        lon1 = siteInfo(3,i);
        lat2 = siteInfo(2,j);
        lon2 = siteInfo(3,j);
        distMt(i,j) = GeoDistanceSimple(lat1, lon1, lat2, lon2);
    end
end
sigma_w = mean(distMt(:));  % sigma for weighting matrix, 
                            % using the mean distances between cities


% 每行选出最小的若干个距离,在标志矩阵上标记，用来画图
conInd = zeros(nSites, nSites); % conInd: connection indicator
for i = 1:nSites
    row = distMt(i,:);
    [dumy, idx]  = sort(row);
    for j = 1:nNeighbours % label the min nNeighbours, say, 6
        conInd(i,idx(j)) = 1;
        conInd(idx(j), i) = 1; % make it symmetric
    end
end

figure;
imagesc(1-conInd); colormap gray;

% 以中国地图为背景显示获得的图
lat = siteInfo(2,:);
lon = siteInfo(3,:);
figure(11);  hold on;
plot(lon,lat,'.r','MarkerSize',20)
for i = 1:nSites
    for j = 1:nSites
        if conInd(i,j) == 1 % 则从site i 到 site j 画一条直线
            line([siteInfo(3,i) siteInfo(3,j)], [siteInfo(2,i) siteInfo(2,j)]);
        end
    end
end
plot_google_map;

% 在图上画上2000年12月的温度
tempVec = zeros(nSites,1);
tempVec = tempData(12,nYears-1,:);
tempVec = tempVec(:);
stem3(lon,lat, tempVec, 'filled', 'sk', 'LineWidth', 2);
view(0,62);
title('Temperature for Dec. 2000');
%axis off;

% 计算加权矩阵矩阵，分解，获得基 
weightMt = exp(-distMt.^2 ./(2*sigma_w^2)).*(conInd); 
degVec = zeros(nSites,1);
degVec = sum(weightMt);
mtD = diag(degVec);
mtL = mtD - weightMt;
[V, D] = eig(mtL);
 
% 画出部分特征向量(比如前10个）
for k = 1:10
    figure;
    plot(lon,lat,'.r','MarkerSize',20);
    for i = 1:nSites
        for j = 1:nSites
            if conInd(i,j) == 1 % 则从site i 到 site j 画一条直线
                line([siteInfo(3,i) siteInfo(3,j)], [siteInfo(2,i) siteInfo(2,j)]);
            end
        end
    end
    hold on; 
    g = V(:,k);
    stem3(lon,lat, g, 'filled', 'k', 'LineWidth', 2);
    view(0,62);
    axis off;
end

% 变换图信号
f = tempVec; %这里暂时只是一个月的数据（one sample)
fTrans = zeros(nSites,1);
 for i = 1:nSites
    fTrans(i) = sum(f .* V(:,i));
 end
 figure; 
 stem(diag(D), fTrans,'filled','LineWidth',2);
 title('Transformed graph signal');
 figure;
 stem(diag(D), abs(fTrans),'filled','LineWidth',2);

%  % 重建
%  fRec = zeros(nSites,1);
%  for i = 1:floor(fracCoeff*nSites)
%     fRec = fRec + fTrans(i).*V(:,i);
%  end
%  sqrt(sum((f-fRec).^2))/sqrt(sum(f.^2))
%  

 % Batch Test
sigMt = zeros(nSites, nYears*12);
 for i = 1:nYears
     for j = 1:12
         temp = tempData(j,i,:); 
         sigMt(:, (i-1)*12+j) = temp(:);
     end
 end
%  err = CalReconErr(fracCoeff, V, D, sigMt);
%  disp('Relative recon. error for method 1');
%  err

% % 去掉均值
% for i = 1: size(sigMt,2)
%     sigMt(:,i) = sigMt(:,i) - mean(sigMt(:,i));
% end


coefMt = V'*sigMt;
save transformed.mat 'coefMt';

[binHeight, binCenter] = hist(coefMt(:),30);
figure; 
bar(binCenter, binHeight, 0.5);
 
ksdensity(coefMt(:));
figure; plot(coefMt(:));
 %===================================================================
 % Method 3: 用前20年数据，以温度距离作为连接的依据
 %===================================================================
 corrMat = zeros(nSites, nSites);
 diffMat = zeros(nSites, nSites);
 sigMat2 = zeros(nYears*12, nSites);
 for i = 1: nSites
    sigMat2(:,i) = reshape(tempData(:,:,i),nYears*12,1 );
 end
 corrMat = corr(sigMat2);
 for i = 1:nSites
     for j = 1:nSites
        diffMat(i,j) = sqrt(sum((sigMat2(:,i)-sigMat2(:,j)).^2))/(nYears*12);
     end
 end
 
 sigma_w = mean(abs(diffMat(:)));
  
% 每行选出最小的6个距离,在标志矩阵上标记，用来画图
conInd = zeros(nSites, nSites); % conInd: connection indicator
for i = 1:nSites
    row = diffMat(i,:);
    [dumy, idx]  = sort(row);
    for j = 1:nNeighbours % label the min nNeighbours
        conInd(i,idx(j)) = 1;
        conInd(idx(j), i) = 1; % make it symmetric
    end
end

figure;
imagesc(1-conInd); colormap gray;

% 在地图上画上连线
lat = siteInfo(2,:);
lon = siteInfo(3,:);
figure(22);  hold on;
plot(lon,lat,'.r','MarkerSize',20)
%set(gca,'LineWidth',6);
for i = 1:nSites
    for j = 1:nSites
        if conInd(i,j) == 1 % 则从site i 到 site j 画一条直线
            line([siteInfo(3,i) siteInfo(3,j)], [siteInfo(2,i) siteInfo(2,j)]);
        end
    end
end
plot_google_map;


distMtExp = exp(-diffMat.^2./(2*sigma_w^2));
weightMt = distMtExp .* conInd;

 
 % degree matrix 要改变
 mtD = diag(sum(weightMt));
 mtL = mtD - weightMt; 
[V, D] = eig(mtL);
coefMt2 = V'*sigMt;
save transformed2.mat 'coefMt2';
 
%======================================================================
% Calculate reconstruction error using given spectral basis
%======================================================================
 function err = CalReconErr(fracCoeff, V, D, sigMt)
 % <inputs>
 %  fracCoeff: fraction of coefficients used in reconstruction, 0<fracCoeff<1
 %  V: matrix of basis vectors, each column is a basis vector
 %  D: spectrum matrix (diagonal)
 %  sigMt: signal matrix, each column is a signal (length=nSites), If totally M
 %         columns, then the error result is averaged over M signals.
 %
 % <outputs>
 %  err: relative reconstruction error averaged over all signals.
 
 [lenSig, nSig] = size(sigMt);
 errEngVec = zeros(nSig,1); % error energy
 sigEngVec = zeros(nSig,1); % signal energy
 for i = 1:nSig
     f = sigMt(:,i);
     % decomposition
     fTran = f'*V;
     % retain only fracCoeff coeff.
     fTran(floor(fracCoeff*lenSig)+1:end) = 0;
     fTran = fTran(:);
     % reconstruction
     fRecon = V * fTran;
     errEngVec(i) = sum((fRecon-f).^2);
     sigEngVec(i) = sum(f.^2);
 end
 err = sum(sqrt(errEngVec))/sum(sqrt(sigEngVec));
