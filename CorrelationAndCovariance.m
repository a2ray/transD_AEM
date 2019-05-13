% clear all
% close all
% clc

Q = load('Residuals.mat');
ndata = size(Q.r,2);

Cr = corr(Q.r);
Cv = cov(Q.r);

Cv2 = zeros(ndata,ndata);
for j=1:ndata
    for k=j:ndata
        Cv2(j,k) = Cr(j,k)*sqrt(Q.sd(j))*sqrt(Q.sd(k));
        Cv2(k,j) = Cv2(j,k);
    end
end

figure(1)
imagesc(Cr)
colorbar
title('Data correlation matrix')

figure(2)
imagesc(log10(abs(Cv2)))
c=colorbar;
ylabel(c,'log10( T/s )')
title('Data covariance matrix')

BinEdges = -5:0.5:5;
y = [-5.5 BinEdges];
x = 1:size(Q.r,2)+1;
[xx,yy] = meshgrid(x,y);

figure(3)
hold on
for j=1:ndata
   rp(:,j) = histc(Q.r(:,j),BinEdges);
   plot(BinEdges,rp(:,j),'LineWidth',2)
end
hold off
xlabel('Normalized residual')
ylabel('Bin count')
title('Normalized residual histograms for all data points')

rp = [rp zeros(size(rp,1),1)];
rp = [zeros(1,size(rp,2)) ; rp];
figure(4)
pcolor(xx,yy,rp)
shading flat
xlabel('Data point number')
ylabel('Normalized residual')
title('Histograms of normalized residuals for each data point')

rpForAll = sum(rp,2);
rpForAll = rpForAll(1:end-1);
rpForAll = rpForAll/((BinEdges(3)-BinEdges(2))*sum(rpForAll));
stdG = (1/sqrt(2*pi)) * exp(-(BinEdges).^2);
figure(5)
plot(BinEdges,rpForAll,'LineWidth',3)
hold on
plot(BinEdges,stdG,'--','LineWidth',2)
xlabel('\sigma')
ylabel('probability density')
legend('estimated distribution','standard Gaussian')
title('Estimated data residual distribution')
