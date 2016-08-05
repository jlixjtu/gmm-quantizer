function p=img_csv_gaupdfexp3(ymean,ymin,ymax,r,fast)
% img_csv_pdfexp3 for optimal gaussian distribution quantizer design
% optimal quantizer cdf3 output for Gaussian point density function
% if already used once, then when fast=1, the function run in fast mode.
% and the gaussian4096cdf3.mat is refered.
% r=4.5 optimal empirical value (4~5) 
% Usage:
% p = img_csv_gaupdfexp3(d) returns an estmate of cdf^£¨1/3£©of the given Gaussian paramerter,
%             where p.cdf3 are the pdf^£¨1/3£© and cdf vectors estimated at
%             k resolutions, respectively and p.h is the bandwidth used for
%             the estimation. 
N=4096;% k=4096 max quantizer cdf3 resolution
yrange=ymax-ymin;
h=yrange/100;% bandwidth
ymax=ymax+3*h;
ymin=ymin-3*h;
dist=[abs(ymax)-abs(ymean),abs(ymin)-abs(ymean)];
sigma=max(dist)/r;
dx=(ymax-ymin)/(N-1);
if fast<1
p.x=ymin+(0:N-1)*dx;
p.pdf = pdf('norm',p.x,ymean,sigma);
p.pdf3=p.pdf.^(1/3);%CAUTION! NOT 3!
p.pdf3=p.pdf3./sum(p.pdf3.*dx);
p.cdf3=cumsum(p.pdf3.*dx);
% cr=range(p.cdf3);
% p.cdf3=p.cdf3/cr;
else
    load gaussian4096cdf3.mat 
    p.x=ymin+(0:N-1)*dx; %fast algorithm
end
end