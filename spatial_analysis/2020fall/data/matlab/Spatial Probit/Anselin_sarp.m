% This file is written for the empirical studies using the data for the
% classical paper by Prof. Anselin

% These codes are from the Matlab manual book on Spatial Econometrics by
% Prof. James P. LeSage

clear
rng(20190611)
load anselin.data;
y = anselin(:,1); [n,junk] = size(y);
x = [ones(n,1) anselin(:,2:3)];
vnames = strvcat('crime','constant','income','hvalue'); 
load Wmat.dat; 
W = sparse(Wmat(:,1),Wmat(:,2),Wmat(:,3));
yc=zeros(n,1);
% now convert the data to 0,1 values
for i=1:n
  if y(i,1) > 40.0, yc(i,1) = 1; end;
end;
ndraw = 1100; nomit = 100;
prior.rval = 7; % logit estimates
results0 = sar(yc,x,W);
prt(results0,vnames);
result1 = sar_g(yc,x,W,ndraw,nomit,prior);
prt(result1,vnames);
% plt(result1,vnames);
% pause;
result2 = sarp_g(yc,x,W,ndraw,nomit,prior);
prt(result2,vnames);
% plt(result2,vnames);
prior.rval = 40; % probit estimates
result3 = sarp_g(yc,x,W,ndraw,nomit,prior);
prt(result3,vnames);
% plt(result3,vnames);

save Anselin_SARP.mat