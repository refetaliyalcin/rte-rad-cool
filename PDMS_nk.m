function [n,k]=PDMS_nk(x)
data=csvread('PDMS_nk.csv',1,0);
n = interp1(data(:,1),data(:,2),x);
k = interp1(data(:,1),data(:,3),x);