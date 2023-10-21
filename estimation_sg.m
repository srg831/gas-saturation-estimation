% estimate Sg by decline line
clc
clear
close all
%%% input data
% bm is bulk modulus
% sm is shear modulus
% phi is porosith
% rho is density
% IL is value of inverse line
load('welldata');
%%% Cubic spline interpolation for IL
m=length(phi);
m_sg=length(IL(1,:));
h=IL(1,m_sg)/(m-1)/10;
sg0=0:h:IL(1,m_sg);
decline=spline(IL(1,:),IL(2,:),sg0);
%%% calculate the vp in sediments with sg=0
Vp0=[];
for j=1:m
    vp1=biot2_vp_sg(phi(j),bm(j),sm(j),rho(j),vs(j));
    Vp0(j)=vp1(1);
end
%%% obtain the sg in free-gas bearing sediments
sg=[];
dec=[];
for i=1:m
    dd=100*(vp(i)-Vp0(i))/Vp0(i);
    if dd>0
        dd=-1.0e-9;
    end
    dec=[dec;dd];
    ddd=abs(decline-dd);
    II=find(ddd==min(ddd));
    sg(i)=sg0(II(end));
end
%%% plot 












