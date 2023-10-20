function vp=biot2_vp_sg(phi0,Ks0,nius0,rhos0,Vs0)
% two-phase biot motion equation
%%%%%%%%%%%%%%%%%%%
phi=phi0;
phis=1-phi0;%骨架占比
Sg0=0;
Sg=Sg0;%流体中气比例
%%%%%%% 组分数据
Ks=Ks0*1.0e9;%Gpa 骨架颗粒体积模量
nius=nius0*1.0e9;%Gpa 骨架颗粒剪切模量
rhow=1000;%Shenhu
rhos=(1000*rhos0-rhow*phi)/(1-phi);
Kw=2.32e9;%减小Shenhu
etaw=1.8e-3;%Shenhu
Kg=1.01325e3;%Shenhu
rhog=0.717;%Shenhu
etag=2.1e-5;
alpha=(1/phi-1)*(nius-rhos*Vs0*Vs0)/2/rhos/Vs0/Vs0;
w=8.0e3;%Shenhu
r=0.5; 
kappa=1e-5;
%%%%%%% 骨架模量公式
gamma=(1+2*alpha)/(1+alpha);
Km=Ks*(1-phi)/(1+alpha*phi);
nium=nius*(1-phi)/(1+alpha*gamma*phi);
%%%%%%%%%%%%
Jn=-1;
c=Km/Ks/phis;
c0=c/24.4;
Jnn=1.51+c0;
%%%%%%  流体模量公式
rr=power(Sg,Jnn);
etaf=etag*power(etaw/etag,1-rr);
rhof=rr*rhog+(1-rr)*rhow;
rhosat=phi*rhof+(1-phi)*rhos;
%%%% coefficient in Biot 应力系数
Kav=power((1-c)*phis*power(Ks,Jn)+phi*(1-rr)*power(Kw,Jn)+phi*rr*power(Kg,Jn),Jn);
A=((1-c)*phis)*((1-c)*phis)*Kav+Km-2*nium/3;
R=phi*phi*Kav;
Q=(1-c)*phis*phi*Kav;
N=nium;%Vs Qs 随含气量增大而增大
%%%%%%%%
%%%%%%%%% 密度系数
ainf=1+(1/phi-1)*r;
rho12=(1-ainf)*phi*rhof;
rho11=(1-phi)*rhos-rho12;
rho22=phi*rhof-rho12;
%%%%%
b=etaf*phi*phi/kappa;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=[rho11,rho12;rho12,rho22];
B=[b,-b;-b, b];
K=[A+N,Q;Q,R];
Niu=[N,0;0,0];
%
Lp = eig(w*w*rho+(1i)*w*B,K);
%Ls = eig(w*w*rho+(1i)*w*A,Niu,'chol');
Ls = eig(w*w*rho+(1i)*w*B,Niu);
%%% 保持正号
sigLp=sign(real(Lp));
sigLs=sign(real(Ls));
Lp=Lp.*sigLp;
Ls=Ls.*sigLs;
vp=w./real(sqrt(Lp));
vs=w./real(sqrt(Ls(1:2)));


