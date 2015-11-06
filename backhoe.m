% this program deals with true data of backhoe,and it is modified
% version of backhoe.m ,however,it is not correct completely
%use polar reformatting method ;elevation angle :0
% azimuth angle 350--100
% this is self-programmed isar imaging of normal target
% it is large angles,wide bandwidth 
% used data is backhoe:Matlab data file for elevation zero, 
% azimuth 350 degrees to 100 degrees to 100 degrees
% Row which has five * is the problem I am confronted with and important 

%%
close all
clear all
load('backhoe_el000_az350to100.mat');
c = 3e8; % speed of light
fc = 10e9; % center frequency
lamadac=c/fc;
phic = (-10+100)/2*pi/180; % center of azimuth look angles
% cut out part of the azimuth angle
N=size(AZ,2);
% index=floor(N*6/11):floor(N*7/11);
% index=floor(N*5/11):floor(N*6/11);
% index=floor(N/3):floor(N*2/3);
index=floor(1):floor(N);
AZ1=AZ(index);
%AZ1=AZ(N/3:N*2/3);
Nr=size(FGHz,1);Na=size(AZ1,2);
% for SAR ,bx is the resolution of azimuth,by the range
df=(FGHz(2)-FGHz(1))*1e9;
by=c/df/2;
%by=by*1.5;
f=fc+[-df*Nr/2:df:df*(Nr/2-1)];
k=2*pi*f/c;% for SAR data ,k=4*pi*f/c
dphi=(AZ(2)-AZ(1))*pi/180;
phi=phic+[-dphi*Na/2:dphi:dphi*(Na/2-1)];
bx=lamadac/dphi/2;
dx=bx/Na;% resolution cell
dy=by/Nr;
x=-bx/2+(0:Na-1)*dx;% cross_range vector
y=-by/2+(0:Nr-1)*dy;
%Es=VV(:,index);
Es=HH(:,index);
%% APPLY 2D IFT to a wideband and large angle 
cro_pro1=fftshift(ifft2(Es));figure;
matplotcorrect(x,-y,abs(cro_pro1),40);
plotti(x,-y,abs(cro_pro1),40);
xlabel('azimuth');ylabel('range');
colormap(1-gray);
colorbar
title('imaging by 2D FFT');
%% polar reformatting method
kx=k'*cos(phi);ky=k'*sin(phi);% compute max and min kx and ky firstly
kxmax=max(max(kx));
kxmin=min(min(kx));
kymax=max(max(ky));
kymin=min(min(ky));
dkx=(kxmax-kxmin)/(Na-1);
dky=(kymax-kymin)/(Nr-1);
kxnew=kxmin:dkx:kxmax;
kynew=kymin:dky:kymax;

kxnew(Na+1)=0;% x is corresponding with Na
kynew(Nr+1)=0;

m=0;n=0;
es=zeros(Nr+1,Na+1);
for u=k
    n=n+1;
    m=0;
    for v=phi
        m=m+1;
        ux=u*cos(v);% corresponding with Na
        uy=u*sin(v);
        inx=floor((ux-kxmin)/dkx)+1;% 
        iny=floor((uy-kymin)/dky)+1;
        
        r1=sqrt(abs(kxnew(inx)-ux)^2+abs(kynew(iny)-uy)^2);
        r2=sqrt(abs(kxnew(inx)-ux)^2+abs(kynew(iny+1)-uy)^2);
        r3=sqrt(abs(kxnew(inx+1)-ux)^2+abs(kynew(iny+1)-uy)^2);
        r4=sqrt(abs(kxnew(inx+1)-ux)^2+abs(kynew(iny)-uy)^2);
        
        R=1/r1+1/r2+1/r3+1/r4;
        
        es(iny,inx)=es(iny,inx)+Es(n,m)/r1/R;
        es(iny+1,inx)=es(iny+1,inx)+Es(n,m)/r2/R;
        es(iny+1,inx+1)=es(iny+1,inx+1)+Es(n,m)/r3/R;
        es(iny,inx+1)=es(iny,inx+1)+Es(n,m)/r4/R;
%         es(inx,iny)=es(inx,iny)+Es(n,m)/r1/R;
%         es(inx+1,iny)=es(inx+1,iny)+Es(n,m)/r2/R;
%         es(inx,iny+1)=es(inx,iny+1)+Es(n,m)/r3/R;
%         es(inx+1,iny+1)=es(inx+1,iny+1)+Es(n,m)/r4/R;
    end
end
esnew=es(1:Nr,1:Na);
% kxnew1=kxnew(1:Na);
% kynew1=kynew(1:Nr); % backscattered field on spatial frequency
% matplotcorrect(kxnew1,kynew1,esnew.',40);% after polar reformatting
% xlabel('kx domain');ylabel('ky domain');
% colormap(1-gray); 
% colorbar
% title('ackscattered field on spatial frequency:after')
%%
bwkx=kxmax-kxmin;
bwky=kymax-kymin;
dxnew=pi/bwkx;% *****dx is decided by bwkx*****
dynew=pi/bwky;% range resolution 
xnew=(-Na/2:Na/2-1)*dxnew;
ynew=(-Nr/2:Nr/2-1)*dynew;
xnew1=(-Na:Na-1)*dxnew;
ynew1=(-Nr:Nr-1)*dynew;

%% FFT processing
isar2=fftshift(ifft2(esnew));figure;
matplotcorrect(-xnew,-ynew,abs(isar2),40);
plotti(-xnew,-ynew,isar2,40);
% colormap(1-gray)
colorbar
title('imaging by polar reformatting');
% axis image
axis([-5 5 -3 3])
xlabel('range');ylabel('cross');