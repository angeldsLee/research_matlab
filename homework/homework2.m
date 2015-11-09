function homework2(M,phi,theta0,sigmax1,sigmax2,beta,sigman)
% This function is to plot beam azimuth shape of regular beamforming
% at azimuth 20 ,and at azimuth 0 and 30 we have two signals whose
% snr respectively 5dB and 10dB ,use regular beamforming to scan 
% form azimuth shape and plot
% its input is homework2(10,pi/2,20*pi/180,10^0.5,10,1,0)
theta=-pi/2:pi/200:pi/2;    % x label;
vtheta=-[sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)*ones(size(theta))];% coodinate of element
% for n=1:M
%     ax(n)=exp(-j*pi*n*sin(theta0));   % atheta0
%    %   ad(n)=exp(-j*pi*n*sin(thetad));   % athetad
% end
% 0 azimuth array manifield vector
for n=1:M
    ax1(n)=exp(-j*pi*n*sin(0));   % atheta0
   %   ad(n)=exp(-j*pi*n*sin(thetad));   % athetad
end
% 30 azimuth array manifield vector
for n=1:M
    ax2(n)=exp(-j*pi*n*sin(30*pi/180));   % atheta0
   %   ad(n)=exp(-j*pi*n*sin(thetad));   % athetad
end

Ex1=sigmax1*ax1'*ax1;  %input signal 5dB means 10^0.5
Ex2=sigmax2*ax2'*ax2;  %input signal 10dB means 10
%Ed=sigmad*ad*ad';    % unit is dB
Rx=beta*(Ex1+Ex2)+sigman;% +Ed+sigman;  % compute covariance matrix dB
out=beam(M,phi,theta0); % plot the beam response 
wtheta=out';
N=length(theta);
for n=1:N
    pth(n)=wtheta(:,n)'*Rx*wtheta(:,n);
end
ptheta1=10*log10(pth/M^2); 
figure;
plot(theta*180/pi,ptheta1);
xlabel('azimuth'); axis([-90 90 -30 40]);
ylabel('Azimuthshape/dB');
