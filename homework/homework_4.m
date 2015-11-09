function homework_4(L,tau)
% This function is used to realize FIR
%  decimal times of sampling frequency 
%  L is the length of FIR, tau is decimal time delay
%  tau is of range -0.5---0.5
N=100;
Ts=0.001; % sampling interval
fs=1/Ts;    % sampling frequency
fd=linspace(0,0.5,N);
f=fd*fs;      % analog frequency
 l=L;
 % non-negative weighting function
l1=length(fd);      % number  of sampling in the passband 
l2=N-length(fd);
lamada=[ones(1,l1),zeros(1,l2)];
hd=[];
% determine the expected time delay
if l/2==1
    d=(l-1)/2;
elseif tau>0
    d=l/2-1;
else
    d=l/2+1;
end
% expected frequency response
hd=exp(-j*2*pi*fd*(d+tau));
% use L2 norm to solve impulse response
fre=[];
h=[];
for n=1:l1
    fre(n)=n*fs/l1;
    for m=1:l
        e(n,m)=exp(-j*2*pi*(m-1)*fre(n)/fs);        % matrix E
    end
end
h=pinv(e)*hd.';         % solve the filter out
hds=e*h;
% plot the amplitude response
figure;plot(fd,10*log10(abs(hd)),'o');title('amplitude');
hold on
plot(fd,10*log10(abs(hds)),'r');
hold off
% plot the phase response
figure;plot(fd,angle(hd),'o');title('phase');
hold on
plot(fd,angle(hds),'r');
hold off
% check the function of the filter
T=0.5;
f0=fs/5;
fl=f0/2;
fh=f0;
m1=500;
k=(fh-fl)/(2*T);
for m=1:m1
    wt2=2*pi*fl*(m/fs-tau)+pi*k*(m/fs-tau)^2;  
    st2(m)=sin(wt2);
 end
figure;
plot(st2,'b');
title('ideal delayed signal');axis([0 500 -1.4 1.4]);
% let the signal pass the filter
for m=1:m1
    wt2=2*pi*fl*(m/fs)+pi*k*(m/fs)^2;  
    st(m)=sin(wt2);        % original signal
end
 hc=[];
% for m=1:l
%     hc(m)=h(l+1-m);
% end
hc= fliplr(h);
s=[];
st1=[];
st1=[zeros(1,l-1),st];
for n=1:m1
    buff=[];
    buff=st1(n:n+14);
   s(n)=buff*hc;
end
s =2* real(s);
ss = zeros(1,m1);
ss(1:end-d) = s(d+1:end);
figure;
plot(-ss);
title('filtered signal');axis([0 500 -1.4 1.4]);
figure;
plot(-ss+st);title('error ');

% ff1=fft(st);
% ff2=fft(2*h,m1);
% ff=ff1.*ff2';
% s=ifft(ff);
% s=real(s);
% ss=zeros(1,m1);
% ss(1:end-d) = s(d+1:end);
% figure;
% plot(s);axis([0 500 -1.4 1.4]);
% figure;
% plot(s-st);



