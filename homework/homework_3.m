function homework_3(M,phi,theta0)
% more detail in array signal processing by Shefeng Yan
% This function is used to realize wide band beamforming with
% DFT in the frequency domain
% M is the number of array element 
% Theta0 is the ideal direction of signal
% total time series length of signal source
% input is simple
N=512;      % sample rate :5 times of original periodicity
T=0.5;
t=linspace(0,T,N*5);   % time range ,unit micro-second
f0=204.8;      % central frequency of chirp signal
fl=f0/2;        % lower bound of chirp
fh=f0; 
fs=5*f0;
k=(fh-fl)/(2*T);
% plot the original signal
 for m=1:length(t)
    wt=2*pi*fl*t(m)+pi*k*(t(m))^2;  
    st(m)=sin(wt);
 end
figure;
plot(t*1000,st,'b');
title('original chirp signal');axis([0 500 -1.4 1.4]);
%plot the samplint signal
m1=length(t)/5;
 for m=1:m1
    wt2=2*pi*fl*m/fs+pi*k*(m/fs)^2;  
    st2(m)=sin(wt2);
 end
 figure;
plot(st2,'b');
title('sampled chirp signal');axis([0 500 -1.4 1.4]);

theta=-pi/2:pi/511:pi/2;    % x label;
vtheta=-[sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)*ones(size(theta))];
y=[];   
c=3.0*1e8;   % speed of light
% position of array
for m=-M/2+1:M/2
        y=[y, [0,m,0]'];
end

% plot the psd function of signal    (indirect way)
nfft=512;
re=xcorr(st,'unbiased');
ref=fft(re,nfft);
pr=abs(ref);
index=0:round(nfft/2-1);
fs=5*f0;
k1=index*fs/nfft;
f1=20*log10(pr(index+1));
                
% plot the psd function of signal    (direct way)
window=boxcar(length(st)); % rectangular window
[Pxx,f]=periodogram(st,window,nfft,fs);
figure;
plot(f/100,10*log10(Pxx));title('psd of signal st');axis([0 2.5 -90 -10]);

tau=[];
tau=vtheta'*y/c;
% plot the signal array element received
wt1=[];
st1=[];
figure;
for n=1:M   
 for m=1:m1
    wt1(m)=2*pi*fl*(m/fs-tau(1,n))+pi*k*(m/fs-tau(1,n))^2;  
    st1(n,m)=sin(wt1(m));
 end
 plot(st1(n,:)+n*2.2);
 hold on;
end
hold off

s1=st1(:,1:256);   % turn the data into two parts without overlap
s2=st1(:,257:512);
% FFT of the data
fs1=[];
fs2=[];
wc=[];
for n=1:M
    fs1(n,:)=fft(s1(n,:));
    fs2(n,:)=fft(s2(n,:));
end
fina=zeros(512,1);       
for n=21:56
    omega=n*fs/256;  % frequency of  sub-band
    for m=1:M
        ak(m)=exp(-j*pi*n*omega/f0*sin(theta0));   % atheta0
    end
    wc(:,n)=ak/M;    % weight vector  regular beamforming
    fina(n)=fs1(:,n)'*wc(:,n);
end
ift1=ifft(fina);
    
for n=21+256:56+256
    omega=n*fs/256;  % frequency of  sub-band
    for m=1:M
        ak(m)=exp(-j*pi*n*omega/f0*sin(theta0));   % atheta0
    end
    wc(:,n)=ak/M;    % weight vector
    fina(n)=fs2(:,n-256)'*wc(:,n-256);
end
ift=ifft(fina);
ift=real(ift);
figure;
plot(2*ift);axis([0 500 -1 1]);
  

