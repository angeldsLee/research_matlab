function homework_5(M,theta1,theta2)
% This function is to realize Iss-music which
% has two signal from theta1 and theta2 
% Input M is the number of array element

% copy from others ,i do not understand it totally
f0=204.8;
fl=f0/2; fh=f0;     % work frequency
fs=3.2*f0;          % sampling frequency
N=256;
len=N*50;          % number of sample data 
T=len/fs;            % time duration
t=(1:len)/fs;
theta=-90:0.1:90;
sx=zeros(41,M,M);% cross-spectral matrix
sxba=sx;
for n=1:50
    for m=1:M;
        tau1=1/2*sin(theta1*pi/180)/fh*(m-1);
        tau2=1/2*sin(theta2*pi/180)/fh*(m-1);
        for k=1:N
        s1(k)=sin(2*pi*(fl +(fh-fl)*(t(k+(n-1)*N)-tau1)/2/T).*(t(k+(n-1)*N)-tau1));
        s2(k)=sin(2*pi*(fl +(fh-fl)*(t(k+(n-1)*N)-tau2)/2/T).*(t(k+(n-1)*N)-tau2));
        end
    x(m,:)=s1+s2;
    X(m,:)=fft(x(m,:),N);
    end
    for k=40:80
    sx(k-39,:,:)=X(:,k)*X(:,k)';
    end
    sxba=sx+sxba;
end
en=zeros(M,10);    % noise subspace
az=0;
sxba1=zeros(M,M);
sxba=sxba/50;
for k=40:80
      q=1;
      for n=1:M
      sxba1(n,:)=sxba(k-39,n,:);
      end
      [v,d]=eig(sxba1); % eigen value decomposition 
      for n=1:10
      en(:,n)=v(:,n);
      end
      for n=-90:0.1:90;
          for m=1:M;
                if k<=N/2
                a(m)=exp(j*pi*fs/fh*k/N*sin(n.*pi/180)*(m-1));
                else
                a(m)=exp(j*pi*fs/fh*(k-N)/N*sin(n.*pi/180)*(m-1));
                end
           end
       music(q)=1/(a*en*en'*a');        % MUSIC azimuth beamshape
       q=q+1;
    end 
    az=music+az;
end
az=10*log10(abs(az/41)/max(abs(az/41)));
figure;
plot(theta,az);
xlabel('azimuth');
ylabel('azimuth shape/dB');
axis([-40,40,-40,0])
 
