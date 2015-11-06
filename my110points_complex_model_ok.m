% created on 2015.10.29
% this is written to test all kinds of conditions,including different
% parameter settings ,velocity ,prf and so on to find the reasons of
% azimuth sidelobes.

% based on paper:
% data level fusion of multilook ISAR images 啊，也不一定是哪篇文章了
%% modified on 2015.11.03 
%% use dongxiao 's range alignment method ,and it does not have
%% multi-images in azimuth. (It seems that a polyfit is neccessary 
%% for range correction)
%% 人为设定几个强散射点，运动模型为老师要求的

%% 注意1:fft 和 ifft 要一致 
%% 2: 用的chirp 信号
%% 总结 ： chirp 信号 , 匹配滤波，包络对齐（简单包络对齐方法）
%% 初相校正 然后对ESS粗包络对齐
%% 得到ESSnew并存储，装载ESSnew，将它分成8个子序列，分别包络对齐，
%% 初相校正，然后方位向加窗去除重复isar图像，得到isarwind1，提取
%% attention：此程序dy代表距离向，dx代表方位向
% 关于circshift 对于列向量：aa1 = circshift(aa,2) ，表示从末端数前移到开头
% 而aa1 = circshift(aa,-2) 表示开头数字移动到末端
% aa = 1     2     3     4     5     6     7
% aa1 =
% 
%      3
%      4
%      5
%      6
%      7
%      1
%      2
clear all;
close all;

r0  = 4e3; % initial radial range 
xsize = 80; %距离向extend
ysize = 80;

c  = 3e8;
fc = 10e9;% carrier frequency
lamada = c / fc;% wave lengths

b  = 3e8;                                                                                                           
fs  = 6e8;% sampling frequency
% M = floor(prfmax / prfmin);
% M0 = 2 * floor( 0.055 / 0.02 * prf); % number of pulses;
M   = 512; 
M1  = 512 * 8;
tp  = 2e-6;
prf = .45e3; % small prf
% prf = 1e3;% big prf
ts  = 1 / fs;
pri = 1 / prf;% pulse interval
N = 2 * floor((tp+2*xsize/c) / 2 / ts);% N sample points of each pulse

load Fighter3 
%%
figure;
plot(Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);grid;
set(gca,'FoNtNaMe', 'Arial', 'FoNtSize',12,'FoNtWeight', 'Bold'); 
axis([-20 20 -20 20])
xlabel('X [M]'); ylabel('Y [M]');

% theta0 = 0 * pi / 180; % initial angle of the target for radar1
% [theta,ra] = cart2pol(Xc,Yc);
% theta = theta + theta0; % add initial angle

dy = c / 2 / fs;% range resolution
y = -dy * N/2:dy:dy * (N/2-1);% range vector

dt = ts;% sampling time interval of ONE PULSE
tf = 2 * r0 / c + (-N/2:N/2-1) * dt; % time vector along chirp pulse:fast time

% f = fc+b/(N-1)*(-N/2:N/2-1);  % frequency vector
fr = fc+ 2*b/(N-1)*(-N/2:N/2-1); % 1103
% tm = zeros(M1,N);
% for p = 1:M1
%      tm(p,1:N) = pri * (p-1) + tf ; 
% end
% ts = 2*r0/c + tp/2:pri:2*r0/c + tp/2 + (M-1)*pri;% slow time
tslow = -N/2*dt:pri:-N/2*dt+(M-1)*pri;% slow time
vtx = 250 ;% from 200 to 150 :azimuth velocity 
vty = 15; % from 30 to 15 : translational velocity
% omega0 = vty / r0;
% omega = (vty+10 * tslow) / r0;

% bdop = 2 * omega * ysize / lamada; % 多普勒带宽
% prfmin = bdop;
prfmax = c / 2 / r0;

kchrp = b / tp;
tau = 2 * r0 / c;
downfreq = exp(-j*2*pi*(fc*(tf-tau)));
s  = exp(j*2*pi*(fc*tf+kchrp/2*(tf.^2)));% original signal
sr = exp(j*2*pi*(fc*(tf-tau)+kchrp/2*((tf-tau).^2)));% replica
H  = conj(fft(sr.*downfreq)/N);   %下变频到基带的matched filter transfer function
F  = (-N/2:N/2-1)*fs/(N-1);
%% 将目标模型的坐标转换到成像坐标下
theta = atan(vty/vtx);
for l = 1 : length(Xc)
    Xcnew(l) = Xc(l) * cos(theta) - Yc(l) * sin(theta);
    Ycnew(l) = Xc(l) * sin(theta) + Yc(l) * cos(theta);
end
figure;
plot(Xcnew,Ycnew,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);grid;
set(gca,'FoNtNaMe', 'Arial', 'FoNtSize',12,'FoNtWeight', 'Bold'); 
axis([-20 20 -20 20])
xlabel('X [M]'); ylabel('Y [M]');
%% 回波形成
for m = 1: M
    Es1(m,1:N) = zeros(1,N);
    for l = 1: length(Xcnew);
        range = sqrt((r0+vty*tslow(m)+Ycnew(l))^2+(vtx*tslow(m)+Xcnew(l))^2);
        phase = 2*pi*(fc*(tf-2 * range / c) + kchrp/2*((tf-2*range/c).^2));
        if mod(l,19) == 0
             Es1(m,1:N) = Es1(m,1:N) + 12* exp(j * phase); 
        else
             Es1(m,1:N) = Es1(m,1:N) + exp(j * phase);  
        end   
    end
    EsF(m,1:N) = fft(Es1(m,1:N).*downfreq)/N;%匹配滤波后再去斜
    ESS(m,1:N) = EsF(m,1:N).*H;%.*fft(exp(-j*2*pi*kchrp/2*((tf-tau).^2)));
    ESS(m,1:N) = ifft(ESS(m,1:N));
end
% save('my110points_complex_model_smlprf_ESS.mat','ESS');% 未包络对齐的距离像
% save('my110points_complex_model_smlprf_ESS_1103_250_15.mat','ESS');% 未包络对齐的距离像

% load my110points_complex_model_smlprf_ESS_1103_250_15.mat;

ESS1 = fftshift(ESS,2);
plotti(1:N,1:M,ESS1,60);title('匹配滤波后的距离像');
plotti(1:N,1:M,fftshift(fft(ESS1),1),60);title('isar image after matched filtering');

%% test RangeAlignment.m and phasealignment written by V.C Chen 
%% 为什么用我的匹配滤波后的数据调用VC chen的方法，成像结果不太对
%% 但是 包络对齐后调用可以

% testRng = RangeAlignment_chen_modify(ESS1*N*N);
% plotti(1:N,1:M,testRng,60);title('range profile after vc_chen rangealign method ');% rangealignment profile
% testISA = fft(testRng);
% % testISA = fft(ESS1);
% plotti(1:512,1:M,testISA,60);title('isar image after vc_chen rangealign method ');% rangealignment isar image
% [chen_rgnpro,phi,ref] = PhaseAdjustment(testRng,4);
% isar_chen = fftshift(fft(chen_rgnpro),1);
% plotti(1:N,1:M,isar_chen,60);title('isar像');

%% small prf data collection
ESSnew = RangeAlignment_cc_center(ESS1, 1, fr, -dt, tslow, 2);  % dongxiao ‘ S range alignment method

% save('my110points_complex_model_smlprf_ESSnew1101_250_15_my110ps.mat','ESSnew');
% save('my110points_complex_model_smlprf_ESSnew0907_200_30.mat','ESSnew'); 

% load my110points_complex_model_smlprf_ESSnew1101_250_15_my110ps.mat;

plotti(1:M,1:N,ESSnew,60);title('range profile after range alignment: center and polyfit');
plotti(1:M,1:N,fftshift(fft(ESSnew),1),60);title('isar image after range alignment: center and polyfit');

% find phase function after range compression at rangecell number
% ,0907_200_0: prf=450,m=512;
% figure;plot(1:M,unwrap(angle(ESSnew(709,:))));

%% apply PhaseAdjustment(X,nref) of VC Chen
[chen_rgnpro,phi,ref] = PhaseAdjustment(ESSnew,5);
% chen_rgnpro = pga_autofocus(ESSnew, 50, 7); % dongxiao's auto-focus program. for big velocity azimuth
isar_chen = fftshift(fft(chen_rgnpro),1); 
win = hanning(M) * ones(1,N); %  window
addwin0 = fftshift(fft(chen_rgnpro.*win),1);
plotti(1:N,1:M,isar_chen,50);title('isar image after phase correction(vc chen method)');
plotti(1:N,1:M,addwin0,50);title('isar image after phase correction(vc chen method)');
% saveas(gcf,'E:\Nutstore\result1\20151101_vc_chen\isarofESSnew1101_250_15_my110ps.png');
%% apply doppler centroid : jiefang yang
S_phasecomed = centriodtrack(ESSnew);
isar_DopCentroid = fftshift(fft(S_phasecomed),1); 
plotti(1:N,1:M,isar_DopCentroid,50);
title('isar image after phase correction(doppler centroid) -- not ok');
%% windowing before keystone
% fftshift(fft2((ifft(chen_rgnpro.')).'.*win),1);
win = hanning(M) * ones(1,N); %  window
addwin = fftshift(fft(chen_rgnpro.*win),1);
plotti(1:N,1:M,isar_chen,50);title('isar before adding window , no keystone');
plotti(1:N,1:M,addwin,50);title(' isar after adding window,no keystone');
%% apply keystone to range aligned range profile: ESSnew
data = ifft(ESSnew.');  % get original backscattered data
% data = ifft(chen_rgnpro.');
delta_f = 2*b / N;% 2 * b : attention here
X_temp = zeros(N,M);
% apply Keystone transform to original backscattered data 
for n = 1 : N
     for m = 1 : M
        m_old = (m - M / 2) * fc / (fc + (n - N / 2) * delta_f) + M / 2;
        if m_old >= 5 & m_old <= (M-5)
            for k = round(m_old)-4 : round(m_old) + 4
                if abs(m_old-k) < 0.01
                    X_temp(n,m) = X_temp(n,m) + data(n,k);
                else
                    X_temp(n,m) = X_temp(n,m) + ...
                              data(n,k) * sin(pi*(m_old-k))/(pi*(m_old-k));
                end
            end
        end
    end
end

%% after keystone of vc chen  , phase correction
rangeprofile = fft(X_temp);
[chen_rgnpro1,phi1,ref1] = PhaseAdjustment(rangeprofile.',3); % apply phase correction

plotti(1:N,1:M,ESSnew,50);title('before keystone: range profile')
plotti(1:N,1:M,chen_rgnpro1,50);title('after keystone: range profile'); % it sames that the two are similar

win_key = hanning(N) * ones(1,M);  % hanning window
isar_chen1 = fftshift(fft(chen_rgnpro1),1); 
addwin1 = fftshift(fft(chen_rgnpro1.*win),1);

plotti(1:N,1:M,isar_chen1,50);title('before win:isar image after keystone , phase correction : vc_chen');
plotti(1:N,1:M,addwin1,50);title('after win:isar after keystone , phase correction : vc_chen');
plotti(1:N,1:M,addwin,50);title(' isar after adding window,no keystone');
%% It seems that vc_chen 's keystone is not appropriate for this simulation
%% using keystone is worse than before

%% dong xiao 's test experiment
rngpro_ks = ifty(Keystone(fty(ESSnew).', fc, delta_f).');
chen_rgnpro2 = pga_autofocus(rngpro_ks, 50, 7); % for big  azimuth velocity : by Dong xiao
[chen_rgnpro22,phi22,ref22] = PhaseAdjustment(rngpro_ks,5); % apply phase correction

plotti(1:N,1:M,ESSnew,50);title('before keystone: range profile')
plotti(1:N,1:M,rngpro_ks,50);title('after keystone: range profile'); % it did straight the range line

isar_chen2 = fftshift(fft(chen_rgnpro2),1); 
isar_chen22 = fftshift(fft(chen_rgnpro22),1); 
addwin2 = fftshift(fft(chen_rgnpro2.*win),1);
addwin22 = fftshift(fft(chen_rgnpro22.*win),1);

plotti(1:N,1:M,isar_chen2,50);title('before win isar after keystone , phase correction : PGA');
plotti(1:N,1:M,addwin2,50);title('windowed isar after keystone , phase correction : PGA');
plotti(1:N,1:M,addwin22,50);title('windowed isar after keystone , phase correction : vcchen');
plotti(1:N,1:M,addwin,50);title(' windowed isar without keystone');
%% It seems that pga method is better than PhaseAdjustment of vc chen 
%% and the keystone achieve the goal of focusing

%% phsae function observation
figure;
plot(unwrap(angle(ESSnew1(:,ref(3)))),'-b','LineWidth',2)
axis tight
ylabel('Phase function')
xlabel('Pulses')
title(sprintf('Phase function at range cell %g before phase adjustment',...
               ref(3)));
drawnow

figure
plot(unwrap(angle(chen_rgnpro(:,ref(3)))),'-b','LineWidth',2)
axis tight
ylabel('Phase function')
xlabel('Pulses')
title(sprintf('Phase function at range cell %g after phase adjustment',...
               ref(3)));
drawnow

%ISAR image after phase adjustment
G1 = fftshift(fft(chen_rgnpro),1);

figure
colormap(jet)
imagesc(1:N,1:M,20*log10(abs(G1)));
xlabel('Range (meter)')
ylabel('Doppler (Hz)')
title('ISAR Image after phase adjustment');
axis('xy')
clim = get(gca,'CLim');
set(gca,'CLim',clim(2) + [-50 0]);
drawnow

% saveas(gcf,'E:\Nutstore\result1\20150914\isarofESSnew0907_200_30_450prf1024m.png');
% prf450 vtx 200 vty 30 azimuth numbers 1024

