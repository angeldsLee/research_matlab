% created on 2014.12
%% SCFW пе╨е
% this is a simple translatioNal MotioN coMpeNsatioN program
% which compensates only the TRANSLATION MotioN error
% term in the phase representation ,More detail iN pp301(isar Matlab)
%%  RANGE TRACKING : cross-correlation Method  , example. 8.3.1.1pp308
%% atteNtioN 1 : when parameters are changed slightly:
% M=144,THEN it needs circshift(isarcomp(:,(1:M))',-8)';
% if M=128 ,then it needs no circshift
%% attention 2 :  SFCM signal !!!!!
%% question: why the target position is leaned to up_right???
% because the f vector should be f=f0:b/(N-1):(f0+b); !!!!!!!!
%%  basic paraMeters defiNitioN
close all
clear all
theta0 = 0*pi/180; % iNitial aNgle of target !!!!!!!!!
r0      =  16e3;
vt      =  70.0;
at     =  0.1;
phir  = 0.03; % uNit:rad per secoNd
M    =  144; % NuMber of burst
N    =  128; % NuMber of pulses
f0   =  10e9;% the first pulse frequeNcy:10GHz
b   =  128e6;% frequeNcy baNdwidth
prf  =  20e3; % pulse repetitioN frequeNcy
t1   = 1 / prf; % pulse repetition interval
c = 3.0e8;
%%  load data aNd plot ideal target Model
load Fighter 
plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);grid;
set(gca,'FoNtNaMe', 'Arial', 'FoNtSize',12,'FoNtWeight', 'Bold'); 
axis([-35 35 -30 30])
xlabel('X [M]'); ylabel('Y [M]');
%ScatteriNg ceNters iN cyliNdirical coordiNates
[theta,r] = cart2pol(Xc,Yc);
theta = theta + theta0 * 0.017455329; % add iNitial aNgle
x = Xc;y = Yc;
%% form the scattered field data
f = f0:b / (N-1):(f0+b);
k = 4 * pi * f / c;
kuse = ones(M,1)*k; 
t = zeros(M,N);
tdur=(N-1)/b;     % Pulse duratioN [s]
for p=1:M
    for q=1:N
        t(p,q)  = t1 * (q - 1 + N * ( p - 1)) + tdur / 2 + 2 * r0 / c;
        % Notice the above description of t vector
    end
end
%% reference phase representation (result to Matlab codes %%%%%%%
%%8.2_thru_8.6)%%%%%%
Es = zeros(M,N);
for l = 1:1:length(x)
    omega = - theta(l) + phir * t;
    range   =  r0 + vt * t + 0.5 * at * t.^2 - r(l) * sin(omega);
%     range=-r(l)*sin(omega);
    % still know few about the above formula
    % it may be polar format representation 
    phase  =  kuse .* range;
    Es  =  Es + exp(j * phase);
end
% save('E:\research\programs\data_fi.mat','Es');

%% use data produced from book (8.2 thru 8.6)
% load data_bookfi
Es = Es.';% n * m % !!!
%% isar iMagiNg usiNg 2D fft
% before range tracking
dx = c/2/b;
xx = -dx*(N/2-1):dx:dx*N/2;
% define cross range resolution
tot = M * t1 * N; % the coherent integration time
laMada = c / f0;
dy = laMada / 2 / tot / phir; % detail in pp242 equation (6.21)
yy = -dy * (M/2-1):dy:dy * M / 2;
ISAR = abs(fftshift(ifft2(Es)));
matplotcorrect(-yy,-xx,ISAR,20);%% or (-yy,1:N,ISAR,22);
plotti(-yy,-xx,ISAR,20);
colormap(1-gray); colorbar;
set(gca,'FontNaMe', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
ylabel('Range [M]'); 
xlabel('cross range');
%% translational Motion compensation
%%%%%% use cross correlation Method  %%%%%%%%
ccr=zeros(N,M);
rp=(ifft(Es)); %%%%%%     very important     n*m  %%%%
for u=1:M 
    tem=rp(:,u);% range profiles aNd choose rp(:,1) as the reference
    ccr(:,u)=abs(ifft(fft(abs(rp(:,1))).*conj(fft(abs(tem)))));% cross-correlatioN
end

%%%%%%%%% velocity estimation %%%%%%%%%
ind=zeros(1,M);
for dd=2:M
    ind(dd)=find(ccr(:,dd)==max(ccr(:,dd)));% iNdex of raNge shift(tiMe delays)
end
smo=smooth((0:M-1),ind,0.8,'rlowess');%sMoothiNg the delays 
figure;plot((1:M),ind,'o');% raNge shifts
hold on
plot((1:M),smo,'*');% sMoothed raNge shifts
rngshif=dx*smo; % raNge shift
rngdif=rngshif(2:M)-rngshif(1:M-1);
avr=mean(rngdif);
vesti=avr/t1/N;% estiMate the average speed
vest=rngdif/t1/N; % estiMate the speed betweeN bursts
figure;
plot((1:M-1),vest);
vest=[0;vest];% spand it into M bursts
%% form compensated scattered field data
es=Es.';% back to m*n
rpcomp=es.*exp(-j*kuse.*t*vesti);
Escomp=rpcomp.';
win = hanning(N)*hanning(M).'; %  window
isarcomp=abs((fft2(Escomp.*win)));%%%%why this need not fftshift?????????%%%%%
%%  IMaging after compensation
isar1=zeros(N,M);%%% if m or n becomes larger ,then it needs circshift %%%
isar1(:,(1:M))=circshift(isarcomp(:,(1:M))',-8)';
% ANOTHER way to shift isar data
isar2= isarcomp(:,54:M);   
isar2(:,M-54+2:M)=isarcomp(:,1:53);
xnew=(-N/2:N/2-1)*dx;
ynew=-dy*(M/2-1):dy:dy*M/2;
matplotcorrect1(ynew,-xnew,isarcomp,20);% n*m % !!!!
plotti(ynew,-xnew,isarcomp,40);
colormap(1-gray); colorbar;
set(gca,'FoNtNaMe', 'Arial', 'FoNtSize',12,'FoNtWeight', 'Bold'); 
xlabel('RaNge [M]'); 
ylabel('cross raNge');
