clear all;
figure(300);clf;
figure(50);clf;
%-------FIR---
% N=3;
% Fs=Fc;% data sampled 
% Fp=200e6;% passband-edge frequency 
% Ap=3; %passband ripple
% Ast=60; %stopband attenuation 
% Rp = (10^(Ap/20) - 1)/(10^(Ap/20) + 1); 
% Rst = 10^(-Ast/20);
% NUM = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge');
% % NUM =floor(10*NUM)/10,
NUM =[1,1];
%----INL DNL----
a=0;
for i=1:63;
   rr(i)=a+1;
   a=rr(i);
end,
RR6bit_ideal(1:63)=rr(1:63)-32;
Fc=6.5e9 ;% RF carrier freq
% noise haping parameters of 3rd order errroe feedback
aa=-3; 
bb=-aa;
cc=-1;
DC=0;
AM=1; % amplitude
nbit=5 ; % 6bit nequest dac;
g1=(1-2^-32)/8; %1/8*(1-4/2^4); %unify input amplitude
fs=Fc*4; % sampling frequnecy
osr=8;% %over sampling ratio
ni=2000; % number of initial poit skiped
nn=2^12; % number of FFT point
n=nn;
RBW=fs/nn;%fft RBW
fbbw=249e6;
fin1=1*fbbw/8;
fin2=2*fbbw/8;
fin3=3*fbbw/8;
fin4=4*fbbw/8;
fin5=5*fbbw/8;
fin6=6*fbbw/8;
fin7=7*fbbw/8;
fin8=8*fbbw/8;

fin1=fs/nn*floor(fin1/fs*nn);
fin2=fs/nn*floor(fin2/fs*nn);
fin3=fs/nn*floor(fin3/fs*nn);
fin4=fs/nn*floor(fin4/fs*nn);
fin5=fs/nn*floor(fin5/fs*nn);
fin6=fs/nn*floor(fin6/fs*nn);
fin7=fs/nn*floor(fin7/fs*nn);
fin8=fs/nn*floor(fin8/fs*nn);


% rng=200;% set randon number seed 
for kk=1:100;
    rng(kk);
derinl2=1*(rand(1,63)-0.5); % DNL with +/- 1 LSB erroer
a=0;

a=0;
for i=1:63;
   rr2(i)=a+derinl2(i)+1;
   a=rr2(i);
end,

n0=1;
n1=63;
p=polyfit(n0:n1,rr2(n0:n1),1);
INLslop=polyval(p,n0:n1);
INL=rr2(n0:n1)-INLslop; %INL
RR6bit(1:63)=rr2(1:63)-32;
rr2(64)=rr2(63)+1;
DNL=rr2(n0+1:n1+1)-rr2(n0:n1)-1;

figure(300);
% clf;
plot(INL,'b-'); hold on;
plot(DNL,'r-'),hold on;
grid on;
grid on;

%-end of ---INL DNL----




% fvtool(NUM,'Fs',Fs);
% ----parameter for simulink model----------

sim('DPA_6bitDAC_MIXER_INL');

figure(50),
% clf,
[p_INL,f] = pwelch(y_INL,nn,nn/8,4*nn,fs); %6 bit with DC SD 
[p_ideal,f] = pwelch(y_ideal,nn, nn/8,4*nn,fs); % perfect 6bit
fmask=[-85,-85, -70 -70,-41.3, -41.3, -65,-65,-85,-85];
fmm=[1.6e9,3.79e9, 3.8e9, 5.99e9, 6e9,8.49e9, 8.5e9, 10.59e9,10.6e9,16e9];
plot(fmm, fmask , 'g-');grid on;hold on,
dbp_INL=10*log10(p_INL)-max(10*log10(p_ideal))-41.3;%unify to band psd 
dbp_ideal=10*log10(p_ideal)-max(10*log10(p_ideal))-41.3;%unify to band psd
maxy=-max(10*log10(p_ideal))-41.3;% calibration 
plot(f,dbp_ideal,'r.-'); hold on;
plot(f,dbp_INL,'b-');hold on;
end
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD (dBm/Mhz)'),

% YDC=mean(y_ideal)*31;
% % % end of the .m file  
