%% Model of a finite length Timoshenko beam with constant cross-section 
%This calculates the FRF of a Timoshenko beam with free, simply supported,
%clamped or roller conditions on either of its two ends.

clc
clear all
close all

%% set parameters

%beam size
L = 2.5; %length 
w = 0.25; %width (either w and h or A and I directly!)
h = 0.2; %heigth

%cross-section parameters (can set any number for arbitrary cross-section)
A = w*h; %surface area 
I = w*h^3/12; %second moment of the area

%force location
x0 = 2.5; %make sure: 0 < x0 < L
%mesh of the beam
dx = 0.01; %element size
x = (0:dx:L);

%material
E = 4.3e10; %Young's modulus
nu = 0.15; %Poisson ratio
rho = 2500; %density
eta = 0.01; %loss factor
kap = 0.83; %shear correction factor

%ballast
kb = 0*200e6; %stiffness
dampb = 0*150e3; %damping coefficient (C) or damping loss factor (eta)
dampf = 2; %1: hysteretic damping [s*(1+i*eta)] 2: viscous damping [s+iw*C]

%frequency
fmin = 10; %minimum frequency
fmax = 10000; %maximum frequency
Nfrq = 400; %number of frequency points
frq = logspace(log10(fmin),log10(fmax),Nfrq);
% frq = linspace(fmin,fmax,Nf);
dt = 0; %differentiation w.r.t time (dt=0:displacement, dt=1:velocity, dt=2:acceleration)

%boundary condition left beam end
BcL = 'F'; %free
% BcL = 'S'; %simply supported
% BcL = 'C'; %clamped
% BcL = 'R'; %roller/slider

%boundary condition right beam end
BcR = 'F'; %free
% BcR = 'S'; %simply supported
% BcR = 'C'; %clamped
% BcR = 'R'; %roller/slider

%% calculate driving point and transfer response

tic
for ifrq = 1:Nfrq
    if dampf == 1 %hysteretic [s = s*(1+i*eta)]
        kb = kb*(1+1i*dampb);
    elseif dampf == 2 %viscous [s = s+iw*C]
        kb = kb+1i*w*dampb;
    end
    [U0(:,:,ifrq),Ux(:,:,:,ifrq)] = TimoshenkoBeam(frq(ifrq),x,x0,L,A,I,E,nu,rho,eta,kap,BcL,BcR,kb,dt);
end
toc

%% plot: driving point FRFs

figure(1)
clf
subplot(3,1,1:2)
loglog(frq,abs(squeeze(U0(1,1,:))),'b')
hold on
loglog(frq,abs(squeeze(U0(2,2,:))),'r--')
loglog(frq,abs(squeeze(U0(1,2,:))),'k--')
% loglog(frq,abs(squeeze(U0(2,1,:))),'g:')
hold off
box on
ylabel('Magnitude, m/N','Fontsize',14)
xlim([10 max(frq)])
legend('u/F','\theta/M','u/M','Location','nw','Fontsize',12)
subplot(3,1,3)
semilogx(frq,180/pi*angle(squeeze(U0(1,1,:))),'b')
hold on
semilogx(frq,180/pi*angle(squeeze(U0(2,2,:))),'r--')
semilogx(frq,180/pi*angle(squeeze(U0(1,2,:))),'g--')
hold off
box on
set(gca,'YTick',-180:90:180)
xlim([min(frq) max(frq)])
ylim([-1,1]*200)
ylabel('Phase, rad/s','Fontsize',14)
xlabel('Frequency, Hz')

%% plot: transfer FRF (at all frequencies)

for ii = 1:Nfrq
% for ii = 250
figure(2)
clf
hold on
plot(x-L/2,real(squeeze(Ux(1,1,:,ii))))
plot(x(x==x0)-L/2,real(squeeze(Ux(1,1,x==x0,ii))),'ko')
hold off
box on
xlim(1.1*[-L/2 L/2])
% ylim(pi*[-1,1])
xlabel('Frequency, Hz')
ylabel('Mobility, m/s/N')
title(num2str(frq(ii)))
pause(0.1)
end

%% plot transfer FRF at resonances (operational mode shape)

figure(3)
clf
[~,fp]= findpeaks(abs(squeeze(U0(1,1,:))),frq);
for ii = 1:9
ind(ii) = find(min(abs(frq-fp(ii)))==abs(frq-fp(ii)));
subplot(3,3,ii)
shape = squeeze(Ux(1,1,:,ind(ii)));
shapenorm = shape/squeeze(U0(1,1,ind(ii)));
plot(x,real(shapenorm),'k','Linewidth',2)
hold on
plot(x0,1,'ko')
title(['f = ', num2str(round(frq(ind(ii)))) ' Hz'])
% axis equal
xlim([min(x),max(x)])
ylim(ceil(max(abs(shapenorm)))*[-1,1])
end
