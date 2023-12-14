function [U0,Ux] = TimoshenkoBeam(frq,xx,x0,L,A,I,E,nu,rho,eta,kap,BcL,BcR,kb,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TimoshenkoBeam calculates the response of a finite length Timoshenko beam
% (bending (u) & rotation (th) DoF) for different combinations of boundary 
% conditions at the ends for a unit force and unit torque and isadapted 
% from the book of DJT [1] with extension of different BC.
%
%   The function inputs are:
%   frq: frequency in [Hz]
%   xx: response position(s)
%   x0: excitation position in [m]
%   L: length of the beam in [m]
%   A: cross-section in [m2]
%   I: second moment of the area in [m]
%   E: Young's modulus of the beam in [Pa]
%   nu: Poisson ratio of the beam in [-]
%   rho: density ratio of the beam in [kg/m3]
%   eta: loss factor of the beam  in [-]
%   kap: shear correction factor in [-]
%   BcL: boundary condition of the left beam end ('S','F','C','R', i.e. simply supported/free/clamped/roller) 
%   BcR: boundary condition of the right beam end ('S','F','C','R')
%   kb: dynamic stiffness of elastic foundation in [N/m]
%   dt: differentiation w.r.t time (0:displacement,1:velocity,2:acceleration)
%
%   The function inputs are:
%   U0: driving point response (x=x0)
%   Ux: transfer response (at all defined x)
%
%   created by Christopher Knuth (C) during PhD in ISVR
%   at University of Southampton in 2020-2024 (ck1g19@soton.ac.uk)
%   [1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% materaial and geometry

%complex Young's modulus
E = E*(1+1i*eta);
%shear modulus
G = E/(2*(1+nu));

%distance to sleeper end from excitation position x0=0
L1 = -x0;
L2 = L-x0;

%circular frequency
w = 2*pi*frq;

%% sleeper parameters

%shear stiffness
GAk = G*A*kap;
%bending stiffness
EI = E*I;
%mass
rA = rho*A;
%inertia
rI = rho*I;

%% ballast parameters

%translation (u)
kb = kb/L;  
%rotation (th)
kbr = 0;

%% calculate wavenumbers

% set up mass and stiffness matrices (vertical bending + rotation + axial compression)
M = [rA,0;0,rI];
Kb =[kb,0;0,kbr];
K0 = [0,0;0,GAk];
K1 = [0,GAk;-GAk,0];
K2 = [-GAk,0;0,-EI];

%rewrite 2nd order QEVP to 1st order EVP
O = zeros(size(M));
I = eye(size(M));
A1 = [K0-w^2*M+Kb,-K1;O,I];
A2 = [O,K2;-I,O];

% calculate wavenumbers
[~,lambda,~] = eig(A1,-1i*A2);
kn = diag(lambda);
N = length(kn);
%relate amplitudes via free vibration
D11 = zeros(N,1);D12 = zeros(N,1);
D21 = zeros(N,1);D22 = zeros(N,1);
for ii=1:N
    DSM(:,:) = K0+Kb-w^2*M-1i*kn(ii)*K1-kn(ii)^2*K2;
    D11(ii) = DSM(1,1);
    D12(ii) = DSM(1,2);
    D21(ii) = DSM(2,1);
    D22(ii) = DSM(2,2);
end

%% calculate wave amplitudes An and Bn

%set up coefficient matrix to calculate wave amplitudes for given boundary
%conditions of the beam
AAA = zeros(16,16);

%boundary condition on left end of the beam x=-x0
if strcmp(BcL,'F') %free
    %shear force equals zero
    AAA(1,0*N+1:1*N) = (-1i*kn).*exp(-1i*kn*L1);
    AAA(1,2*N+1:3*N) = -exp(-1i*kn*L1);
    %bending moment equals zero
    AAA(2,2*N+1:3*N) = (-1i*kn).*exp(-1i*kn*L1);
elseif strcmp(BcL,'C') %clamped
    %displacement equals zero
    AAA(1,0*N+1:1*N) = exp(-1i*kn*L1);
    %rotation equals zero
    AAA(2,2*N+1:3*N) = exp(-1i*kn*L1);
elseif strcmp(BcL,'S') %simply supported
    %displacement equals zero
    AAA(1,0*N+1:1*N) = exp(-1i*kn*L1);
    %bending moment equals zero
    AAA(2,2*N+1:3*N) = (-1i*kn).*exp(-1i*kn*L1);
elseif strcmp(BcL,'R') %roller (or slider) support
    %shear force equals zero
    AAA(1,0*N+1:1*N) = (-1i*kn).*exp(-1i*kn*L1);
    AAA(1,2*N+1:3*N) = -exp(-1i*kn*L1);
    %rotation equals zero
    AAA(2,2*N+1:3*N) = exp(-1i*kn*L1);
end

%boundary condition on right side of the beam
if strcmp(BcR,'F') %free
    %shear force equals zero
    AAA(3,1*N+1:2*N) = (-1i*kn).*exp(-1i*kn*L2);
    AAA(3,3*N+1:4*N) = -exp(-1i*kn*L2);
    %bending moment equals zero
    AAA(4,3*N+1:4*N) = (-1i*kn).*exp(-1i*kn*L2);
elseif strcmp(BcR,'C') %clamped
    %displacement equals zero
    AAA(3,1*N+1:2*N) = exp(-1i*kn*L2);
    %rotation equals zero
    AAA(4,3*N+1:4*N) = exp(-1i*kn*L2);
elseif strcmp(BcR,'S') %simply supported
    %displacement equals zero
    AAA(3,1*N+1:2*N) = exp(-1i*kn*L2);
    %bending moment equals zero
    AAA(4,3*N+1:4*N) = (-1i*kn).*exp(-1i*kn*L2);
elseif strcmp(BcR,'R') %roller (or slider) support
    %shear force equals zero
    AAA(3,1*N+1:2*N) = (-1i*kn).*exp(-1i*kn*L2);
    AAA(3,3*N+1:4*N) = -exp(-1i*kn*L2);
    %rotation equals zero
    AAA(4,3*N+1:4*N) = exp(-1i*kn*L2);
end

%continuity of displacement at x=x0 (excitation position)
AAA(5,0*N+1:1*N) =  1;
AAA(5,1*N+1:2*N) =  -1;
%continuity of rotation at x=x0 (excitation position)
AAA(6,2*N+1:3*N) = 1;
AAA(6,3*N+1:4*N) = -1;

%equilibrium of shear force at x=x0 (excitation position)
AAA(7,0*N+1:1*N) = (-1i*kn);
AAA(7,1*N+1:2*N) = -(-1i*kn);
AAA(7,2*N+1:3*N) = -1;
AAA(7,3*N+1:4*N) = 1;
%equilibrium of bending moment at x=x0 (excitation position)
AAA(8,2*N+1:3*N) = (-1i*kn);
AAA(8,3*N+1:4*N) = -(-1i*kn);

%relate wave amplitudes An with Bn via free vibration equation
AAA(2*N+1:3*N,0*N+1:1*N) = diag(D11);
AAA(2*N+1:3*N,2*N+1:3*N) = diag(D12);
AAA(3*N+1:4*N,1*N+1:2*N) = diag(D21);
AAA(3*N+1:4*N,3*N+1:4*N) = diag(D22);

%set up force vector to calculate wave amplitudes
F = zeros(4*N,1); F(2*N-1) = 1/GAk;
M = zeros(4*N,1); M(2*N) = 1/EI;

%re-scale (may be necessary for either very big or very small numbers?)
fac = 1;
F = F./fac;
M = M./fac;
AAA = AAA./fac;

%calculate wave amplitudes
AmpF = AAA\F;
AmpM = AAA\M;

An(:,1) = AmpF(1:2*N);
An(:,2) = AmpM(1:2*N);
Bn(:,1) = AmpF(2*N+1:4*N);
Bn(:,2) = AmpM(2*N+1:4*N);

%% calculate response (driving point and transfer)

%driving point

%calculate driving point response (u,th) per unit force F
u0F  = (1i*w)^dt*sum(An(1:N,1));
th0F = (1i*w)^dt*sum(Bn(1:N,1));
%calculate driving point response (u,th) per unit torque M
u0M  = (1i*w)^dt*sum(An(1:N,2));
th0M = (1i*w)^dt*sum(Bn(1:N,2));

%combine in matrix
U0 = [u0F,u0M;th0F,th0M];

%transfer
x = xx-x0;
Ux = zeros(2,2,length(x));
%divide beam in two sections left and right to excitation
x1 = x(x<=0&x>=L1);
x2 = x(x>0&x<=L2);
%calculate transfer response per unit force and unit torque
if ~isempty(x1) && ~isempty(x2) %both values x1 and x2 nonzero
    uF  = (1i*w)^dt*[sum(An(1:N,1).*exp(-1i*kn.*x1),1),sum(An(N+1:2*N,1).*exp(-1i*kn.*x2),1)];
    thF = (1i*w)^dt*[sum(Bn(1:N,1).*exp(-1i*kn.*x1),1),sum(Bn(N+1:2*N,1).*exp(-1i*kn.*x2),1)];
    uM  = (1i*w)^dt*[sum(An(1:N,2).*exp(-1i*kn.*x1),1),sum(An(N+1:2*N,2).*exp(-1i*kn.*x2),1)];
    thM = (1i*w)^dt*[sum(Bn(1:N,2).*exp(-1i*kn.*x1),1),sum(Bn(N+1:2*N,2).*exp(-1i*kn.*x2),1)];
elseif isempty(x2) %only x1 exists
    uF  = (1i*w)^dt*[sum(An(1:N,1).*exp(-1i*kn.*x1),1)];
    thF = (1i*w)^dt*[sum(Bn(1:N,1).*exp(-1i*kn.*x1),1)];
    uM  = (1i*w)^dt*[sum(An(1:N,2).*exp(-1i*kn.*x1),1)];
    thM = (1i*w)^dt*[sum(Bn(1:N,2).*exp(-1i*kn.*x1),1)];
elseif isempty(x1) %only x2 exists
    uF  = (1i*w)^dt*[sum(An(N+1:2*N,1).*exp(-1i*kn.*x2),1)];
    thF = (1i*w)^dt*[sum(Bn(N+1:2*N,1).*exp(-1i*kn.*x2),1)];
    uM  = (1i*w)^dt*[sum(An(N+1:2*N,2).*exp(-1i*kn.*x2),1)];
    thM = (1i*w)^dt*[sum(Bn(N+1:2*N,2).*exp(-1i*kn.*x2),1)];
end

%combine response in matrix
for ix = 1:length(x)
    Ux(:,:,ix) = [uF(ix),uM(ix);thF(ix),thM(ix)];
end


end