%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m 
%
% Created May, 2024
% Elyes Balti
% The University of Texas at Austin
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the following papers: 
%
% E. Balti and B. K. Johnson, 
% "Tractable Approach to MmWaves Cellular Analysis With FSO Backhauling Under Feedback Delay and Hardware Limitations,"
% in IEEE Transactions on Wireless Communications, vol. 19, no. 1, pp. 410-422, Jan. 2020
% 
% E. Balti and M. Guizani, 
% "Mixed RF/FSO Cooperative Relaying Systems With Co-Channel Interference,"
% in IEEE Transactions on Communications, vol. 66, no. 9, pp. 4014-4027, Sept. 2018
%
% Contact email: ebalti@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script produces the analytical and Monte Carlo simulations of the Cumulative Distribution Function (CDF)
% of the Double Generalized Gamma distribution with the pathloss and the pointing errors.
%% Parameters 
% T: Threshold
% alpha: positive parameter related to the effective number of large-scale cells of the scattering process
% beta: natural number and it stands for the amount of fading parameter
% Mc: number of Monte Carlo iterations
% SNR: average SNR in linear scale
% ur: average electrical SNR
% SNRFSO: instantaneous SNR
% a1, a2, omega1, omega2: distribution parameters of the Double Generalized Gamma model
% p, q: positive parameters that satisfy p/q = a1/a2
% m1, m2: shaping parameters modeling the severity of fading
% X: small scale fluctuations
% Y: large scale fluctuations
% Im: optical irradiance (atmospheric turbulences)
% Cn: refractive index of the medium
% eta: electrical to optical conversion
% lambda: wavelength
% k: wave number
% L: link in meter
% rytov: Rytov variance
% a: radius of the receiver aperture
% sigma_s: jitter variance
% F0: radius of the curvature
% w0: beam waist at the transmitter
% Ip: pointing errors
% sigma: pathloss attenuation
% Il: pathloss
% r: detection metric (r=1: Homodyne, r=2: IM/DD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc
T = -20:1:10;
T = db2pow(T);
SNR_dB = 5;
SNR = db2pow(SNR_dB);
Mc = 1e4;
r = 2;
eta = 1;  


%% Atmospheric Turbulences
m1 = 3;
m2 = 2;
p = 2;
q = 1;
a1 = 4;
a2 = 2;

omega1 = (gamma(m1)/gamma(m1+1/a1))^a1 *m1;
omega2 = (gamma(m2)/gamma(m2+1/a2))^a2 *m2;

X = gamrnd(m1,omega1/m1,Mc,1);
Y = gamrnd(m2,omega2/m2,Mc,1);
X = X.^(1/a1);
Y = Y.^(1/a2);
Im = X .* Y;

%% Pointing Errors
Cn=sqrt(2*10^(-14)); 
L = 1000;
lambda=785*10^(-9);
k=2*pi/lambda;
rytov=sqrt(1.23*Cn^2*k^(7/6) * L^(11/6));
a=5*10^(-2);
sigma_s = .01;
F0 = -10;
w0=0.005;
Theta0=1-L/F0;
V0 = 2*L/(k*w0^2);
V1=V0/( Theta0^2 + V0^2 );
wL=w0 * sqrt( (Theta0 + V0) * (1 + 1.63*rytov^(12/5) * V1 ) );
v=sqrt(pi/2)*a/wL;
A0=abs(erf(v))^2;
wLeq = sqrt (wL^2 * sqrt(pi) * erf(v) / (2*v*exp(-v.^2)) );
xi=wLeq/(2*sigma_s);
R = raylrnd(sigma_s,Mc,1);
Ip = A0*exp(-2*abs(R).^2/wLeq^2);
h = xi^2/(xi^2+1);

%% Pathloss
sigma=5*10^(-4); 
Il=exp(-sigma*L);

%% Unified Double Generalized Gamma
I = Im.*Il.*Ip;


%% Initialize the CDF vectors
CDFa = zeros(length(T),1);
CDFm = zeros(length(T),1);

       if r==1
           ur = eta*SNR;
       else
           ur = h^2*A0^2*Il^2*eta^2*SNR^2;
       end

for ii=1:length(T)
%% Monte Carlo    
   for kk=1:Mc
       SNRFSO = I(kk)^r * ur;
      if SNRFSO < T(ii)
          CDFm(ii) = CDFm(ii) + 1;
      end      
   end
%% Analytical    
CDFa(ii) = UnifiedDoubleGeneralizedGammaCDF(m1,m2,p,q,a1,a2,Il,A0,xi,r,eta,SNR,T(ii));
end

CDFm = CDFm/Mc;% Averaging over the number of Monte Carlo iterations
T = pow2db(T);

figure; hold on
plot(T,CDFa,'linewidth',2)
plot(T,CDFm,'*','linewidth',2)
xlabel('Threshold (dB)')
ylabel('Cumulative Distribution Function (CDF)')
title('Unified Double Generalized Gamma Distribution')
legend('Analytical','Monte Carlo','location','northwest')
