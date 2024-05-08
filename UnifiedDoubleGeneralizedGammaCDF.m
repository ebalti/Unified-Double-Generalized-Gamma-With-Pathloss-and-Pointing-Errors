%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnifiedDoubleGeneralizedGammaCDF.m 
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
% This script returns the analytical expression of the Cumulative
% Distribution Function (CDF) of the Double Generalized Gamma distribution
% with the pathloss and the pointing errors
%% Parameters 
% alpha: positive parameter related to the effective number of large-scale cells of the scattering process
% beta: natural number and it stands for the amount of fading parameter
% Mc: number of Monte Carlo iterations
% a1, a2, omega1, omega2: distribution parameters of the Double Generalized Gamma model
% p, q: positive parameters that satisfy p/q = a1/a2
% m1, m2: shaping parameters modeling the severity of fading
% Il: pathloss
% SNR: average SNR in linear scale
% ur: average electrical SNR
% r: detection metric (r=1: Homodyne, r=2: IM/DD)
% xi: pointing errors parameter
% x: threshold (argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = UnifiedDoubleGeneralizedGammaCDF(m1,m2,p,q,a1,a2,Il,A0,xi,r,eta,SNR,x)
                                              
h = xi^2/(xi^2+1);
if r ==1
ur = eta*SNR;
else
ur = h^2*A0^2*Il^2*eta^2*SNR^2;
end
        
omega1 = (gamma(m1)/gamma(m1+1/a1))^a1 *m1;
omega2 = (gamma(m2)/gamma(m2+1/a2))^a2 *m2;

A = xi^2 * p^(m2-3/2) * q^(m1-1/2) * (2*pi)^(1-(p+q)/2) / (gamma(m1)*gamma(m2)*a2); 
B = (p*omega2/m2)^p * (q*omega1/m1)^q * (A0*Il)^(a2*p);
z = B* (ur/x)^(a2*p/r);

an=[Delta(a2*p,1-xi^2) Delta(q,1-m1) Delta(p,1-m2)];
ap=[Delta(a2*p,1)];
bm=[Delta(a2*p,0)];
bq=[ Delta(a2*p,-xi^2)];
M=meijerG( an, ap, bm, bq, z );

output = A * M;



end
