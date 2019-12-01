function SP_ASYM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha)
% SP_ASYM  : Asymptotic OP for Shorstest Path Protocol
% OP: Outage Probability
OP = zeros(1, length(PdB));
%
for aa = 1 : length(PdB)
    OP(aa) = SPfunc(PdB(aa),QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
end
%
OP
%
plot(PdB,OP,'b--'); grid on;hold on;
end
%
function out = SPfunc(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha)
% PdB       : Transmit power of beacons
% QdB       : Interference Constraints
% MM        : Number of Paths
% LL        : a vectors including the number of intermediate nodes on each path
% PL        : Path-Loss
% RR        : Target Rate
% xB, yB    : co-ordinates of Beacons
% xP, yP    : co-ordinates of Primary Users
% xE, yE    : co-ordinates of Eavesdopper
% eta       : energy harvesting efficiency
% alpha     : fraction of time for energy harvesting
% Num_Trial : Number of Trials
% From dB to Watt
PP          = 10.^(PdB/10);
QQ          = 10.^(QdB/10);
% Define kappa
kp          = eta*alpha/(1-alpha);
% Select shortest path, Lmin is the number of hops
Lmin        = min(LL) + 1;
% Define rho
rho         = 2^(Lmin*RR/(1 - alpha)) - 1;
% OP: Outage Probability
out         = 1;
%
for bb = 1 : Lmin    
    % Parameter of data links: Lambda_D and Omega_D
    LD     = (1/Lmin)^PL;
    % Parameter of energy harvesting links: Lambda_B and Omega_B
    LB     = sqrt(((bb-1)/Lmin - xB)^2 + yB^2)^PL;
    OMB    = LB/PP/kp;
    % Parameter of interference links: Lambda_P and Omega_P
    LP     = sqrt(((bb-1)/Lmin - xP)^2 + yP^2)^PL;
    OMP    = LP*QQ;
    % Parameter of eavesdopping links: Lambda_E and Omega_E
    LE     = sqrt(((bb-1)/Lmin - xE)^2 + yE^2)^PL;
    OME    = LE*rho;
    %     
    out    = out*( 1 - LD*rho/(LD*rho + OMP) -  LD*rho/(LD*rho + OME) + LD*rho/(LD*rho + OMP + OME) );    
end
%
out = 1 - out;
end





