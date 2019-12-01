function h_out = BP_ASYM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha)
% BP_ASYM  : Asymptotic OP for Best Path Protocol
% OP: Outage Probability
OP_ASYMP = zeros(length(alpha), length(PdB));
%
for nn = 1 : length(alpha)
for aa = 1 : length(PdB)
   OP_ASYMP(nn,aa) = BPfunc(PdB(aa),QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha(nn));
    if nn ==1
    h3 = semilogy(PdB,OP_ASYMP ,'r--'); grid on;hold on;
    else
    semilogy(PdB,OP_ASYMP,'r--'); grid on;hold on;  
   end
end
end
%
OP_ASYMP 

h_out = h3;
end
%
function out = BPfunc(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha)
% PdB       : Transmit power of beacons
% QdB       : Interference Constraints
% MM        : Number of Beacons
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
% OP: Outage Probability
out         = 1;
%
for aa = 1 : length(LL)
    % Number of Hops
    Hop         = LL(aa) + 1;
    % Define rho
    rho         = 2^(Hop*RR/(1 - alpha)) - 1;
    %
    OP_hop      = 1; 
    for bb = 1 : Hop
        % Parameter of data links: Lambda_D and Omega_D
        LD     = (1/Hop)^PL;
        % Parameter of energy harvesting links: Lambda_B and Omega_B
        LB     = sqrt(((bb-1)/Hop - xB)^2 + yB^2)^PL;
        OMB    = LB/PP/kp;
        % Parameter of interference links: Lambda_P and Omega_P
        LP     = sqrt(((bb-1)/Hop - xP)^2 + yP^2)^PL;
        OMP    = LP*QQ;
        % Parameter of eavesdopping links: Lambda_E and Omega_E
        LE     = sqrt(((bb-1)/Hop - xE)^2 + yE^2)^PL;
        OME    = LE*rho;
        %
        OP_hop = OP_hop*( 1 - LD*rho/(LD*rho + OMP) -  LD*rho/(LD*rho + OME) + LD*rho/(LD*rho + OMP + OME) );
    end
    out = out*(1 - OP_hop);
end
%
end





