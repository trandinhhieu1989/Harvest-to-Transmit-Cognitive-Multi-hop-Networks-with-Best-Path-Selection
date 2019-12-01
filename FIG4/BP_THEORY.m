function h_out = BP_THEORY(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha)
% BP_THEORY  : Theory of Best Path Protocol
% OP: Outage Probability
OP = zeros(length(PdB), length(xE));
%
for nn = 1 : length(PdB)
for aa = 1 : length(xE)
    OP(nn,aa) = BPfunc(PdB(nn),QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE(aa),yE,eta,alpha);
    if nn==1
    h2= semilogy(xE,OP,'r-'); grid on;hold on;
    else
    semilogy(xE,OP,'r-'); grid on;hold on;
    end
end
end
%
OP
%

h_out = h2;
%plot(alpha,OP,'r-'); grid on;hold on;
end
%
function out = BPfunc(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha)
% PdB       : Transmit power of beacons
% QdB       : Interference Constraints
% NN        : Number of Beacons
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
    %
    % the number of hops
    Hop        = LL(aa) + 1;
    % Define rho at each hop
    rho         = 2^(Hop*RR/(1 - alpha)) - 1;
    % Outage at each hop
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
        hs     = 0;
        for nnn = 0 : NN - 1
            gt1    = 2/factorial(nnn)*(OMB*LD*rho)^((nnn+1)/2)*besselk(nnn-1,2*sqrt(OMB*LD*rho));
            gt2    = 2/factorial(nnn)*(OMB)^((nnn+1)/2)*LD*rho*(LD*rho + OMP)^((nnn-1)/2)*besselk(nnn-1,2*sqrt(OMB*(LD*rho + OMP)));
            gt3    = 2/factorial(nnn)*(OMB)^((nnn+1)/2)*LD*rho*(LD*rho + OME)^((nnn-1)/2)*besselk(nnn-1,2*sqrt(OMB*(LD*rho + OME)));
            gt4    = 2/factorial(nnn)*(OMB)^((nnn+1)/2)*LD*rho*(LD*rho + OMP + OME)^((nnn-1)/2)*besselk(nnn-1,2*sqrt(OMB*(LD*rho + OMP + OME)));
            hs     = hs + gt1 - gt2 - gt3 + gt4;
        end
        %
        OP_hop = OP_hop*hs;
    end
    out = out*(1 - OP_hop);
end
end





