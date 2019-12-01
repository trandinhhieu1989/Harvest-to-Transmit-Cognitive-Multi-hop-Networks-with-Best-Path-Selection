function h_out = BP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial)
% BP_ASYM  : Asymptotic OP for Best Path Protocol
% OP: Outage Probability
OP = zeros(length(PdB), length(xE));
%
for nn = 1 : length(PdB)
for aa = 1 : length(xE)
    OP(nn,aa) = BPfunc(PdB(nn),QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE(aa),yE,eta,alpha,Num_Trial);
    if nn ==1
    h1 = semilogy(xE,OP,'ro'); grid on;hold on;
    else
    semilogy(xE,OP,'ro'); grid on;hold on;
    end
end
end
%
OP
%

h_out = h1;
end
function out =BPfunc(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial)
% BP_SIM    : Simulation of Shorstest Path Protocol
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
% OP: Outage Probability
OP          = zeros(1,length(xE));
% Define kappa
kp          = eta*alpha/(1-alpha);
%
for aa = 1 : length(xE)
    fprintf('Running %d per %d \n',aa,length(xE));
    for bitnum   =  1 : Num_Trial
        % Find best path with highest capacity
        CC_max  = 0;
        for bb = 1 : length(LL)
            % Number of hops
            Hop  = LL(bb) + 1;
            % rho at each hop
            rho  = 2^(Hop*RR/(1 - alpha)) - 1;
            % Find end-to-end SNR
            SNRmin   = inf;
            %
            for cc = 1 : Hop
                % Parameter of data links: Lambda_D
                LD    = (1/Hop)^PL;
                % Parameter of energy harvesting links: Lambda_B
                LB    = sqrt(((cc-1)/Hop - xB)^2 + yB^2)^PL;
                % Parameter of interference links: Lambda_P
                LP    = sqrt(((cc-1)/Hop - xP)^2 + yP^2)^PL;
                % Parameter of eavesdopping links: Lambda_E
                LE    = sqrt(((cc-1)/Hop - xE(aa))^2 + yE^2)^PL;
                % Rayleigh fading channel of data links
                hD    = sqrt(1/LD/2)*(randn(1,1)+j*randn(1,1));               
                % Rayleigh fading channel of interference links
                hP    = sqrt(1/LP/2)*(randn(1,1)+j*randn(1,1));
                % Rayleigh fading channel of eavesdopping links
                hE    = sqrt(1/LE/2)*(randn(1,1)+j*randn(1,1));
                % Power from energy harvesting
                P_EH  = 0;
                for dd = 1 : NN
                    % Rayleigh fading channel of energy harvesting links
                    hB   = sqrt(1/LB/2)*(randn(1,1)+j*randn(1,1));
                    P_EH = P_EH + kp*PP*abs(hB)^2;
                end
                % Power constrained by interference
                P_Int = QQ/abs(hP)^2;
                % Power constrained due to the presence of eavesdopper
                P_Eav = rho/abs(hE)^2;
                % Transmit power
                P_t   = min(P_EH, min (P_Int, P_Eav));
                % Signal-to-noise ratio (SNR)
                SNR   = P_t*abs(hD)^2;
                if (SNR < SNRmin)
                    SNRmin = SNR;
                end
            end
            % End-to-end channel capacity
            CC       = (1-alpha)/Hop*log2(1 + SNRmin);
            if (CC > CC_max)
                CC_max = CC;
            end
        end
        %
        if (CC_max <= RR)
            OP(aa) = OP(aa) + 1;
        end
    end
end
%
out = OP./Num_Trial;
%
%OP
%
%plot(xE,OP,'ro'); grid on;hold on;
end





