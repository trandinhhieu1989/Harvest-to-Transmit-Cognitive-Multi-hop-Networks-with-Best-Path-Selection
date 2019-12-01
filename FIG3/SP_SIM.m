function SP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial)    
% SP_SIM    : Simulation of Shorstest Path Protocol
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
OP          = zeros(1,length(alpha));
% Define kappa
kp          = eta*alpha./(1-alpha);
%
for aa = 1 : length(alpha)
    fprintf('Running %d per %d \n',aa,length(alpha));
    % Select shortest path, Lmin is the number of hops
    Lmin    = min(LL) + 1;
    % Define rho
    rho     = 2^(Lmin*RR/(1 - alpha(aa))) - 1;     
    for bitnum   =  1 : Num_Trial       
        % Find end-to-end SNR
        SNRmin   = inf;
        for bb = 1 : Lmin
            % Parameter of data links: Lambda_D                         
            LD    = (1/Lmin)^PL;
            % Parameter of energy harvesting links: Lambda_B
            LB    = sqrt(((bb-1)/Lmin - xB)^2 + yB^2)^PL;            
            % Parameter of interference links: Lambda_P
            LP    = sqrt(((bb-1)/Lmin - xP)^2 + yP^2)^PL;
            % Parameter of eavesdopping links: Lambda_E
            LE    = sqrt(((bb-1)/Lmin - xE)^2 + yE^2)^PL;
            % Rayleigh fading channel of data links
            hD    = sqrt(1/LD/2)*(randn(1,1)+j*randn(1,1));
            % Power from energy harvesting
            P_EH  = 0;
            for dd = 1 : NN
                % Rayleigh fading channel of energy harvesting links
                hB   = sqrt(1/LB/2)*(randn(1,1)+j*randn(1,1));
                P_EH = P_EH + kp(aa)*PP*abs(hB)^2;
            end                                     
            % Rayleigh fading channel of interference links
            hP    = sqrt(1/LP/2)*(randn(1,1)+j*randn(1,1));
            % Rayleigh fading channel of eavesdopping links
            hE    = sqrt(1/LE/2)*(randn(1,1)+j*randn(1,1));            
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
        CC       = (1-alpha(aa))/Lmin*log2(1 + SNRmin);                         
        %      
        if (CC <= RR)
            OP(aa) = OP(aa) + 1;
        end
    end
end
%
OP = OP./Num_Trial;
%
OP
%
plot(alpha,OP,'bs'); grid on;hold on;
end





