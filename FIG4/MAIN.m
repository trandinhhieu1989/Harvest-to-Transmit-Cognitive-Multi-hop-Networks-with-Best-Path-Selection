clear all;  clc; closereq;
tic;
PdB           = [1 3 5];
QdB           = 5;
LL            = [3 4 5];
PL            = 3;
RR            = 0.5;
%
NN            = 1;
xB            = 0.5;
yB            = 0.1;
%
xP            = 0.5;
yP            = -1;
%
xE            = 0 : 0.1 : 1;
yE            = 0.4;
%
eta           = 1;
alpha         = 0.1;
%
Num_Trial     = 5*10^5;
%

% BPS Protocol
h1= BP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial);
h2= BP_THEORY(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
legend([h1(1,1) h2(1,1)],{'BRS (Sim)','BRS (Exact)'});
xlabel('Q (dB)');
ylabel('Outage probability (OP)');
toc