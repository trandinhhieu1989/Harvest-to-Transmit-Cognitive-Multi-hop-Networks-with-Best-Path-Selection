clear all;  clc; closereq;
PdB           = [5 7 9];
QdB           = 7.5;
LL            = [3 4 5];
PL            = 3;
RR            = 0.5;
%
NN            = 4;
xB            = 0.5;
yB            = 0.2;
%
xP            = 0.5;
yP            = -1;
%
xE            = 0.5;
yE            = 0.7;
%
eta           = 0.5;
alpha         = 0.02 : 0.02 : 0.06;
%
Num_Trial     = 10^6;
%

% BPS Protocol
h1= BP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial);
h2= BP_THEORY(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
%h3 = BP_ASYM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
legend([h1(1,1) h2(1,1)],{'BRS (Sim)','BRS (Exact)'});
xlabel('Q (dB)');
ylabel('Outage probability (OP)');
toc