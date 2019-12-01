clear all;  clc; closereq;
tic;
PdB           = -5:5:25;
QdB           = 5;
LL            = [3 4 5];
PL            = 3;
RR            = 0.5;
%
NN            = 2;
xB            = 0.5;
yB            = 0.25;
%
xP            = 0.5;
yP            = -1;
%
xE            = 0.5;
yE            = 0.5;
%
eta           = 0.1;
% time fraction
alpha         = [0.1 0.3 0.5];
%
Num_Trial     = 10^5;
%
% RPS Protocol
%RP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial);
%RP_THEORY(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
%RP_ASYM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
% SPS Protocol
%SP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial);
%SP_THEORY(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
%SP_ASYM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
% BPS Protocol
h1= BP_SIM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha,Num_Trial);
h2= BP_THEORY(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
h3 = BP_ASYM(PdB,QdB,LL,PL,RR,NN,xB,yB,xP,yP,xE,yE,eta,alpha);
legend([h1(1,1) h2(1,1) h3(1,1) ],{'BRS (Sim)','BRS (Exact)','BRS (Asymp)'});
xlabel('Q (dB)');
ylabel('Outage probability (OP)');
toc
