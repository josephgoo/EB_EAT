dimA = 4;
dimB = 2;
ketPlus = 1/sqrt(2)*[1;1];
ketMinus = 1/sqrt(2)*[1;-1];
pz=0.9;
Q=0.19;

%Bell states
phiPlus = 1/sqrt(2)*[1;0;0;1];
phiMinus = 1/sqrt(2)*[1;0;0;-1];
psiPlus = 1/sqrt(2)*[0;1;1;0];
psiMinus = 1/sqrt(2)*[0;-1;1;0];



APOVM={};
for i=1:dimA
    APOVM={APOVM{:},zket(4,i) * zket(4,i)'};
end

BPOVM={pz*zket(2,1)*zket(2,1)',pz*zket(2,2)*zket(2,2)',(1-pz)*ketPlus*ketPlus',(1-pz)*ketMinus*ketMinus'};



M0=(1-pz)^2*kron(zket(2,1)*zket(2,1)',ketMinus*ketMinus')+(1-pz)^2*kron(zket(2,2)*zket(2,2)',ketPlus*ketPlus');
M1=(1-pz)^2*kron(zket(2,1)*zket(2,1)',ketPlus*ketPlus')+(1-pz)^2*kron(zket(2,2)*zket(2,2)',ketMinus*ketMinus');
rho_sim = (1-(3/2)*Q)*(phiPlus*phiPlus')+(Q/2)*((phiMinus*phiMinus')+(psiPlus*psiPlus')+(psiMinus*psiMinus'));


M0=(1/(1-pz)^2)*M0;
M1=(1/(1-pz)^2)*M1;
trace(rho_sim*M0);
trace(rho_sim*M1);
M0+M1;

rho=(1/2)^2*kron(zket(4,1)*zket(4,1)',zket(2,1)*zket(2,1)')+...
    (1/2)^2*kron(zket(4,2)*zket(4,2)',zket(2,2)*zket(2,2)')+...
    (1/2)^2*kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+...
    (1/2)^2*kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');

real(trace(rho*kron()))






