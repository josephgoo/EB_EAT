dimA = 4;
dimB = 2;
ketPlus = 1/sqrt(2)*[1;1];
ketMinus = 1/sqrt(2)*[1;-1];
pz=0.9;
Q=0.3;

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

zeroPOVM=0;
onePOVM=0;
abortPOVM=0;
for a0=0:1
    for a1=0:1
        for b0=0:1
            for b1=0:1
                if a0==0 && a1==0
                    a=1;
                elseif a0==0 && a1==1
                    a=2;
                elseif a0==1 && a1==0
                    a=3;
                elseif a0==1 && a1==1
                    a=4;
                end

                if b0==0 && b1==0
                    b=1;
                elseif b0==0 && b1==1
                    b=2;
                elseif b0==1 && b1==0
                    b=3;
                elseif b0==1 && b1==1
                    b=4;
                end

                %calculate zeroPOVM
                if a0==b0 && a1~=b1 &&a0==1
                    zeroPOVM=zeroPOVM+kron(APOVM{a},BPOVM{b1+3});
                end
                
                %calculate onePOVM
                if a0==b0 && a1==b1 &&a0==1
                    onePOVM=onePOVM+kron(APOVM{a},BPOVM{a});
                end
            
                %calculate abortPOVM
                if a0==0 || b0==0
                    abortPOVM=abortPOVM+kron(APOVM{a},BPOVM{b});
                end
            end 
        end    
    end    
end         

testxpovm=zeroPOVM+onePOVM+abortPOVM;

krausOpsp={};
for a=1:4
    for b=1:4
        krausOpsp={krausOpsp{:},sqrtm(kron(APOVM{a},BPOVM{b}))};
    end
end


krausOp_00=0;
krausOp_01=0;
krausOp_10=0;
krausOp_11=0;
for a1=0:1
    for b1=0:1
        krausOp_00=krausOp_00+kron(APOVM{a1+1},BPOVM{b1+1});
        krausOp_01=krausOp_01+kron(APOVM{a1+1},BPOVM{b1+3});
        krausOp_10=krausOp_10+kron(APOVM{a1+3},BPOVM{b1+1});
        krausOp_11=krausOp_11+kron(APOVM{a1+3},BPOVM{b1+3});
    end
end
krausOp_00=sqrtm(krausOp_00);
krausOp_01=sqrtm(krausOp_01);
krausOp_10=sqrtm(krausOp_10);
krausOp_11=sqrtm(krausOp_11);
krausOp_p={krausOp_00,krausOp_01,krausOp_10,krausOp_11};


M0=(1-pz)^2*kron(zket(2,1)*zket(2,1)',ketMinus*ketMinus')+(1-pz)^2*kron(zket(2,2)*zket(2,2)',ketPlus*ketPlus');
M1=(1-pz)^2*kron(zket(2,1)*zket(2,1)',ketPlus*ketPlus')+(1-pz)^2*kron(zket(2,2)*zket(2,2)',ketMinus*ketMinus');
rho_sim = (1-(3/2)*Q)*(phiPlus*phiPlus')+(Q/2)*((phiMinus*phiMinus')+(psiPlus*psiPlus')+(psiMinus*psiMinus'));


M0=(1/(1-pz)^2)*M0;
M1=(1/(1-pz)^2)*M1;
real(trace(rho_sim*M0));
M0+M1;

rho_sim = (1-3/2*Q)*phiPlus*phiPlus'+Q/2*(phiMinus*phiMinus'+psiPlus*psiPlus'+psiMinus*psiMinus');

%8 dimension
rho=(pz)*kron(zket(4,1)*zket(4,1)',zket(2,1)*zket(2,1)')+...
    (pz)*kron(zket(4,2)*zket(4,2)',zket(2,2)*zket(2,2)')+...
    (1-pz)*kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+...
    (1-pz)*kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');
rho=rho/2;
trace(rho*(zeroPOVM+onePOVM))-(1-pz)^2;

M0=(1-pz)*kron(zket(4,3)*zket(4,3)',ketMinus*ketMinus')+(1-pz)*kron(zket(4,4)*zket(4,4)',ketPlus*ketPlus');
M1=(1-pz)*kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+(1-pz)*kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');
trace(rho*(M0+M1));
isequal(M0,zeroPOVM);

%4 dimension
% rho=kron(zket(2,1)*zket(2,1)',ketPlus*ketPlus')+...
%     kron(zket(2,2)*zket(2,2)',ketMinus*ketMinus');
% rho=rho/2;
% trace(rho);
% 
% M0=kron(zket(2,1)*zket(2,1)',ketMinus*ketMinus')+kron(zket(2,2)*zket(2,2)',ketPlus*ketPlus');
% M1=kron(zket(2,1)*zket(2,1)',ketPlus*ketPlus')+kron(zket(2,2)*zket(2,2)',ketMinus*ketMinus');
% M0+M1;

%8d imension modified
M0bar=(1-pz)*kron(zket(4,3)*zket(4,3)',ketMinus*ketMinus')+(1-pz)*kron(zket(4,4)*zket(4,4)',ketPlus*ketPlus');
M1bar=(1-pz)*kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+(1-pz)*kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');
Mperpbar=abortPOVM;


A=[];
A(1:2,1:2)=eye(2);
A(3:4,3:4)=rho_sim(1:2,1:2);
A(5:6,5:6)=eye(2);
A(7:8,7:8)=rho_sim(3:4,3:4);

A(3:4,7:8)=rho_sim(3:4,1:2);
A(7:8,3:4)=rho_sim(1:2,3:4);
trace(A*M0);

B=[];
B(1:4,1:4)=zeros(4);
B(5:8,5:8)=rho_sim;




%check positive semidefinite
X=M1;
symmetry = issymmetric(X);
[~,D]=eig(X);
eigenvalues = diag(D);
if all(eigenvalues>=0) & symmetry
    disp('Positive semi-definite matrix.');
else
     disp('Non positive semi-definite matrix.');
end




%8 dimen px=1
M0=kron(zket(4,3)*zket(4,3)',ketMinus*ketMinus')+kron(zket(4,4)*zket(4,4)',ketPlus*ketPlus');
M1=kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');
Mperp=eye(8)-M0-M1;
trace(rho*(M0+M1))
trace(rho*(zeroPOVM+onePOVM))
