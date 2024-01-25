dimA = 4;
dimB = 2;
ketPlus = 1/sqrt(2)*[1;1];
ketMinus = 1/sqrt(2)*[1;-1];
pz=0.1;
Q=0.1;

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

% rho0 = eye(prod(protocolDescription.dimensions));
% ksp=protocolDescription.krausOp_sp;
% kp=protocolDescription.krausOp_p;
% primalDf(rho0,ksp,kp)

for i=1:numel(krausOp_p)
    A=krausOp_p{i};
    isequal(A*rho0*A',A'*rho0*A);
end


A=[1,2]
A(1)

