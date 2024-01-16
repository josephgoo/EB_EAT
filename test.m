dimA = 4;
dimB = 2;
ketPlus = 1/sqrt(2)*[1;1];
ketMinus = 1/sqrt(2)*[1;-1];

APOVM={};
for i=1:dimA
    APOVM={APOVM{:},zket(4,i) * zket(4,i)'};
end

BPOVM={1/2*zket(2,1)*zket(2,1)',1/2*zket(2,2)*zket(2,2)',1/2*ketPlus*ketPlus',1/2*ketMinus*ketMinus'};

onePOVM=0;
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
                
                if a0==b0 && a1==b1
                    onePOVM=onePOVM+kron(APOVM{a},BPOVM{a})
                end
            end 
        end    
    end    
end                    

zeroPOVM=0;
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
                
                if a0==b0 && a1~=b1
                    if a0==0
                        zeroPOVM=zeroPOVM+kron(APOVM{a},BPOVM{b1+1})
                    elseif a0==1
                        zeroPOVM=zeroPOVM+kron(APOVM{a},BPOVM{b1+3})
                    end
                end
            end 
        end    
    end    
end   

abortPOVM=0
for a0=0:1
    for a1=0:1
        for b0=0:1
            for b1=0:1
                if a0==0 && a1==0
                    a=1;
                    b=1;
                elseif a0==0 && a1==1
                    a=2;
                    b=2;
                elseif a0==1 && a1==0
                    a=3;
                    b=3;
                elseif a0==1 && a1==1
                    a=4;
                    b=4;
                end
                
                if a0~=b0
                    abortPOVM=abortPOVM+kron(APOVM{a},BPOVM{b})
                end
            end 
        end    
    end    
end 

testsum=zeroPOVM+onePOVM+abortPOVM





