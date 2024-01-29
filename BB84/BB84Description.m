%% FUNCTION NAME: BB84Description
% Simple description for entanglement-based BB84. 
%%

function protocolDescription = BB84Description(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    dimA = 4;
    dimB = 2;
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];

    % Alice and Bob POVM
    aPOVM={};
    for i=1:dimA
        aPOVM={aPOVM{:},zket(4,i) * zket(4,i)'};
    end
    bPOVM={pz*zket(2,1)*zket(2,1)',pz*zket(2,2)*zket(2,2)',(1-pz)*ketPlus*ketPlus',(1-pz)*ketMinus*ketMinus'};
    
    %kraus operator for joint secret and public information
    krausOp_sp={};
    for a=1:4
        for b=1:4
            krausOp_sp={krausOp_sp{:},sqrtm(kron(aPOVM{a},bPOVM{b}))};
        end
    end


    %kraus operator for public information
    krausOp_00=0;
    krausOp_01=0;
    krausOp_10=0;
    krausOp_11=0;
    for a1=0:1
        for b1=0:1
            krausOp_00=krausOp_00+kron(aPOVM{a1+1},bPOVM{b1+1});
            krausOp_01=krausOp_01+kron(aPOVM{a1+1},bPOVM{b1+3});
            krausOp_10=krausOp_10+kron(aPOVM{a1+3},bPOVM{b1+1});
            krausOp_11=krausOp_11+kron(aPOVM{a1+3},bPOVM{b1+3});
        end
    end
    krausOp_00=sqrtm(krausOp_00);
    krausOp_01=sqrtm(krausOp_01);
    krausOp_10=sqrtm(krausOp_10);
    krausOp_11=sqrtm(krausOp_11);

    krausOp_p={krausOp_00,krausOp_01,krausOp_10,krausOp_11};

    % Constraints xPOVM
    % zeroPOVM=0;
    % onePOVM=0;
    % abortPOVM=0;
    % for a0=0:1
    %     for a1=0:1
    %         for b0=0:1
    %             for b1=0:1
    %                 if a0==0 && a1==0
    %                     a=1;
    %                 elseif a0==0 && a1==1
    %                     a=2;
    %                 elseif a0==1 && a1==0
    %                     a=3;
    %                 elseif a0==1 && a1==1
    %                     a=4;
    %                 end
    % 
    %                 if b0==0 && b1==0
    %                     b=1;
    %                 elseif b0==0 && b1==1
    %                     b=2;
    %                 elseif b0==1 && b1==0
    %                     b=3;
    %                 elseif b0==1 && b1==1
    %                     b=4;
    %                 end
    % 
    %                 %calculate zeroPOVM
    %                 if a0==b0 && a1~=b1 && a0==1
    %                     zeroPOVM=zeroPOVM+kron(aPOVM{a},bPOVM{b1+3});
    %                 end
    % 
    %                 %calculate onePOVM
    %                 if a0==b0 && a1==b1 && a0==1
    %                     onePOVM=onePOVM+kron(aPOVM{a},bPOVM{a});
    %                 end
    % 
    %                 %calculate abortPOVM
    %                 if a0==0 || b0==0
    %                     abortPOVM=abortPOVM+kron(aPOVM{a},bPOVM{b});
    %                 end
    %             end
    %         end
    %     end
    % end
    % 

    zeroPOVM=kron(zket(4,3)*zket(4,3)',ketMinus*ketMinus')+kron(zket(4,4)*zket(4,4)',ketPlus*ketPlus');
    onePOVM=kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');

    addObservables(zeroPOVM);
    addObservables(onePOVM);
    %addObservables(abortPOVM);
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    protocolDescription.observables = observables;
    protocolDescription.krausOp_sp = krausOp_sp;
    protocolDescription.krausOp_p = krausOp_p;
    protocolDescription.dimensions = [dimA,dimB];

end