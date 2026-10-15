
function W2 = apply_wax_er_stage_states_exactD( ...
    waxDeltaDMatrix,alkaneCarbonNumber,traits)
%APPLY_WAX_ER_STAGE_STATES_EXACTD Architecture-v2 effective stage adapter.

if nargin<3 || isempty(traits)
    traits=wax_er_stage_state_traits();
end

D=waxDeltaDMatrix;
C=alkaneCarbonNumber(:);

if size(D,1)~=numel(C)
    error('Row count of waxDeltaDMatrix must match alkaneCarbonNumber.');
end

if ~traits.enabled || ...
        (traits.stage6Delta_permil==0 && traits.stage7Delta_permil==0)
    W2=pack_output(D,D,C,traits);
    return
end

Rstd=traits.Rstd_VSMOW;
Hstage=traits.hydrogenAtomsPerStage;

Rbase=Rstd.*(1+D./1000);
xbase=Rbase./(1+Rbase);
xnew=xbase;

for i=1:numel(C)
    nH=2*C(i)+2;
    xi=xbase(i,:);

    if C(i)>=traits.stage6MinimumCarbonNumber
        xi=apply_one_stage( ...
            xi,nH,Hstage,traits.stage6Delta_permil);
    end
    if C(i)>=traits.stage7MinimumCarbonNumber
        xi=apply_one_stage( ...
            xi,nH,Hstage,traits.stage7Delta_permil);
    end

    xnew(i,:)=xi;
end

Rnew=xnew./(1-xnew);
Dnew=1000.*(Rnew./Rstd-1);

W2=pack_output(Dnew,D,C,traits);
W2.DatomFraction=xnew;
W2.ratioMatrix=Rnew./Rstd;
end

function xout=apply_one_stage(xin,nH,nStage,deltaStage)
Rin=xin./(1-xin);
alphaStage=1+deltaStage/1000;
Rstage=Rin.*alphaStage;
xstage=Rstage./(1+Rstage);
xout=((nH-nStage).*xin+nStage.*xstage)./nH;
end

function W=pack_output(Dnew,Dbase,C,T)
W=struct;
W.deltaDMatrix=Dnew;
W.baselineDeltaDMatrix=Dbase;
W.deltaFromBaseline_permil=Dnew-Dbase;
W.alkaneCarbonNumber=C;
W.traits=T;
end
