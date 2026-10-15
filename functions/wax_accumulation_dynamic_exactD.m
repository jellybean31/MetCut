
function A = wax_accumulation_dynamic_exactD( ...
    timeDays,productDeltaD_permil,alkaneCarbonNumber,T)
%WAX_ACCUMULATION_DYNAMIC_EXACTD
% Integrate standing leaf-wax mass and D inventory from instantaneous
% newly synthesized wax isotope composition.
%
% INPUTS
% ------
% timeDays
%   1 x nTime or nTime x 1 time vector in days.
%
% productDeltaD_permil
%   nHomologue x nTime instantaneous newly synthesized wax deltaD.
%
% alkaneCarbonNumber
%   nHomologue x 1 or 1 x nHomologue homologues (e.g. [27 29 31 33]).
%
% T
%   accumulation traits from wax_accumulation_default_traits().
%
% OUTPUT
% ------
% A.productDeltaD_permil
%   unchanged instantaneous biosynthetic product.
%
% A.poolDeltaD_permil
%   isotope composition of the accumulated standing wax pool.
%
% A.mass
%   standing wax mass in arbitrary internally consistent units.
%
% A.Dinventory
%   D atom inventory in the same mass-normalized units.
%
% A.productionRate
%   new wax production rate, mass d^-1.
%
% A.lossRate
%   first-order standing-pool loss rate, d^-1.

if nargin<4 || isempty(T)
    T=wax_accumulation_default_traits();
end

timeDays=timeDays(:).';
C=alkaneCarbonNumber(:);

Y=productDeltaD_permil;
if size(Y,2)~=numel(timeDays) && size(Y,1)==numel(timeDays)
    Y=Y.';
end

nH=numel(C);
nT=numel(timeDays);

if size(Y,1)~=nH || size(Y,2)~=nT
    error(['productDeltaD_permil must be nHomologue x nTime and agree ' ...
        'with alkaneCarbonNumber and timeDays.']);
end

if any(diff(timeDays)<=0)
    error('timeDays must be strictly increasing.');
end

S=wax_accumulation_production_schedule(timeDays,nH,T);

xProd=deltaD_to_Datomfraction(Y,T.Rstd_VSMOW);

M=NaN(nH,nT);
Q=NaN(nH,nT);
xPool=NaN(nH,nT);

M0=expand_homologue_vector(T.initialMass,nH,'initialMass');

if all(isnan(T.initialDeltaD_permil))
    x0=NaN(nH,1);
else
    d0=expand_homologue_vector(T.initialDeltaD_permil,nH, ...
        'initialDeltaD_permil');
    x0=deltaD_to_Datomfraction(d0,T.Rstd_VSMOW);
end

M(:,1)=M0;

for h=1:nH
    if M0(h)>T.massTolerance
        if ~isfinite(x0(h))
            error(['initialDeltaD_permil must be finite when initialMass ' ...
                'is nonzero.']);
        end
        Q(h,1)=M0(h)*x0(h);
        xPool(h,1)=x0(h);
    else
        Q(h,1)=0;
    end
end

for i=1:nT-1
    dt=timeDays(i+1)-timeDays(i);

    % Trapezoidal mass influx and D-atom influx.
    Pbar=0.5*(S.productionRate(:,i)+S.productionRate(:,i+1));
    DfluxBar=0.5*( ...
        S.productionRate(:,i).*xProd(:,i) + ...
        S.productionRate(:,i+1).*xProd(:,i+1));

    xIn=NaN(nH,1);
    idx=Pbar>0;
    xIn(idx)=DfluxBar(idx)./Pbar(idx);

    kbar=0.5*(S.lossRate(:,i)+S.lossRate(:,i+1));

    for h=1:nH
        k=kbar(h);
        p=Pbar(h);

        if k>0
            e=exp(-k*dt);
            gain=(1-e)/k;

            M(h,i+1)=M(h,i)*e + p*gain;

            if p>0
                Q(h,i+1)=Q(h,i)*e + p*xIn(h)*gain;
            else
                Q(h,i+1)=Q(h,i)*e;
            end
        else
            M(h,i+1)=M(h,i)+p*dt;

            if p>0
                Q(h,i+1)=Q(h,i)+p*xIn(h)*dt;
            else
                Q(h,i+1)=Q(h,i);
            end
        end

        if M(h,i+1)>T.massTolerance
            xPool(h,i+1)=Q(h,i+1)/M(h,i+1);
        end
    end
end

poolDeltaD=Datomfraction_to_deltaD(xPool,T.Rstd_VSMOW);

A=struct;
A.timeDays=timeDays;
A.alkaneCarbonNumber=C;
A.productDeltaD_permil=Y;
A.productDatomFraction=xProd;
A.poolDeltaD_permil=poolDeltaD;
A.poolDatomFraction=xPool;
A.mass=M;
A.Dinventory=Q;
A.productionRate=S.productionRate;
A.lossRate=S.lossRate;
A.referenceMatureMass=S.referenceMatureMass;
A.traits=T;

end

function x=deltaD_to_Datomfraction(d,Rstd)
R=Rstd.*(1+d./1000);
x=R./(1+R);
end

function d=Datomfraction_to_deltaD(x,Rstd)
R=x./(1-x);
d=1000.*(R./Rstd-1);
d(~isfinite(x))=NaN;
end

function v=expand_homologue_vector(x,n,name)
x=x(:);
if isscalar(x)
    v=repmat(x,n,1);
elseif numel(x)==n
    v=x;
else
    error('%s must be scalar or have one value per homologue.',name);
end
end
