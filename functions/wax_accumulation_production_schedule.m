
function S = wax_accumulation_production_schedule(timeDays,nHomologue,T)
%WAX_ACCUMULATION_PRODUCTION_SCHEDULE
% Construct wax production and loss schedules.
%
% Output:
%   S.productionRate       [nHomologue x nTime] mass d^-1
%   S.lossRate             [nHomologue x nTime] d^-1
%   S.referenceMatureMass  [nHomologue x 1]
%
% -------------------
% If developmentWindowDays <= 0, juvenile production is identically zero.
% This makes a zero development window a true no-production condition and
% avoids a spurious t=0 production impulse under trapezoidal integration.

timeDays=timeDays(:).';
nT=numel(timeDays);

prodMult=expand_homologue_vector(T.productionMultiplier,nHomologue, ...
    'productionMultiplier');
lossMult=expand_homologue_vector(T.lossMultiplier,nHomologue, ...
    'lossMultiplier');

peak=T.peakProductionRate_perDay .* prodMult;

if T.developmentWindowDays>0
    referenceMatureMass=peak .* T.developmentWindowDays;
else
    referenceMatureMass=zeros(nHomologue,1);
end

% Juvenile production envelope.
if T.developmentWindowDays<=0
    juvenile=zeros(1,nT);
else
    switch string(T.developmentShape)
        case "flat"
            juvenile=double(timeDays>=0 & ...
                timeDays<=T.developmentWindowDays);

        case "sine"
            juvenile=zeros(1,nT);
            idx=timeDays>=0 & timeDays<=T.developmentWindowDays;
            juvenile(idx)=sin(pi*timeDays(idx)/T.developmentWindowDays);

        otherwise
            error('Unknown developmentShape: %s',string(T.developmentShape));
    end
end

P=zeros(nHomologue,nT);

switch string(T.productionMode)
    case "continuous"
        P=peak*ones(1,nT);

    case "developmentGated"
        P=peak*juvenile;

    case "developmentPlusMaintenance"
        P=peak*juvenile;

        % Mature maintenance only exists when there is a defined mature
        % reference pool from a positive juvenile development window.
        if T.developmentWindowDays>0
            mature=timeDays>T.developmentWindowDays;
            matureRate=(T.matureReplacementFractionPerDay .* ...
                referenceMatureMass);
            P(:,mature)=matureRate*ones(1,sum(mature));
        end

    otherwise
        error('Unknown productionMode: %s',string(T.productionMode));
end

K=(T.lossRatePerDay.*lossMult)*ones(1,nT);

S=struct;
S.productionRate=P;
S.lossRate=K;
S.referenceMatureMass=referenceMatureMass;
S.juvenileShape=juvenile;

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
