function Q = apply_nuisance_transform(parameterID,value,direction)
%APPLY_NUISANCE_TRANSFORM
%
% Apply one production Stage-N nuisance transform.
%
% INPUT
% -----
% parameterID : B08, W01, I01, I02, I07, I08, or W02
% value       : numeric scalar/vector
% direction   : "physical_to_sampler" or "sampler_to_physical"
%
% OUTPUT
% ------
% Q.value
% Q.valid
% Q.logAbsDxDz
% Q.parameterID
% Q.transform
% Q.physicalLower / Q.physicalUpper
% Q.samplerLower / Q.samplerUpper
%
% ------
% - No clipping.
% - Invalid physical support returns valid=false and NaN transformed value.
% - Invalid bounded sampler support returns valid=false and NaN physical value.
% - Scaled-logit sampler values are allowed on all finite R.
% - Exact physical endpoints for scaled-logit are invalid because they map to
%   +/-Inf; finite sampler coordinates map strictly inside physical support.


S=nuisance_transform_spec();

id=string(parameterID);
j=find(S.parameterID==id,1);

if isempty(j)
    error('Unknown nuisance parameterID: %s',id);
end

x=value;
if ~isnumeric(x)
    error('value must be numeric.');
end

dir=string(direction);
tr=S.transform(j);

lo=S.physicalLower(j);
hi=S.physicalUpper(j);
slo=S.samplerLower(j);
shi=S.samplerUpper(j);

out=nan(size(x));
valid=isfinite(x);
logJac=nan(size(x));

switch dir

    case "physical_to_sampler"

        switch tr
            case "scaled_logit"
                valid = valid & x>lo & x<hi;

                if any(valid,'all')
                    xv=x(valid);
                    u=(xv-lo)./(hi-lo);
                    out(valid)=log(u)-log1p(-u);

                    % dx/dz = (hi-lo) * u * (1-u)
                    logJac(valid)= ...
                        log(hi-lo)+log(u)+log1p(-u);
                end

            case "log"
                valid = valid & x>=lo & x<=hi & x>0;

                if any(valid,'all')
                    out(valid)=log(x(valid));

                    % Existing loguniform prior is represented directly in
                    % y=log(x), so Stage N does not add a Jacobian term.
                    logJac(valid)=0;
                end

            otherwise
                error('Unsupported transform.');
        end

    case "sampler_to_physical"

        switch tr
            case "scaled_logit"
                % Any finite z is valid.
                valid=isfinite(x);

                if any(valid,'all')
                    z=x(valid);

                    % Numerically stable logistic.
                    u=zeros(size(z));
                    pos=z>=0;
                    u(pos)=1./(1+exp(-z(pos)));
                    ez=exp(z(~pos));
                    u(~pos)=ez./(1+ez);

                    out(valid)=lo+(hi-lo).*u;

                    % dx/dz.
                    logJac(valid)= ...
                        log(hi-lo)+log(u)+log1p(-u);
                end

            case "log"
                valid = isfinite(x) & x>=slo & x<=shi;

                if any(valid,'all')
                    out(valid)=exp(x(valid));
                    logJac(valid)=0;
                end

            otherwise
                error('Unsupported transform.');
        end

    otherwise
        error('direction must be physical_to_sampler or sampler_to_physical.');
end

Q=struct;
Q.version="N2_nuisance_transform_production_v1";
Q.parameterID=id;
Q.direction=dir;
Q.transform=tr;
Q.value=out;
Q.valid=valid;
Q.logAbsDxDz=logJac;
Q.physicalLower=lo;
Q.physicalUpper=hi;
Q.samplerLower=slo;
Q.samplerUpper=shi;
Q.jacobianRequired=S.jacobianRequired(j);

end
