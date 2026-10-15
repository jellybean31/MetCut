function M = map_geometric_two_surface_to_anatomical( ...
    geometricValues,upperAnatomicalSide,inclinationFromHorizontal_deg)
% MAP_GEOMETRIC_TWO_SURFACE_TO_ANATOMICAL
%
% Generic production mapper from geometric two-surface order
%
%     [LOWER, UPPER]
%
% to the model anatomical order
%
%     [ABAXIAL, ADAXIAL].
%
% INPUTS
% ------
% geometricValues
%     Numeric, real, finite two-element vector in [LOWER,UPPER] order.
%
% upperAnatomicalSide
%     "ADAXIAL", "ABAXIAL", or "UNRESOLVED".
%
% inclinationFromHorizontal_deg
%     Scalar in [0,90], with 0 horizontal and 90 vertical.
%
% Surface semantics:
%     ADAXIAL upper:
%         geometric [LOWER,UPPER] -> anatomical [ABAXIAL,ADAXIAL]
%         by identity.
%
%     ABAXIAL upper:
%         geometric [LOWER,UPPER] -> anatomical [ABAXIAL,ADAXIAL]
%         by swap.


M=struct();
M.success=false;
M.code="UNSET";
M.sourceOrder=["LOWER","UPPER"];
M.outputOrder=["ABAXIAL","ADAXIAL"];
M.upperAnatomicalSide="";
M.facingResolved=false;
M.verticalDegenerate=false;
M.geometricValues=[NaN NaN];
M.anatomicalValues=[NaN NaN];
M.permutationMatrix=NaN(2);

if ~(isnumeric(geometricValues) && isreal(geometricValues) && ...
        numel(geometricValues)==2 && all(isfinite(geometricValues(:))))
    M.code="INVALID_GEOMETRIC_VALUES";
    return
end
x=reshape(double(geometricValues),1,2);

if ~(isnumeric(inclinationFromHorizontal_deg) && ...
        isreal(inclinationFromHorizontal_deg) && ...
        isscalar(inclinationFromHorizontal_deg) && ...
        isfinite(inclinationFromHorizontal_deg) && ...
        inclinationFromHorizontal_deg>=0 && inclinationFromHorizontal_deg<=90)
    M.code="INVALID_INCLINATION";
    return
end

try
    state=upper(string(upperAnatomicalSide));
catch
    M.code="INVALID_FACING_STATE";
    return
end
if ~isscalar(state)
    M.code="INVALID_FACING_STATE";
    return
end

M.geometricValues=x;
M.upperAnatomicalSide=state;

Pidentity=[1 0;0 1];
Pswap=[0 1;1 0];

switch state
    case "ADAXIAL"
        M.permutationMatrix=Pidentity;
        M.anatomicalValues=x*Pidentity;
        M.facingResolved=true;
        M.code="OK";
        M.success=true;

    case "ABAXIAL"
        M.permutationMatrix=Pswap;
        M.anatomicalValues=x*Pswap;
        M.facingResolved=true;
        M.code="OK";
        M.success=true;

    case "UNRESOLVED"
        % The unresolved state is admissible only at the exactly vertical
        % endpoint and only when the geometric two-surface quantity is
        % numerically degenerate.
        tolAngle=64*eps(max(abs(inclinationFromHorizontal_deg),1));
        if abs(inclinationFromHorizontal_deg-90)>tolAngle
            M.code="NONVERTICAL_FACING_REQUIRED";
            return
        end

        tolValue=1024*eps(max([abs(x) 1]));
        if abs(x(1)-x(2))>tolValue
            M.code="VERTICAL_VALUES_NOT_DEGENERATE";
            return
        end

        % Symmetric roundoff cleanup only. Anatomical facing remains
        % unresolved and is explicitly retained as such.
        common=0.5*(x(1)+x(2));
        M.anatomicalValues=[common common];
        M.permutationMatrix=NaN(2);
        M.facingResolved=false;
        M.verticalDegenerate=true;
        M.code="VERTICAL_DEGENERATE_UNRESOLVED_FACING";
        M.success=true;

    otherwise
        M.code="INVALID_FACING_STATE";
        return
end

M.provenance=struct();
M.provenance.scope="GEOMETRIC_TO_ANATOMICAL_TWO_SURFACE_MAPPING";
M.provenance.implementationRole="PRODUCTION_GENERIC";
M.provenance.B2020ConditioningEmbedded=false;
M.provenance.defaultFacingEmbedded=false;
end
