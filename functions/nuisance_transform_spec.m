function S = nuisance_transform_spec()
%NUISANCE_TRANSFORM_SPEC
%
% Production nuisance transform specification.
%
%   B08  [0,0.5]         uniform    -> scaled_logit
%   W01  [7,21]          uniform    -> scaled_logit
%   I01  [0.5,5]         loguniform -> log
%   I02  [2,20]          loguniform -> log
%   I07  [0.03,0.09]     uniform    -> scaled_logit
%   I08  [10,70]         loguniform -> log
%   W02  [0.0009,0.0109] loguniform -> log


parameterID=["B08";"W01";"I01";"I02";"I07";"I08";"W02"];
physicalLower=[0;7;0.5;2;0.03;10;0.0009];
physicalUpper=[0.5;21;5;20;0.09;70;0.0109];
distribution=["uniform";"uniform";"loguniform";"loguniform"; ...
              "uniform";"loguniform";"loguniform"];
transform=["scaled_logit";"scaled_logit";"log";"log"; ...
           "scaled_logit";"log";"log"];
source=["gsa_isotope_parameter_definitions"; ...
        "run_J1_integrated_cached_gsa"; ...
        "gsa_isotope_parameter_definitions"; ...
        "gsa_isotope_parameter_definitions"; ...
        "gsa_isotope_parameter_definitions"; ...
        "gsa_isotope_parameter_definitions"; ...
        "run_J1_integrated_cached_gsa"];

samplerLower=nan(size(physicalLower));
samplerUpper=nan(size(physicalUpper));
jacobianRequired=false(size(physicalLower));

for i=1:numel(parameterID)
    switch transform(i)
        case "scaled_logit"
            samplerLower(i)=-Inf;
            samplerUpper(i)=Inf;
            jacobianRequired(i)=true;

        case "log"
            samplerLower(i)=log(physicalLower(i));
            samplerUpper(i)=log(physicalUpper(i));
            jacobianRequired(i)=false;

        otherwise
            error('Unsupported transform in production specification.');
    end
end

S=table(parameterID,physicalLower,physicalUpper,distribution,transform, ...
    samplerLower,samplerUpper,jacobianRequired,source);

end
