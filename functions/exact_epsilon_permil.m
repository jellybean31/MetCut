
function epsPermil = exact_epsilon_permil(deltaProduct,deltaWater)
%EXACT_EPSILON_PERMIL Exact isotope fractionation in ratio space.
epsPermil = 1000 .* ( ...
    (1 + deltaProduct./1000) ./ ...
    (1 + deltaWater./1000) - 1);
end
