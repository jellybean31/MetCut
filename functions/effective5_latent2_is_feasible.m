
function tf = effective5_latent2_is_feasible(L,B)
%EFFECTIVE5_LATENT2_IS_FEASIBLE

E=latent2_to_effective5(L,B);

lo=[B.parameterDefinitions.low].';
hi=[B.parameterDefinitions.high].';

tol=1e-12;
tf=all(E>=lo-tol & E<=hi+tol);

end
