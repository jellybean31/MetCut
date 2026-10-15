
function Y = full_chain_wax_only_candidate(Mreference,Upstream,L,B)
%FULL_CHAIN_WAX_ONLY_CANDIDATE
% Evaluate candidate wax trajectory for one I02 upstream state and latent
% biochemical state L=[L1;L2].

E=latent2_to_effective5(L(:),B);
Bio=call_effective5_biochemistry_current(Mreference,Upstream.model,E);
Y=extract_wax_matrix_from_bio(Bio);

end
