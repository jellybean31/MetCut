
function anglesDeg = subspace_principal_angles_deg(A,B)
%SUBSPACE_PRINCIPAL_ANGLES_DEG Principal angles between column subspaces.
%
% A and B may have non-orthonormal columns.

[QA,~]=qr(A,0);
[QB,~]=qr(B,0);

sv=svd(QA.'*QB);
sv=min(max(sv,-1),1);

anglesDeg=acosd(sv);

end
