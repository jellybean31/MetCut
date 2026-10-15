function R = test_METCUT_reaction_incorporation_exactD()
fprintf('\n============================================================\n');
fprintf('METCUT UNIT: REACTION INCORPORATION EXACT-D\n');
fprintf('============================================================\n');

Rstd=155.76e-6;
xToDelta=@(x)1000.*(((x./(1-x))./Rstd)-1);

t=(0:1/24:1).';
I=struct;
I.deltaD.plastid=-75.*ones(size(t));
I.deltaD.cytosol=-85.*ones(size(t));

Ts=wax_hydrogen_source_traits_postmalonate_null();
Ts.acetylStarterH.source='prescribed';
Ts.acetylStarterH.deltaD_permil=-120;
Ts.acetylStarterH.ratioNormalized=[];
Ts.retainedMalonylH.plastid.source='prescribed';
Ts.retainedMalonylH.plastid.deltaD_permil=-150;
Ts.retainedMalonylH.plastid.ratioNormalized=[];
Ts.retainedMalonylH.plastid.waterExchangeFraction=0.25;
Ts.retainedMalonylH.er.source='prescribed';
Ts.retainedMalonylH.er.deltaD_permil=-170;
Ts.retainedMalonylH.er.ratioNormalized=[];
Ts.retainedMalonylH.er.waterExchangeFraction=0.40;
Ts.cer3AldehydeH.mode='erNADPH';
Ppost=wax_hydrogen_source_pools_postmalonate_exactD(t,I,Ts);

details=struct;

S29=alkane_biosynthesis_h_stoichiometry_reaction_explicit(29,struct('C16',0,'C18',1));
H=S29.compounds(1).mixed.counts;
reported=[H.acetylStarterH,H.plastidRetainedMalonylH,H.erRetainedMalonylH, ...
    H.plastidKARNADPHHydrideH,H.plastidENRNADHHydrideH,H.plastidENRProtonH, ...
    H.erKCRNADPHHydrideH,H.erECRNADPHHydrideH,H.erECRProtonH,H.cer3AldehydeHydrideH];
expected=[3 8 6 8 8 8 6 6 6 1];
details.stoichiometryResidual=max(abs(reported-expected));
assert(details.stoichiometryResidual<1e-12 && abs(sum(reported)-60)<1e-12);

Tr=wax_reaction_incorporation_traits_null();

T2=Tr; T2.plastid.KAR.alphaHydride=0.90;
P2=wax_apply_reaction_incorporation_exactD(Ppost,T2);
r0=Ppost.classes.plastidNADPHHydrideH.ratio;
details.singleReactionResidual=max(abs(P2.classes.plastidKARNADPHHydrideH.ratio-0.90.*r0));
assert(details.singleReactionResidual<1e-12);

T3=Tr; T3.er.KCR.alphaHydride=0.92; T3.er.ECR.alphaHydride=0.97;
P3=wax_apply_reaction_incorporation_exactD(Ppost,T3);
rER=Ppost.classes.erNADPHHydrideH.ratio;
details.sharedSourceResidual=max([ ...
    max(abs(P3.classes.erKCRNADPHHydrideH.ratio-0.92.*rER)), ...
    max(abs(P3.classes.erECRNADPHHydrideH.ratio-0.97.*rER))]);
details.KCRECRSeparation=max(abs(P3.classes.erKCRNADPHHydrideH.deltaD-P3.classes.erECRNADPHHydrideH.deltaD));
assert(details.sharedSourceResidual<1e-12 && details.KCRECRSeparation>0);

T4=Tr;
T4.plastid.KAR.alphaHydride=0.96;
T4.plastid.ENR.alphaHydride=0.98;
T4.plastid.ENR.alphaProton=1.01;
T4.er.KCR.alphaHydride=0.94;
T4.er.ECR.alphaHydride=0.97;
T4.er.ECR.alphaProton=1.02;
T4.cer3.alphaAldehydeHydride=0.95;
P4=wax_apply_reaction_incorporation_exactD(Ppost,T4);
W4=wax_combine_hydrogen_classes_exactD(S29,P4);

names=S29.sourceClasses;
xExpected=zeros(1,numel(t));
for j=1:numel(names)
    nm=names{j};
    xExpected=xExpected+(H.(nm)./60).*P4.classes.(nm).DAtomFraction;
end
details.mixResidual=max(abs(W4.compounds(1).deltaD-xToDelta(xExpected)));
assert(details.mixResidual<1e-12);

details.upstreamResidual=max([ ...
    max(abs(P4.upstreamPostMalonatePools.classes.plastidNADPHHydrideH.deltaD-Ppost.classes.plastidNADPHHydrideH.deltaD)), ...
    max(abs(P4.upstreamPostMalonatePools.classes.erNADPHHydrideH.deltaD-Ppost.classes.erNADPHHydrideH.deltaD)), ...
    max(abs(P4.upstreamPostMalonatePools.classes.plastidRetainedMalonylH.deltaD-Ppost.classes.plastidRetainedMalonylH.deltaD))]);
assert(details.upstreamResidual<1e-12);

T6=Tr; T6.plastid.KAR.alphaHydride=linspace(0.90,1.00,numel(t));
P6=wax_apply_reaction_incorporation_exactD(Ppost,T6);
expectedRatio=T6.plastid.KAR.alphaHydride.*Ppost.classes.plastidNADPHHydrideH.ratio;
details.timeVaryingResidual=max(abs(P6.classes.plastidKARNADPHHydrideH.ratio-expectedRatio));
assert(details.timeVaryingResidual<1e-12);

P0=wax_apply_reaction_incorporation_exactD(Ppost,Tr);
Sall=alkane_biosynthesis_h_stoichiometry_reaction_explicit([27 29 31 33]);
W0=wax_combine_hydrogen_classes_exactD(Sall,P0);
frac=W0.sourceFractionTable{:,2:end};
details.sourceFractionResidual=max(abs(sum(frac,2)-1));
assert(details.sourceFractionResidual<1e-12);

R=struct('name','test_METCUT_reaction_incorporation_exactD','pass',true,'details',details);

fprintf('Stoichiometry residual:          %.3e H atoms\n',details.stoichiometryResidual);
fprintf('Single-reaction residual:        %.3e\n',details.singleReactionResidual);
fprintf('Independent mix residual:        %.3e per mil\n',details.mixResidual);
fprintf('Source-fraction residual:        %.3e\n',details.sourceFractionResidual);
fprintf('PASS\n');
fprintf('============================================================\n');
end
