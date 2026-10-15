
function C = run_full_chain_upstream_I02_candidate(Mreference,baseSample,I01,I02)
%RUN_FULL_CHAIN_UPSTREAM_I02_CANDIDATE
% Exact current upstream chain through intracellular water using the GSA
% production wrapper. Hydraulics are reused from Mreference.

sample=baseSample;
sample=set_gsa_sample_parameter(sample,'I01',I01);
sample=set_gsa_sample_parameter(sample,'I02',I02);

Ustate=gsa_build_upstream_state(Mreference,sample);
Mup=extract_model_from_upstream_state(Ustate);

C=struct;
C.I01=I01;
C.I02=I02;
C.sample=sample;
C.Ustate=Ustate;
C.model=Mup;
C.intracellular=Mup.intracellular;

end
