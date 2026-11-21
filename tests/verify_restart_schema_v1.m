function report=verify_restart_schema_v1()
%VERIFY_RESTART_SCHEMA_V1 Contract test for lossless, self-describing restarts.
solver_root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(solver_root, 'matlab'));
p=sma.params('niti_zhang2022');
s=sma.initial_state_1d(p);
s.eps_tr=0.0123; s.xi=0.17; s.eps=0.04; s.xi_r=0.11;
if p.enable_explicit_r_phase_state, s.xi_R=0.62; end
packet=sma.pack_restart_1d_v1(p,s);
restored=sma.unpack_restart_1d_v1(p,packet);
payload_error=max(abs(sma.pack_state_1d(p,restored)-packet.state_vector));
assert(payload_error==0,'Restart payload changed after round trip.');

% Exercise every optional payload branch with nonzero history so that the
% contract is not only proven for the 15-value legacy prefix.
p2=p; p2.enable_explicit_r_phase_state=true;
p2.Y_activity_population_mix=0.5; p2.zhang_direct_Yp_coupling=1.0;
s2=sma.initial_state_1d(p2); s2.xi=0.21; s2.xi_R=0.58;
s2.Y_population_activity(:)=0.07; s2.zhang_xi_Mr=0.02;
s2.zhang_rho_M=1.2e13; s2.zhang_rho_AR=1.4e13;
s2.zhang_D_M=2.5e12; s2.zhang_D_AR=1.1e13;
s2.zhang_Yp=0.3; s2.zhang_Omega=1.1; s2.zhang_Lambda_M=0.04;
packet2=sma.pack_restart_1d_v1(p2,s2);
restored2=sma.unpack_restart_1d_v1(p2,packet2);
optional_payload_error=max(abs(sma.pack_state_1d(p2,restored2)-packet2.state_vector));
assert(optional_payload_error==0,'Optional restart payload changed after round trip.');

errors=struct();
q=packet; q.model_case_id=999;
errors.card=expect_rejection(@()sma.unpack_restart_1d_v1(p,q),'CardMismatch');
q=packet; q.parameter_signature.EA=q.parameter_signature.EA+1;
errors.parameters=expect_rejection(@()sma.unpack_restart_1d_v1(p,q),'ParameterMismatch');
q=packet; q.feature_flags.enable_explicit_r_phase_state=~q.feature_flags.enable_explicit_r_phase_state;
errors.flags=expect_rejection(@()sma.unpack_restart_1d_v1(p,q),'FeatureMismatch');
q=packet; q.state_width=q.state_width+1;
errors.width=expect_rejection(@()sma.unpack_restart_1d_v1(p,q),'Width');
q=rmfield(packet,'semantic_contract');
errors.metadata=expect_rejection(@()sma.unpack_restart_1d_v1(p,q),'MissingMetadata');

report=struct('status','PASS_RESTART_SCHEMA_V1_LOSSLESS', ...
    'payload_error',payload_error,'state_width',packet.state_width, ...
    'optional_payload_error',optional_payload_error, ...
    'optional_state_width',packet2.state_width, ...
    'card_rejection',errors.card,'feature_rejection',errors.flags, ...
    'parameter_rejection',errors.parameters, ...
    'width_rejection',errors.width,'metadata_rejection',errors.metadata, ...
    'production_forward_connected',false,'experimental_data_read',false);
disp(report);
end

function msg=expect_rejection(fn,needle)
try
    fn();
    error('Expected rejection did not occur.');
catch e
    if contains(e.identifier,needle) || contains(e.message,needle)
        msg=e.message;
    else
        rethrow(e);
    end
end
end
