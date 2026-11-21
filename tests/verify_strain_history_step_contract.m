function report = verify_strain_history_step_contract()
%VERIFY_STRAIN_HISTORY_STEP_CONTRACT Contract and smoke test for explicit strain history.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab'));
p = sma.params('niti_zhang2022');
s = sma.initial_state_1d(p);

% Synthetic elastic-scale targets exercise the interface only; they are not a
% Yuan protocol and are not compared with observations.
targets = [1.0e-4; 2.0e-4];
durations = [0.01; 0.01];
[trial,aux] = sma.strain_history_step_1d(p,s,targets,durations,p.Tamb);

errors = struct();
errors.endpoint_target = abs(trial.eps-targets(end));
errors.receipt_endpoint = abs(aux.endpoint_control_residual);
errors.duration = abs(aux.total_duration_s-sum(durations));
errors.control_mode = double(~strcmp(aux.control_mode,'explicit_strain_history'));
errors.target_count = double(numel(aux.target_strain)~=numel(targets));
errors.duration_count = double(numel(aux.step_duration_s)~=numel(durations));

invalid_empty = false;
try
    sma.strain_history_step_1d(p,s,[],[],p.Tamb);
catch ME
    invalid_empty = strcmp(ME.identifier,'sma:strain_history_step_1d:History');
end
errors.empty_history_rejected = double(~invalid_empty);

limits = struct('endpoint_target',1.0e-8, ...
    'receipt_endpoint',1.0e-8,'duration',1.0e-14, ...
    'control_mode',0,'target_count',0,'duration_count',0, ...
    'empty_history_rejected',0);
report = struct('status','PASS_EXPLICIT_STRAIN_HISTORY_CONTRACT', ...
    'synthetic_only',true,'targets',targets,'durations_s',durations, ...
    'trial',trial,'aux',aux,'errors',errors,'limits',limits, ...
    'experimental_acceptance',false,'production_protocol_modified',false);
out = fullfile(root,'outputs','residual_audit');
if ~exist(out,'dir'), mkdir(out); end
save(fullfile(out,'strain_history_step_contract.mat'),'report');
names = string(fieldnames(errors)); vals = zeros(numel(names),1);
for k=1:numel(names), vals(k)=errors.(char(names(k))); end
writetable(table(names,vals,'VariableNames',{'Check','Error'}), ...
    fullfile(out,'strain_history_step_contract.csv'));
assert(errors.endpoint_target<=limits.endpoint_target && ...
    errors.receipt_endpoint<=limits.receipt_endpoint && ...
    errors.duration<=limits.duration && all(struct2array(errors)<= ...
    struct2array(limits)), ...
    'Explicit strain-history contract failed.');
disp(report.status);
end
