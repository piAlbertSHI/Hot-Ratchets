function [trial,receipt] = candidate_yuan_strain_rate_event_step_1d( ...
    p,state,signed_rate_1_per_s,event_stress_MPa,Tamb,dt,event_tolerance_MPa)
%SMA.CANDIDATE_YUAN_STRAIN_RATE_EVENT_STEP_1D
% Advance a signed engineering-strain-rate step and stop at a stress event.
%
% A crossed event is located by bisection in elapsed time.  Every probe
% starts from the same committed STATE; the stress endpoint is never assigned
% to the state and no observation is used to repair the constitutive result.

if nargin<7 || isempty(event_tolerance_MPa), event_tolerance_MPa=1e-6; end
assert(isfinite(signed_rate_1_per_s) && ...
    any(signed_rate_1_per_s==[0.15,-0.15]), ...
    'sma:candidate_yuan_strain_rate_event_step_1d:Rate', ...
    'Yuan direct-eCE candidate rate must be signed +/-0.15 1/s.');
assert(isfinite(event_stress_MPa) && isfinite(event_tolerance_MPa) && ...
    event_tolerance_MPa>0 && isfinite(dt) && dt>0, ...
    'sma:candidate_yuan_strain_rate_event_step_1d:Event', ...
    'Event target, tolerance and positive time step are required.');

control = struct('mode','engineering_strain_rate', ...
    'signed_rate_1_per_s',signed_rate_1_per_s);
old_residual = state.sigma-event_stress_MPa;
[full_trial,full_receipt] = sma.candidate_yuan_reduced_material_step_1d( ...
    p,state,control,Tamb,dt);
full_residual = full_trial.sigma-event_stress_MPa;
crossed = old_residual*full_residual<=0;
iterations = 0;
accepted_dt = dt;
trial = full_trial;
step_receipt = full_receipt;
if crossed && abs(full_residual)>event_tolerance_MPa && ...
        abs(old_residual)>event_tolerance_MPa
    lo = 0;
    hi = dt;
    rlo = old_residual;
    for iterations=1:80
        mid = 0.5*(lo+hi);
        [candidate,candidate_receipt] = ...
            sma.candidate_yuan_reduced_material_step_1d( ...
            p,state,control,Tamb,mid);
        residual = candidate.sigma-event_stress_MPa;
        trial = candidate;
        step_receipt = candidate_receipt;
        accepted_dt = mid;
        if abs(residual)<=event_tolerance_MPa
            break
        end
        if rlo*residual<=0
            hi = mid;
        else
            lo = mid;
            rlo = residual;
        end
    end
end
event_residual = trial.sigma-event_stress_MPa;
event_reached = crossed && abs(event_residual)<=event_tolerance_MPa;
receipt = step_receipt;
receipt.schema_id = 'yuan_strain_rate_event_step_v1';
receipt.event_stress_MPa = event_stress_MPa;
receipt.event_tolerance_MPa = event_tolerance_MPa;
receipt.event_residual_MPa = event_residual;
receipt.event_crossed = crossed;
receipt.event_reached = event_reached;
receipt.requested_dt_s = dt;
receipt.accepted_dt_s = accepted_dt;
receipt.event_iterations = iterations;
receipt.endpoint_assignment_used = false;
receipt.observation_accessed = false;
end
