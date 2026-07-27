function value = update_training_activity(p,state,trial,branch)
%SMA.UPDATE_TRAINING_ACTIVITY Commit normalized defect activity per loop.

value = state.zeta_p;
if ~p.enable_mt_dissipation_evolution || branch>=0 || ...
        state.xi<=state.xi_res+p.mt_training_completion_tol || ...
        trial.xi>trial.xi_res+p.mt_training_completion_tol
    return;
end
peak_fraction = max(state.xi_r,state.xi);
normalized_peak = min(peak_fraction/p.mt_training_reference_fraction,1.0);
value = value+p.mt_training_cycle_increment*normalized_peak;
end
