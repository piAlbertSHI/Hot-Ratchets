function state = commit_zhang_direct_history(state,history)
%SMA.COMMIT_ZHANG_DIRECT_HISTORY Store a converged direct history trial.

state.zhang_xi_Mr=history.xi_Mr;
state.zhang_rho_M=history.rho_M; state.zhang_rho_AR=history.rho_AR;
state.zhang_D_M=history.D_M; state.zhang_D_AR=history.D_AR;
state.zhang_DeltaEtaP=history.DeltaEtaP;
state.zhang_Omega=history.Omega; state.zhang_Yp=history.Yp;
state.zhang_Lambda_M=history.Lambda_M;
end
