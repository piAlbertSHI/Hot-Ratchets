# Model notes

## Material updates

`sma.params` constructs a material card. `sma.initial_state_1d` and `sma.initial_state_3d` construct committed states for that card. A constitutive update takes the committed state and a control target and returns the next state and an auxiliary structure. Advance the committed state after checking the returned residuals for your application. Keep the material card and the time sequence with a saved history.

| Quantity | Convention |
| --- | --- |
| Stress | MPa. Tensor stress control uses a symmetric Kirchhoff stress tensor. |
| Strain | Dimensionless logarithmic strain in the uniaxial solver. The Cu–Al–Mn engineering-strain-rate interface converts its target with `log1p`. |
| Temperature and time | K and s. |
| Specific energy and source rate | J/kg and W/kg in the material heat balance. |
| Cylindrical band geometry | Radius in m and dimensionless length-to-radius ratio. The band-energy gradient is returned in MPa m². |

The material state includes martensite fraction, recoverable transformation strain, transformation-induced plastic strain (TRIP), austenite/R-phase plastic strain (ARP), and history variables for cyclic hardening and retained martensite. States with an explicit R phase carry a separate phase fraction. The restart helpers `sma.pack_restart_1d_v1` and `sma.unpack_restart_1d_v1` bind the state payload to its parameter card and mechanism settings.

## Implementations

- **NiTi** uses a Zhang-based multimechanism material card. Standard uniaxial and tensor updates provide the public entry points. The `candidate_*` routines include source-level phase, plastic history, and thermal residual experiments. [Zhang et al. (2022)](https://doi.org/10.1016/j.ijsolstr.2022.111638)
- **Cu–Al–Mn** uses reduced material cards for annealed and cold-drawn microwires. The reduced Yuan provider carries a normalized reaction extent and explicit strain-rate control. The example exercises the annealed card. [Yuan et al. (2020)](https://doi.org/10.1016/j.ijrefrig.2020.02.020)
- **Sun band extension** evaluates the elastic shape factor and the derivative of band energy with respect to length for an axisymmetric cylindrical inclusion. These kernels support band-energy studies. [Zhong, Sun, and Tong (2000)](https://doi.org/10.1016/S0020-7683(99)00269-3) · [Sun and Zhong (2000)](https://doi.org/10.1016/S0749-6419(00)00006-1)

## Verification scope

The public checks exercise state compatibility, explicit control histories, composition-dependent transformation properties, a coupled stress-control transaction, the reduced Cu–Al–Mn step, and the band-energy derivative. Each example uses a short control history with its inputs defined in the script. Full experimental-cycle validation and the spatial coupling of band motion to material kinetics remain research work.
