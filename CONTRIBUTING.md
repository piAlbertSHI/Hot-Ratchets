# Contributing

Run the numerical checks from the repository root before proposing a change.

```sh
matlab -batch "addpath('tests'); run_tests"
```

For a constitutive change, include the material card, complete control history, units, and a small reproducible test. State whether the comparison checks a numerical identity, a synthetic simulation, or an experimental result. Link to published sources for equations and material parameters.

Use the `matlab/+sma` package for the solver, `examples` for runnable demonstrations, and `tests` for contract checks. Use English filenames and documentation. Keep generated files in `outputs/` and cite source documents by reference.
