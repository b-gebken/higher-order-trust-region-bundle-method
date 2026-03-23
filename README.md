<h1>Higher-order trust-region bundle method</h1>

An implementation of a trust-region bundle method using higher-order cutting-plane models for the solution of nonsmooth, nonconvex optimization problems. It was used for the numerical experiments in [GU2026a] and [GU2026b].

[GU2026a] Gebken, Ulbrich (2026): Superlinear convergence in nonsmooth optimization via higher-order cutting-plane models (To be submitted)<br/>
[GU2026b] Gebken, Ulbrich (2026): Enclosing minima in nonsmooth optimization via trust regions of higher-order cutting-plane models (To be submitted)<br/>

The correspondence between paper examples and folders is as follows:<br/>
| Paper    | Folder |
| -------- | ------- |
| [GU2026a]  | GU2026a_Local_experiments  |
| Fig. 1 | Sketch_model |
| Ex. 6.1/Fig. 2(a) | Example_1d_symbolic |
| Ex. 6.2/Fig. 2(b)/Fig. 3 | Example_LW2019_85 |
| Ex. 6.3/Fig. 4 | Example_LW2019_eigval |
| Ex. 6.4/Fig. 5 | Example_halfhalf |
| Ex. 6.5/Fig. 6 | Example_max_root |

| Paper    | Folder |
| -------- | ------- |
| [GU2026b]  | GU2026b_Global_experiments  |
| Ex. 4.1/Fig. 1 | Example_polynomial_growth_not_sufficient |
| Ex. 5.1/Fig. 2 | Example_LW2019_84_global |
| Ex. 5.2/Fig. 3 | Example_LW2019_85_global |
| Ex. 5.3/Fig. 4(a) | Example_LW2019_eigval_global |
| Rem. 5.1 | Remark_toy_example |
| Ex. 5.4/Fig. 4(b) | Example_LW2019_eigval_local |
| Ex. 5.5/Fig. 4(c) | Example_LW2019_84_local |

The trust-region subproblems are solved via the Matlab interface _mexIPOPT_ (https://github.com/ebertolazzi/mexIPOPT) of _IPOPT_ (https://github.com/coin-or/Ipopt). The runtimes reported in the papers were achieved with a slight modification to _mexIPOPT_, where the check whether Matlab or Octave is used is disabled. See "Algorithms/solve_subproblem_IPOPT.m" for details.

The figures are saved via _export_fig_ (https://github.com/altmany/export_fig).

<h1>Acknowledgements</h1>
This research was funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Projektnummer 545166481.
