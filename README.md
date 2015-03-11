%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CWM, 2014-06-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This is a stand-alone suite of MATLAB functions that implement the analytical solution to a nonlinear, hyperbolic PDE as described in [MacMinn, Szulczewski, and Juanes, JFM 2010](http://dx.doi.org/10.1017/S0022112010003319 "MacMinn et al., JFM 2010").

The user-facing functions are
  - hyp_params: Calculate characteristic scales and dimensionless parameters from aquifer and fluid properties etc.
  - hyp_crits: Calculate the critical times and positions
  - hyp_chars: Draw the full characteristics diagram
  - hyp_plume: Calculate the plume shape at desired times
  - hyp_eff: Calculate the efficiency factor (Eq. 4.1-4.2)

The functions in hyp_toolbox are helper functions. You should never need to interact with these directly.

The standard input parameters (all dimensionless) are:
  - M: Mobility ratio (below Eq. 2.9)
  - gamma: Capillary trapping number (below Eq. 2.2)
  - Nf: Flow number (Eqs. 2.10)
  - Ns: Slope number (Eqs. 2.10)
  - hmin: Optional. This is the minimum (dimensionless) plume thickness that you care about. For working with analytical solutions, usually set this to zero (the default). For comparing with experimental data or numerical solutions, usually use a small, nonzero number (ie, bead diameter / H).

The scripts fig_* will use hyp_crits, hyp_chars, hyp_plume, and hyp_eff to draw the figures in the paper.
  - Figure 4: fig_chars_inj
  - Figure 5: fig_fluxfns
  - Figure 6: fig_fluxes
  - Figure 7: fig_chars_collisions_shock
  - Figure 8: fig_chars_collisions_peak
  - Figure 11: fig_chars_migr_case1
  - Figure 12: fig_chars_migr_case2
  - Figure 13: fig_chars_migr_case3
  - Figure 14: fig_chars_migr_case4
  - Figure 15: fig_chars_migr_case6
  - Figure 16: fig_chars_migr_case0
  - Figure 17: fig_effs

Note that these scripts will all try to save figures to the subfolder ./figures, and will produce an error if no such folder exists.

To get started, see EXAMPLE.