# Method and source audit

## Huang, Ling, and Kanso reference method

Supplementary Section S.4 assumes that the reduced variable is a stationary,
continuous Markov process.  It bins the reduced coordinate, collects short
trajectory segments beginning in each bin, and obtains drift and diffusion
from linear slopes of the conditional ensemble mean and variance.  With the
convention used here,

\[
\partial_t p=-\partial_q(Vp)+\partial_q^2(Dp),\qquad
D=\frac{1}{2}\lim_{\tau\to0}\frac{\operatorname{Var}(\Delta q\mid q)}{\tau}.
\]

The zero-current effective potential is

\[
\Phi(q)=\log D(q)-\int^q V(s)/D(s)\,ds+C.
\]

The paper's `K` is continuous by construction and follows a well-populated
one-dimensional manifold in its original two-dimensional order-parameter
space.  Neither property is assumed for chamber occupancy in this audit.

## Problems in the archived `FK_potential` scripts

### `effective_potential.m`

- It allocates 1,000 states on `[0,1]` but writes and integrates `phi` through
  indices 1,987--2,000.  Those indices do not refer to the allocated drift and
  diffusion arrays.
- Its hard-coded integration spacing (`1e-3`) is not derived from the actual
  grid object.
- Empty bins remain numerical zeros.  Dividing `V./D` therefore creates
  undefined values, while hard-coded index exclusions conceal rather than
  report the support problem.
- Only the first two conditional samples are used even though ten are
  collected.  This does not implement the regression procedure described in
  the supplement.

### `effective_potential_linear_regression.m`

- It declares `dt=1e-2`, but constructs its regression time axis in `1e-3`
  increments.
- `linspace(0,1,2001)` has spacing `5e-4`, not the commented `1e-3`; the
  tolerance and later potential integration use yet another implied spacing.
- It assumes an externally existing `Collector`, with no check that its grid
  matches the newly allocated 2,000-state grid.
- Unsupported regions are again selected through fixed indices rather than
  measured sample support.

### `effective_potential_pdf.m`

- Construction of `Collector` is commented out, so the script depends on
  mutable MATLAB workspace state.
- Column 10 of a ten-point segment beginning at time zero is at `9*dt`, but the
  conditional moments are divided by `1e-2`.
- MATLAB histogram values with `Normalization='probability'` are probability
  masses.  Applying `trapz` to them without a consistent grid measure is not a
  correct expectation; a weighted sum is required.
- The valid integration interval is hard-coded, and sample counts and
  uncertainties are not retained.

None of these archived scripts performs an independent-realization holdout,
Chapman--Kolmogorov test, stationary-current check, lag-stability audit, or
first-passage validation.

## Corrections in the standalone gap-50 analysis

- The primary grid is exactly the 101 occupancy states separated by `2/N=0.02`.
- Physical lag is always calculated from the stored `dt=0.001`.
- Conditional increment moments and short-window regressions are both retained.
- Every coefficient carries its origin count; unsupported or non-positive-`D`
  states are excluded, never filled with artificial zeros.
- Potential integration uses the measured coordinate spacing and only a
  contiguous supported component.
- Seeds 0--14 train the model; seeds 15--19 validate stationary statistics,
  hidden-current closure, and first-passage times.
- A continuous coordinate uses `ell=0.5`, the rounded divider-tip radius, with
  `ell=0.25` and `1.0` as sensitivity controls.  It is not substituted for the
  discrete coordinate merely to force a smooth potential.
- The neck current is the advective derivative of the `ell=0.5` smoothed
  occupancy,

  \[
  J_{\rm neck}=N^{-1}\sum_i v_{x,i}\,\operatorname{sech}^2(x_i/\ell)/\ell.
  \]

- The backward equation is solved with state-dependent diffusion; no
  constant-diffusion Kramers formula is used.

## Dataset source verification

The 20 raw files have `state=(3,100,1000002)` and
`rdot=(2,100,1000002)`.  Stored parameters are `R=7.1589105316`,
`gap=7.1589105316`, `gap/R=1.0`, `gap/(2R)=0.5`, `N=100`, `Ip=4`, `In=0.5`, `Iw=0.94`,
`If=0.01`, `dt=0.001`, and `T=1000`.  In
`BEM_bichamber_parallel_instance.py`, `GetHydro` and `GetHydroBEM` are both
called and their translational and rotational contributions enter `rDot` and
`thetaDot`.  Hydrodynamics is therefore source- and metadata-verified rather
than inferred from filenames.

The task described this as `gap/R=0.5`.  Under the stored parameter definitions,
that literal ratio would require an `R50/gap25` ensemble; none was found.  The
available `R50/gap50` ensemble corresponds to 0.5 only under diameter
normalization.  This analysis does not silently substitute one convention for
the other.
