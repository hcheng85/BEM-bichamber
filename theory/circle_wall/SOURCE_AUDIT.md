# Source audit: circular N=1 wall-only specialization

Audited production source (outside this publication checkout):
`archive_one_R_60/BEM_circular_parallel_instance.py`, SHA-256
`2a415d26ad83cf42a2f3cb0f4ba156b8dd6c11b32972ce8305499eebced22da4`.

- Lines 266--350 solve both circle/ray roots, reject negative roots, and retain
  the smallest admissible positive distance. The circle is centered at `(0,0)`.
- Lines 600--619 wrap `theta - theta_hit` into `(-pi,pi]` and use
  `Iw*sign(phi)/d`. Here `theta_hit` is the outward-normal angle of the hit
  point, so the sign in the analytical reduction is correct.
- Lines 614--617 multiply the wall torque by `1 + 1/d**5` when `d <= 0.1`,
  then clip it at `pi/0.001`. The latter is hard-coded, even if `para["dt"]`
  is changed.
- Lines 592--597 use `In*sqrt(dt)*N(0,1)`, hence the angular diffusion
  coefficient is `D_theta=In**2/2`.
- Lines 706--707 apply a second total-torque clip at `pi/para["dt"]`.
- Lines 724--749 perform an additional positional safety update: if the trial
  point is within radial clearance 0.1 and moved outward relative to the stored
  trial distance, translation is frozen, while deterministic angular motion
  and angular noise are retained.
- The polar chart is singular at `r=0`; the Cartesian source is not. At the
  exact center, a radial heading has zero signed incidence and therefore zero
  wall torque until translation breaks the degeneracy.

The last item changes the theory's applicability condition. It is not enough
that the orbit's forward-ray distance satisfies `d_star>0.1`; the source-exact
smooth orbit must also satisfy `R-r_star>0.1`. For `R=1`, this requires

```text
Iw > sqrt(1/0.9**2 - 1) = 0.4843221048.
```

Thus the formal `Iw=0.25` orbit lies inside the production sticky layer and is
an intentional negative-control case, not a valid test of the smooth branch.
