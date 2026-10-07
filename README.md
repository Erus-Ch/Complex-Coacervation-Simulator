# LLPS Complex Coacervation Simulator

A standalone, local Voorn–Overbeek (VO) simulator for exploring complex-coacervation phase behavior. The GUI includes a symmetric binodal with tie-lines, an asymmetric ion-resolved 2D spinodal map, and an interactive 3D spinodal surface over charge ratio, total polyion molecule concentration, and added-salt concentration.

Equations use GitHub Markdown math delimiters: `$...$` for inline notation and `$$...$$` for display equations. See [GitHub’s mathematical-expression guide](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).

## Quick start

Requires Python 3.10 or newer. From this directory:

~~~bash
python3 -m venv .venv
source .venv/bin/activate                 # Windows: .venv\\Scripts\\activate
python -m pip install --upgrade pip
python -m pip install -e .
vo-sim                                    # opens http://127.0.0.1:8770
~~~

On macOS, double-click `launch_vo_simulator.command` to start the local GUI and open it in a browser. The simulator binds to loopback only; calculations run locally. No molecular database, experiment table, fitted model, or image-analysis component is included.

## What is included

- **Symmetric VO binodal:** dilute and dense coexistence branches, tie-lines, a feed lever-rule readout, an approximate critical marker, and a reduced-free-energy/common-tangent view.
- **Asymmetric 2D spinodal:** at fixed charge-equivalent feed ratio, scans total polyion molecule concentration against selected-salt formula concentration; shows counterions, salt ions, ionic strength, electroneutrality residual, and the smallest constrained-Hessian eigenvalue.
- **Higher-dimensional 3D spinodal:** axes are exactly $R_q$ (dimensionless), $c_P+c_Q$ (M of polyion molecules), and $c_s$ (M of selected salt formula units). Each axis has independently editable bounds and physical units. Drag rotates; Shift-drag pans; wheel/buttons zoom; reset restores the camera. A selectable 2D slice is shown below the surface.
- **Salt presets:** NaCl, KCl, MgCl₂, CaCl₂, Na₂SO₄, K₂SO₄, and MgSO₄. The interface reports the salt's conventional ionic-strength factor separately from formula-unit concentration.
- **Adjustable model parameters:** chain segment counts, charged-site fractions, lattice spacing, Bjerrum length, and optional polymer–water interaction parameter.
- **Local interactive plots:** phase diagrams and energy profiles support zoom/pan and CSV/SVG export where available.

## Model definitions

The interface has two related VO calculations. The symmetric page computes a global coexistence binodal within a symmetric reduction. The asymmetric pages resolve counterions and salt ions, allow non-stoichiometric feeds, and compute local stability limits. A spinodal is not a binodal: it marks loss of local stability and does not by itself give the full two-phase coexistence envelope or tie-lines.

### Symmetric reduction and binodal

Let $p$ be the total volume fraction of positive and negative polyion segments, split equally as $\phi_P=\phi_Q=p/2$. Let $n$ be the volume fraction of **each** mobile-ion sign, including counterions and added monovalent 1:1 salt, and let $\phi_w=1-p-2n$ be the solvent fraction. A chain has $N$ lattice segments; $q$ is the fraction of charged segments. The reduced Helmholtz free-energy density is

$$
\begin{aligned}
\tilde{f}
&= \frac{f a^3}{k_{\mathrm B}T} \\
&= \frac{p}{N}\ln\left(\frac{p}{2}\right)
   + 2n\ln n + \phi_w\ln\phi_w \\
&\quad - \alpha(qp+2n)^{3/2} + \chi_{pw}p\phi_w, \\
\alpha &= \frac{2\sqrt{\pi}}{3}\left(\frac{l_B}{a}\right)^{3/2}, \\
\kappa^2 &= \frac{4\pi l_B}{a^3}(qp+2n).
\end{aligned}
$$

Here $a$ is the equal-volume lattice length and $l_B$ is the Bjerrum length. The first three terms are ideal polymer, mobile-ion, and solvent mixing entropies. The negative term is the Debye–Hückel (DH) electrostatic correlation contribution. $\chi_{pw}$ is an optional short-range polymer–water Flory term; $\chi_{pw}=0$ retains the main VO mixing and DH terms. The screening sum counts charged **sites** ($qp$) and mobile monovalent ions ($2n$), rather than treating a whole polyion as a single point charge.

For a symmetric feed with total polymer-segment molarity $c_p$, each sign of mobile ion has

$$
c_{\mathrm{ion}} = c_s + \frac{q c_p}{2},
\qquad
c_{\mathrm{lattice}} = \frac{1}{0.602214076\,a^3}\,\mathrm{M}.
$$

where $c_s$ is added 1:1 salt formula molarity and $a$ is in nm. Thus $p c_{\mathrm{lattice}}$ is total polymer-segment molarity and $n c_{\mathrm{lattice}}$ is mobile-ion molarity of one sign. The symmetric construction fixes equal polyion segment concentrations and charge-equivalent ratio $R_q=1$; it does not represent an arbitrary feed ratio.

**Common-tangent solver.** At a chosen mobile-ion chemical potential $\mu_n$, the code first minimizes over the mobile-ion fraction:

$$
g(p;\mu_n) = \min_n\left[\tilde{f}(p,n)-\mu_n n\right].
$$

It finds the candidate non-convex interval using a monotone-chain lower convex hull on a polymer-fraction grid. The mobile-ion chemical-potential equation is solved at all grid points by 45 vectorized bisection iterations. SciPy `optimize.root` then refines the dilute and dense endpoints in log-polymer and logit-dense coordinates by enforcing equality of polymer chemical potential and semi-grand tangent intercept. The resulting common tangent defines the coexistence pair; the GUI connects pairs with tie-lines. The high-salt end is bracketed/refined with 15 bisection steps. The reported critical marker is the approximate midpoint of the last converging pair, not a dedicated critical-point root solve. The energy-profile plot shows $g-g_{\mathrm{tangent}}$ for the selected pair.

### Asymmetric, ion-resolved extension

The asymmetric pages use molecule molarities $c_P$ and $c_Q$, not segment molarities, as the feed concentration coordinates. Let $N_P$, $N_Q$ be segment counts and $q_P$, $q_Q$ charged-site fractions. The GUI's charge-equivalent feed ratio and total polyion molecule concentration are

$$
R_q = \frac{q_P N_P c_P}{q_Q N_Q c_Q},
\qquad
c_{\mathrm{tot}} = c_P+c_Q,
\qquad
\phi_P = \frac{N_Pc_P}{c_{\mathrm{lattice}}},
\quad
\phi_Q = \frac{N_Qc_Q}{c_{\mathrm{lattice}}}.
$$

Any positive $R_q$ is allowed within the displayed scan range; global electroneutrality does not require $R_q=1$. Each polyion-derived counterion is monovalent in this implementation. Feed counterion fractions are $\phi_X=q_P\phi_P$ for $X^-$ and $\phi_Y=q_Q\phi_Q$ for $Y^+$. For a selected salt with formula concentration $c_s$, ionic charges $z_+$, $z_-$, and formula stoichiometries $\nu_+$, $\nu_-$, its ion fractions are

$$
\phi_{S^+} = \frac{\nu_+c_s}{c_{\mathrm{lattice}}},
\qquad
\phi_{S^-} = \frac{\nu_-c_s}{c_{\mathrm{lattice}}}.
$$

and charge neutrality is maintained by the salt stoichiometry $\nu_+z_+ + \nu_-z_-=0$. The model's total bulk charge constraint is

$$
q_P\phi_P-q_Q\phi_Q-\phi_X+\phi_Y
+z_+\phi_{S^+}+z_-\phi_{S^-}=0.
$$

The asymmetric reduced free-energy density is

$$
\begin{aligned}
\tilde{f}
&= \frac{\phi_P}{N_P}\ln\phi_P
 + \frac{\phi_Q}{N_Q}\ln\phi_Q
 + \sum_{i\in\{X,Y,S^+,S^-\}}\phi_i\ln\phi_i
 + \phi_w\ln\phi_w \\
&\quad - \alpha I_{\mathrm{DH}}^{3/2}
 + \chi_{pw}(\phi_P+\phi_Q)\phi_w, \\
\phi_w &= 1-\phi_P-\phi_Q-\phi_X-\phi_Y-\phi_{S^+}-\phi_{S^-}, \\
I_{\mathrm{DH}}
&= q_P\phi_P+q_Q\phi_Q+\phi_X+\phi_Y
 + z_+^2\phi_{S^+}+z_-^2\phi_{S^-}, \\
\alpha &= \frac{2\sqrt{\pi}}{3}\left(\frac{l_B}{a}\right)^{3/2}.
\end{aligned}
$$

The DH screening variable includes charged polyion sites, their monovalent counterions, and each explicit added-salt ion weighted by charge squared. It is distinct from conventional salt ionic strength. For the selected salt,

$$
I_s = \frac{1}{2}
\left(\nu_+z_+^2+\nu_-z_-^2\right)c_s.
$$

so $I_s/c_s$ is 1 for 1:1 salts, 3 for 2:1 and 1:2 salts, and 4 for 2:2 salts. The plot axis remains formula-unit concentration $c_s$ in M; the GUI labels/reports $I_s$ separately.

**Asymmetric spinodal solver.** At each feed point the code maps molecule concentrations, selected-salt ions, and polyion counterions into lattice fractions and checks the packing constraint $\phi_w>0$. It builds the free-energy Hessian analytically in independent composition coordinates. A constraint basis eliminates one counterion direction to enforce local electroneutrality; the Hessian is projected into this neutral subspace, and its smallest eigenvalue is the local-stability diagnostic. The ideal-mixing, DH, and polymer–water terms each contribute their analytic second derivatives. Negative values indicate local instability; the spinodal is the zero-eigenvalue boundary. In 3D, the implementation scans salt concentration at each $(R_q,c_{\mathrm{tot}})$ pair, detects a negative-to-nonnegative crossing, and linearly interpolates the salt threshold. The surface is therefore a finite-grid local-stability estimate, not a global phase-equilibrium calculation.

The default 3D scan uses 19 $R_q$ × 24 $c_{\mathrm{tot}}$ × 61 $c_s$ grid points. GUI range inputs control the physical minimum and maximum of every axis; changing those ranges recalculates the surface over those coordinates. Narrow features can be missed by a finite grid. Increase scan density in the UI/API only when needed, since cost rises with the grid size.

## Parameters and defaults

All concentrations shown below are molar (`M`).

| Parameter | Default | GUI/API range | Meaning |
| --- | ---: | ---: | --- |
| Symmetric $N$ | 80 | 1–400 | Lattice segments per polyion |
| Symmetric $q$ | 0.75 | 0.10–1.00 | Charged fraction of polymer segments |
| $N_P$, $N_Q$ | 80, 80 | 1–400 each | Polycation/polyanion segment count |
| $q_P$, $q_Q$ | 0.75, 0.75 | 0.05–1.00 each | Charged-site fraction for each polyion |
| $a$ | 0.60 nm | 0.56–1.20 nm | Equal-volume lattice spacing |
| $l_B$ | 0.70 nm | 0.20–1.00 nm | Bjerrum length; code also requires $l_B/a\leq 1.8$ |
| $\chi_{pw}$ | 0 | -0.40–0.50 | Optional polymer–water interaction parameter |
| Asymmetric $R_q$ | 0.50 in 2D page | 0.10–5.00 | Feed charge-equivalent ratio; 3D default scan is 0.25–4.00 |
| 2D $c_{\mathrm{tot}}$ maximum | 0.025 M | UI up to 0.5 M | Maximum total polyion molecule concentration |
| 2D $c_s$ maximum | 1.00 M | UI up to 5 M | Maximum selected-salt formula concentration |
| 3D $c_{\mathrm{tot}}$ range | 0.0001–0.025 M | independently editable | Total polyion molecule concentration axis |
| 3D $c_s$ range | 0–1.00 M | independently editable | Selected-salt formula concentration axis |

The composition builder rejects negative/non-finite concentrations, $R_q$ outside 0.1–5, and feeds that violate lattice packing. Some combinations of nominally valid slider values may therefore return a packing-limit error at high total concentration or salt.

## What this model leaves out

This is a mean-field lattice/DH simulator. It assumes equal-volume lattice sites and ideal mixing for small-ion species. The asymmetric calculation represents the polymer-derived counterions as monovalent and computes a local-stability Hessian at feed composition. Neither calculation includes ion-specific activities, ion pairing/binding, counterion-release equilibria, finite-size ion packing beyond the lattice solvent constraint, unequal species volumes, polymer conformational entropy beyond the Flory chain term, explicit spacer connectivity, or solvent-dependent dielectric response. The symmetric binodal and asymmetric spinodal are different reductions/outputs; compare them as complementary views, not interchangeable boundaries. Treat predictions as theory-model behavior under the stated assumptions, especially for concentrated solutions and multivalent small molecules.

The symmetric formulation follows the VO phase-diagram discussion in Sing and Perry, *Soft Matter* 16, 2885–2914 (2020), [doi:10.1039/D0SM00001A](https://doi.org/10.1039/D0SM00001A). The general VO formulation originates with Overbeek and Voorn, *Journal of Cellular and Comparative Physiology* 49, 7–26 (1957), [doi:10.1002/jcp.1030490404](https://doi.org/10.1002/jcp.1030490404).

## Development and package commands

~~~bash
python -m pip install -e .
vo-sim --port 8770
vo-sim --port 8770 --no-browser
python -m pip wheel --no-deps . --wheel-dir dist
~~~

No project-wide source license has been selected; the bundled MathJax assets retain their separate Apache-2.0 license, included with the assets.

The local HTTP API is served on `127.0.0.1` and exposes the same calculations used by the GUI: `/api/simulate`, `/api/profile`, `/api/asymmetric/2d`, and `/api/asymmetric/high-dimensional`. MathJax assets are bundled locally so formulas render without a CDN connection.
