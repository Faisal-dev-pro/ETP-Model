# Paper C — Complete Draft Package

> **Working title:** A Coupled Electro-Thermal-Pressure Model for Prediction of Internal
> Gas Pressure Evolution During Thermal Runaway in 21700 NMC Lithium-Ion Cells
>
> **Target journal:** Journal of Power Sources Advances (Elsevier, open access)
>
> **Validation data:** Gulsoy et al. (2024/2025), WMG open-access dataset
>
> **Base model:** Ostanek (2025) Simulink implementation of Hatchard–Kim 0D lumped-body TR model


---

## PART 1 — FULL PAPER OUTLINE

### 1. Introduction
- LIB thermal runaway safety problem and the need for predictive models
- Internal gas pressure as an under-modelled but critical observable
- Gap: most Hatchard–Kim implementations predict temperature only; pressure sub-models
  are rare and typically untethered from reaction kinetics
- Contribution statement (see below)

### 2. Literature Review

#### 2.1 Hatchard–Kim Reaction Kinetics
- Hatchard et al. (1999); Kim et al. (2007) — canonical four-reaction abuse model
- Ostanek (2025) — open-source Simulink re-implementation, 0D lumped body, 18650 focus
- Parhizi, Ostanek & Jeevarajan (2023) — numerical stiffness alleviation for the same framework

#### 2.2 Internal Pressure Modelling in TR
- Coman et al. (2017) — early coupled electro-thermal + gas generation model (18650 NMC)
- Kong et al. (2022), Mao et al. (2022) — venting and jet-fire dynamics, chamber-pressure validation
- Chen, Gulsoy, Barai et al. (2025) — lumped TR model with vapour–liquid equilibrium pressure
  sub-model for 21700 NMC, validated against Gulsoy et al. dataset [KEY COMPARATOR]

#### 2.3 Experimental Datasets
- Gulsoy, Chen, Briggs, Vincent, Sansom & Marco (2024) — first simultaneous in-situ
  internal T + P measurement in 21700 cells (J. Power Sources, 617, 235147)
- Gulsoy, Briggs, Ngo, Faraji Niri & Marco (2025) — open-access dataset release
  (Data in Brief, 63, 112190)

### 3. Model Formulation
- 3.1 Thermal model (Hatchard–Kim four-reaction, lumped body) [from Phase 2]
- 3.2 Pressure sub-model (this work's contribution) [from Phase 3]
- 3.3 Numerical implementation (ode15s, pure MATLAB, no Simulink dependency)

### 4. Calibration Methodology
- 4.1 Phase 2: thermal parameter fitting (training cell 2, locked)
- 4.2 Phase 3: pressure parameter fitting (training cell 2, cascaded)
- 4.3 Pre-vent windowing strategy

### 5. Results and Discussion
- 5.1 Thermal model performance (Phase 2 recap)
- 5.2 Pressure model performance (Phase 3 results, validation table)
- 5.3 Comparison with Ostanek (2025) base model [NEW SECTION — see Part 2]
- 5.4 Comparison with Chen et al. (2025) [NEW SECTION — see Part 2]
- 5.5 Limitations and future work

### 6. Conclusions

### Data Availability
- GitHub repository link
- Reference to Gulsoy et al. open-access Mendeley dataset

### Appendix
- Full parameter table
- Sensitivity of pressure to nu_e, nu_O2, V_int (optional)


---

## PART 2 — OSTANEK COMPARISON DISCUSSION (Section 5.3)

### 5.3  Comparison with the Ostanek (2025) Base Model

The thermal model in the present work is derived from the open-source Simulink
implementation published by Ostanek (2025) through the Purdue University Research
Repository. That implementation provides a 0D lumped-body formulation of the
Hatchard–Kim four-reaction abuse mechanism (SEI decomposition, anode–electrolyte
reaction, cathode decomposition, and electrolyte decomposition), targeting 18650-format
NMC cells under oven-test conditions.

The present work extends the Ostanek base model in four respects:

**1. Cell format and parameterisation.** The Ostanek model is parameterised for
18650-format cells (~2.5 Ah). The present implementation is re-parameterised for
21700-format NMC cells (~5.0 Ah) using the cell-level properties reported in the
Gulsoy et al. (2024) experimental campaign, including revised values for cell mass,
heat capacity, external surface area, and active-material inventories.

**2. Addition of a gas-pressure state variable.** The Ostanek model solves six coupled
ODEs (four reaction extents, temperature, and SEI thickness). The present model
appends a seventh state variable — the moles of gas accumulated in the sealed cell
headspace — and computes internal pressure via the ideal gas law. This extension
enables direct comparison with the in-situ pressure measurements of Gulsoy et al.,
which are not addressable by the original thermal-only formulation.

**3. Kinetic-driven gas generation.** Rather than prescribing gas generation as an
empirical function of temperature (as in some pressure models in the literature), the
gas generation rate in the present model is coupled to the existing Hatchard–Kim
reaction rates. Specifically, the electrolyte decomposition rate (dc_e/dt, Reaction 4)
and cathode conversion rate (dalpha/dt, Reaction 3) each contribute to the total
gas production through fitted stoichiometric yield coefficients (nu_e and nu_O2).
This approach ensures thermodynamic consistency between the thermal and pressure
predictions: the same reactions that generate heat also generate gas, with no
additional empirical functions.

**4. Pure-MATLAB ODE implementation.** The Ostanek model is distributed as a
Simulink block diagram. The present implementation re-codes the entire system
as a pure-MATLAB function using ode15s, which eliminates the Simulink licence
dependency, enables straightforward batch execution for parameter fitting, and
reduces per-simulation wall time by approximately 20x on the same hardware. This
is particularly relevant for the iterative optimisation required in Phase 3 pressure
fitting, where thousands of forward model evaluations are needed.

The Ostanek model, being thermal-only, cannot predict internal pressure and therefore
cannot be directly benchmarked against the Gulsoy et al. pressure data. The comparison
is therefore structural rather than quantitative: the present work preserves the validated
Hatchard–Kim kinetic framework while adding the pressure observable as a new model
output.


---

## PART 3 — CHEN ET AL. COMPARISON (Section 5.4)

### 5.4  Comparison with Chen et al. (2025)

Chen, Gulsoy, Barai, Nakhanivej, Loveridge and Marco (2025) published a lumped
thermal runaway model with an internal pressure sub-model for 21700 NMC cells,
validated against the same Gulsoy et al. experimental dataset used in the present work.
Their model, published in the Journal of Energy Storage (vol. 116, 116066), represents
the most direct comparator to the present study. The key differences are summarised
below.

| Feature | Chen et al. (2025) | Present work |
|---|---|---|
| Kinetic mechanism | Custom lumped (not Hatchard–Kim) | Hatchard–Kim four-reaction |
| Gas sources | Vapour–liquid equilibrium + decomposition | Decomposition only (vapourisation disabled) |
| Vent mechanism | Included (vent area parametric study) | Not included (pre-vent only) |
| Pressure validation range | Full trace including post-vent | Pre-vent window only |
| Thermal validation | Internal thermocouple | Internal thermocouple |
| Software | Not specified | Pure MATLAB (ode15s), open-source |
| Reproducibility | Closed code | GitHub repository |

**Complementary strengths.** The Chen et al. model includes a vent-release mechanism
and vapour–liquid equilibrium, enabling prediction of post-vent pressure dynamics
(depressurisation, secondary pressure peaks). The present model does not include
venting but provides a more physically transparent gas-generation pathway by coupling
directly to the well-established Hatchard–Kim reaction kinetics. This makes the model
more interpretable for sensitivity analysis (e.g., which reaction dominates pre-vent gas
production?) and more portable to other cell chemistries where Hatchard–Kim parameters
are available.

**Key finding confirmed by both models.** Both the present work and Chen et al. (2025)
find that the pre-vent internal pressure in 21700 NMC cells reaches approximately
20 bar before safety-vent activation. The present model achieves a training-cell pre-vent
pressure RMSE of 1.08 bar, demonstrating that the Hatchard–Kim reaction framework —
originally developed for temperature prediction only — can be extended to capture
the pressure observable with quantitative accuracy when appropriate stoichiometric
yield coefficients are introduced.

**Relevance to the WMG research programme.** The present work validates against and
extends the experimental methodology developed at WMG by Gulsoy et al. (2024, 2025),
and provides an independent modelling perspective that complements the Chen et al. (2025)
study from the same research group. The open-source nature of both the Ostanek (2025)
base model and the present extension facilitates community reproduction and further
development.


---

## PART 4 — CONTRIBUTION STATEMENT (for Section 1)

The specific contributions of this work are:

1. Extension of the Ostanek (2025) open-source Hatchard–Kim thermal runaway model
   with a gas-pressure state variable that couples internal pressure prediction to
   the existing reaction kinetics, without introducing additional empirical temperature
   functions.

2. Re-parameterisation from 18650 to 21700 NMC format and calibration against the
   first-of-its-kind simultaneous internal temperature and pressure dataset of
   Gulsoy et al. (2024, 2025).

3. Demonstration that the Hatchard–Kim reaction framework can predict pre-vent
   internal gas pressure with a training RMSE of 1.08 bar (validation: 2.4–2.5 bar)
   using only three additional fitted parameters (two stoichiometric yields and the
   internal free volume).

4. Open-source release of the complete model implementation in pure MATLAB,
   including parameter fitting scripts and validated parameter sets, via a public
   GitHub repository.


---

## PART 5 — PROF. MARCO EMAIL

**Context:** This email is sent AFTER the preprint is posted (arXiv or SSRN) and
the manuscript is submitted to JoPSA. The email references Faisal's three confirmed
research outputs only (Paper A Energies submission, Paper B Sobol UQ study, Paper C
this work). No withdrawn ICHEV papers.

---

**To:** J.Marco@warwick.ac.uk
**From:** u2796240@uel.ac.uk
**Subject:** PhD enquiry — ETP thermal runaway modelling (validated against Gulsoy et al. dataset)

Dear Professor Marco,

I am writing to enquire about PhD supervision opportunities within your group at
WMG, starting October 2026. My research background is in electric vehicle powertrain
modelling and lithium-ion battery safety, and I believe there is strong alignment with
the thermal runaway characterisation work led by your group.

I have recently completed a coupled Electro-Thermal-Pressure (ETP) model for
prediction of internal gas pressure evolution during thermal runaway in 21700 NMC
cells, validated against the open-access dataset published by Gulsoy et al. (2024, 2025)
from your group. The model extends the Ostanek (2025) Hatchard–Kim Simulink
implementation with a kinetic-driven gas-pressure state variable and achieves a
pre-vent pressure RMSE of 1.08 bar on the training cell. A preprint is available at
[PREPRINT URL] and the manuscript has been submitted to the Journal of Power Sources
Advances.

I am aware of the recent paper by Chen et al. (2025) from your group, which presents
a complementary approach using vapour–liquid equilibrium and a vent mechanism. My
work differs in coupling gas generation directly to the Hatchard–Kim reaction rates,
which I believe offers advantages for parametric sensitivity analysis and cross-chemistry
portability. I would welcome the opportunity to discuss how these approaches might
be combined or extended in a doctoral programme.

My research outputs to date are:

1. "Cross-Cycle Energy-Consumption Characterisation of a Tesla Model 3-Class BEV
   Powertrain Under Field-Oriented Control: NEDC, WLTP and UDDS" — submitted to
   Energies (MDPI), under review. First and corresponding author.

2. Global sensitivity analysis (Sobol indices) of BEV powertrain energy consumption
   under parametric uncertainty — manuscript in preparation.

3. "A Coupled Electro-Thermal-Pressure Model for Prediction of Internal Gas Pressure
   Evolution During Thermal Runaway in 21700 NMC Lithium-Ion Cells" — submitted to
   Journal of Power Sources Advances, preprint at [URL]. First and corresponding author.

I hold an MSc in Electric Vehicle Engineering (First Class Honours) from the
University of East London, supervised by Dr Thayaalan Sutharssan, and a BTech in
Mechanical Engineering. I am currently on a UK Student Visa (Graduate Route) and
available to start in October 2026.

I would be very grateful for the opportunity to discuss potential PhD projects. I am
happy to provide any additional materials.

Kind regards,
Faisal Shah
MSc Electric Vehicle Engineering (Distinction)
University of East London
u2796240@uel.ac.uk
[PHONE NUMBER — update before sending]


---

## PART 6 — PREPRINT & SUBMISSION CHECKLIST

### Before submission to JoPSA:
- [ ] Convert pressure sections (from this document) to LaTeX in Overleaf
- [ ] Merge with existing thermal model sections (Phase 2)
- [ ] Insert publication-quality figures (phase3_prevent_cell1/2/3.pdf)
- [ ] Finalise parameter table (Appendix)
- [ ] Write abstract (combine thermal + pressure contributions)
- [ ] Add GenAI disclosure (Claude, same format as Paper A)
- [ ] Suggested reviewers: Chen/Gulsoy/Marco (WMG), Ostanek (Purdue),
      Coman (if still active), Parhizi (UL/ESRI)
- [ ] Cover letter referencing Chen et al. (2025) as complementary work

### Preprint posting:
- [ ] Post to SSRN (Elsevier preprint server, free, DOI assigned)
      OR arXiv (eess.SY or physics.app-ph)
- [ ] Update GitHub README with preprint link
- [ ] Wait 1-2 days for indexing before emailing Prof. Marco

### Prof. Marco email:
- [ ] Insert preprint URL
- [ ] Update phone number
- [ ] Send from UEL email (u2796240@uel.ac.uk)
- [ ] CC Dr Sutharssan (optional, shows supervisor awareness)

### Suggested reviewers for JoPSA:
1. Prof. James Marco (WMG/Warwick) — dataset author, group leader
   J.Marco@warwick.ac.uk
2. Prof. Jason Ostanek (Purdue) — base model author
   jostanek@purdue.edu
3. Dr Begum Gulsoy (WMG/Warwick) — experimental dataset first author
4. Dr Mohammad Parhizi (UL/ESRI) — numerical methods for TR models
5. Prof. Joris de Hoog (if reachable) — alternative TR modelling

### Note on suggesting Marco/Gulsoy as reviewers:
This is a judgment call. Suggesting them signals confidence and transparency,
but the editor may exclude them due to proximity to the validation data.
Ostanek and Parhizi are safer primary suggestions. List Marco/Gulsoy as
"non-preferred but acknowledged" in the cover letter if JoPSA allows that format.


---

*End of Paper C complete draft package.*
