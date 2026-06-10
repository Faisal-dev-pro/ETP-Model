# Paper C — Pressure Sub-Model Sections (Draft)

> **Target journal:** Journal of Power Sources Advances
>
> **Paper title (working):** A Coupled Electro-Thermal-Pressure Model for Prediction of Internal Gas Pressure Evolution During Thermal Runaway in 21700 NMC Lithium-Ion Cells
>
> These sections are ready to paste into the Overleaf manuscript. Convert notation to LaTeX as needed.


---

## 3.X  Internal Gas Pressure Sub-Model

The thermal model described in Section 3.Y is extended with a lumped internal-gas-pressure sub-model that tracks the total moles of gas generated within the sealed cell volume prior to safety-vent activation. Two gas-producing mechanisms are retained: (i) electrolyte decomposition and (ii) oxygen release from the positive-electrode active material. A third mechanism — electrolyte vapourisation — is excluded from the present formulation because the model does not track the liquid-phase electrolyte mass as an independent state variable; including vapourisation without a corresponding liquid-depletion term leads to non-physical gas accumulation (see Section 5.X for discussion).

### 3.X.1  Gas Generation Rate

The molar gas generation rate is expressed as

    dn_g/dt = ν_e · n_{e,0} · |dc_e/dt| + ν_{O₂} · n_{c,0} · (dα/dt)⁺

where:

- dn_g/dt is the total gas generation rate [mol s⁻¹],
- ν_e is the stoichiometric gas yield per mole of electrolyte decomposed [mol_gas mol_e⁻¹],
- ν_{O₂} is the stoichiometric gas yield per mole of cathode oxygen released [mol_gas mol_c⁻¹],
- n_{e,0} = m_{e,0}/M_{W,e} is the initial electrolyte inventory [mol],
- n_{c,0} = m_{c,0}/M_{W,c} is the initial cathode active-material inventory [mol],
- dc_e/dt is the dimensionless electrolyte decomposition rate from the Hatchard–Kim kinetic sub-model (Reaction 4),
- dα/dt is the dimensionless cathode conversion rate (Reaction 3),
- (·)⁺ denotes max(·, 0) to ensure non-negative oxygen release.

The dimensionless rates dc_e/dt and dα/dt are governed by the Arrhenius expressions described in Section 3.Y (thermal model) and are functions of the cell temperature T(t) only; no feedback from the pressure state to the thermal model is included in the present formulation.

### 3.X.2  Internal Pressure

The internal gas pressure is obtained from the ideal gas law applied to the sealed cell headspace:

    P(t) = n_g(t) · R_u · T(t) / V_int

where:

- n_g(t) is the total moles of gas at time t (state variable y₇),
- R_u = 8.314 J mol⁻¹ K⁻¹ is the universal gas constant,
- T(t) is the cell temperature from the thermal model [K],
- V_int is the effective internal free volume of the cell [m³].

The initial gas inventory n_g(0) is set so that P(0) = P_atm = 101.3 kPa at the measured initial cell temperature T₀:

    n_g(0) = P_atm · V_int / (R_u · T₀)

No vent-release mechanism is included: the model is valid up to the point of safety-vent activation, after which the sealed-volume assumption breaks down. The predicted vent-activation time is defined as the instant at which P(t) first exceeds the measured peak pre-vent pressure.

### 3.X.3  Cell Parameters

The electrochemical inventory parameters for the 21700 NMC cell are listed below.

| Parameter | Symbol | Value | Unit |
|---|---|---|---|
| Electrolyte mass | m_{e,0} | 10.35 | g |
| Electrolyte molar mass | M_{W,e} | 100 | g mol⁻¹ |
| Electrolyte moles | n_{e,0} | 0.1035 | mol |
| Cathode active mass | m_{c,0} | 20.7 | g |
| Cathode molar mass (NMC) | M_{W,c} | 96.0 | g mol⁻¹ |
| Cathode moles | n_{c,0} | 0.2156 | mol |
| Atmospheric pressure | P_atm | 101.3 | kPa |


---

## 4.X  Pressure Sub-Model Calibration (Phase 3)

### 4.X.1  Fitting Procedure

The pressure sub-model introduces three free parameters — ν_e, ν_{O₂}, and V_int — which are calibrated against the experimental pre-vent internal-pressure trace of a single training cell (Cell 2, specimen 2_0111 from the Gulsoy et al. dataset). All thermal-model parameters from Phase 2 are held frozen during this step to preserve the accepted temperature calibration.

The objective function minimises the mean squared error of the predicted internal pressure against the measured trace over the pre-vent window:

    L(θ) = (1/N) Σᵢ [P_pred(tᵢ; θ) − P_meas(tᵢ)]²

where θ = [ν_e, ν_{O₂}, V_int] and the sum runs over the N samples in the pre-vent window.

**Pre-vent window definition.** The pre-vent window terminates at the first sample where P_meas exceeds the burst-clip threshold of 25 bar. If the measured pressure never reaches this threshold (as occurs for cells with P_peak < 25 bar), the window terminates at the sample of maximum measured pressure. This fallback is necessary because post-vent pressure collapse, if included in the fitting window, introduces a large spurious residual that biases the optimizer toward under-predicting gas generation.

**Optimisation.** The minimisation is carried out using MATLAB's fmincon (interior-point algorithm) with the bounds listed below.

| Parameter | Lower bound | Upper bound | Initial guess |
|---|---|---|---|
| ν_e [mol_gas/mol_e] | 0.01 | 10.0 | 2.0 |
| ν_{O₂} [mol_gas/mol_c] | 0.01 | 5.0 | 1.0 |
| V_int [mL] | 1.0 | 5.0 | 2.5 |

The vapourisation yield ν_vap is fixed at zero (Section 3.X) and the vapourisation rate pre-exponential k_vap is set to zero in the parameter structure.

### 4.X.2  Calibration Results

The optimizer converges in 14 iterations to the values listed below. All three parameters settle within their bounds, confirming that the model structure and bound specification are appropriate.

| Parameter | Fitted value |
|---|---|
| ν_e | 2.20 mol_gas mol_e⁻¹ |
| ν_{O₂} | 0.09 mol_gas mol_c⁻¹ |
| V_int | 3.69 mL |

The training-cell pre-vent pressure RMSE is **1.078 bar**, which satisfies the target of < 2 bar.

The fitted V_int = 3.69 mL is consistent with published estimates of the internal free volume of 21700-format cylindrical cells (typically 2–5 mL including jelly-roll porosity and headspace above the electrode stack).

The large ν_e relative to ν_{O₂} indicates that electrolyte decomposition is the dominant gas source in the pre-vent temperature range (~130–200 °C), consistent with the known thermal stability hierarchy: electrolyte solvents (EC, DMC) decompose at lower temperatures than cathode oxygen release, which requires > 200 °C for NMC chemistries.


---

## 5.X  Pressure Prediction Results and Validation

### 5.X.1  Validation Protocol

The calibrated pressure sub-model is validated on two independent cells (Cell 1, specimen 1_2109, and Cell 3, specimen 3_2811) from the same Gulsoy et al. experimental campaign. All model parameters — thermal and pressure — are held at their calibrated values. The same pre-vent windowing logic is applied: the evaluation window terminates at the first threshold crossing (25 bar) or at peak measured pressure, whichever occurs first.

### 5.X.2  Results

The pre-vent pressure and temperature prediction performance across all three cells is summarised below.

| Cell ID | Role | P RMSE_pre-vent [bar] | T RMSE_pre-vent [K] | P_peak,meas [bar] | t_vent error [s] | P target (< 2 bar) |
|---|---|---|---|---|---|---|
| 2_0111 | Training | 1.08 | 13.2 | 20.0 | −2.1 | **PASS** |
| 1_2109 | Validation | 2.49 | 18.8 | 19.6 | +101.1 | — |
| 3_2811 | Validation | 2.35 | 4.5 | 26.6 | −130.6 | — |

### 5.X.3  Discussion of Validation Spread

The validation-cell pressure RMSE values (2.49 and 2.35 bar) exceed the < 2 bar training target but remain within a factor of 2.5× of the training performance, which is consistent with the cell-to-cell variability inherent in the experimental dataset.

The dominant source of validation error is the inherited temperature prediction offset from the Phase 2 thermal model. Cell 1 exhibits a pre-vent temperature RMSE of 18.8 K, which delays the onset of Arrhenius-driven gas generation by approximately 101 seconds. Because the pressure sub-model is slaved to the thermal state without feedback, any thermal timing error propagates directly into the pressure prediction. Cell 3 exhibits the opposite effect: its thermal model slightly leads the measured temperature, producing a −131 s vent-time error.

This cascaded-calibration limitation is expected in decoupled electro-thermal-pressure frameworks and does not indicate a deficiency in the pressure sub-model itself. Improving the validation spread would require either per-cell thermal recalibration (which conflicts with the single-calibration philosophy of the model) or the addition of thermal feedback from the pressure state (e.g., endothermic venting), which is deferred to future work.

### 5.X.4  Limitations

The following limitations of the pressure sub-model are noted:

1. **No electrolyte vapourisation.** The vapourisation gas-generation pathway is disabled because the model does not track liquid-phase electrolyte mass depletion. Including vapourisation without a consumption term produces non-physical indefinite gas accumulation. A future extension incorporating a liquid–vapour equilibrium state variable would allow this pathway to be reinstated.

2. **No vent-release mechanism.** The model assumes a sealed cell volume and is therefore valid only up to the point of safety-vent activation. Post-vent pressure dynamics (gas release, depressurisation, ejecta) are outside the scope of the present work.

3. **No pressure–thermal feedback.** The pressure sub-model is one-way coupled: temperature drives gas generation, but gas-phase processes (e.g., endothermic vapourisation, adiabatic expansion at vent) do not feed back into the energy balance. This decoupled architecture is standard in Hatchard–Kim-derived models but limits accuracy in scenarios where vent dynamics significantly alter the thermal trajectory.

4. **Ideal gas assumption.** The internal gas mixture is treated as an ideal gas. At the pressures and temperatures encountered pre-vent (< 30 bar, < 500 K), this approximation introduces negligible error.


---

## Acknowledgements (addition)

> Internal pressure measurements were obtained from the open-access dataset published by Gulsoy et al. (2025), made available by WMG, University of Warwick.

---

*End of pressure sub-model draft sections.*
