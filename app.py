import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

def simulate(
    start_age=55,
    use_ssa=True,
    use_il7=True,
    use_rapa=True,
    T0=1.0, kT=0.035,
    kI=0.015,
    N0=1.0, I0=0.0,
    d0=0.05, aI=0.25,
    b_ssa=0.50,
    r_il7=0.25,
    r_rapa=0.40,
    max_age=90
):
    ages = np.arange(0, max_age + 1)
    T = np.zeros_like(ages, dtype=float)
    I = np.zeros_like(ages, dtype=float)
    N = np.zeros_like(ages, dtype=float)

    I[0] = I0
    N[0] = N0

    kI_current = kI

    for idx, t in enumerate(ages):
        T[idx] = T0 * np.exp(-kT * t)

        if use_ssa and t >= start_age:
            T[idx] *= (1.0 + b_ssa)

        d = d0 + aI * I[idx]

        if use_il7 and t >= start_age:
            d *= (1.0 - r_il7)

        if idx < len(ages) - 1:
            N[idx + 1] = max(0.0, N[idx] + T[idx] - d * N[idx])

            if use_rapa and t >= start_age:
                kI_current = kI * (1.0 - r_rapa)
            else:
                kI_current = kI

            I[idx + 1] = max(0.0, I[idx] + kI_current)

    return pd.DataFrame({"age": ages, "T": T, "I": I, "N": N})

def immune_age(baseline_df, treated_df, target_age=80):
    Nt = float(treated_df.loc[treated_df["age"] == target_age, "N"].values[0])
    diffs = (baseline_df["N"] - Nt).abs()
    best_idx = int(diffs.idxmin())
    return int(baseline_df.loc[best_idx, "age"])

st.title("Stacked Immune Rejuvenation Simulator (Toy Model)")

st.markdown(
"""
This model simulates thymic decline, inflammaging, and naive T-cell dynamics.
Interventions can be toggled to observe theoretical immune age shifts.
"""
)

start_age = st.slider("Intervention Start Age", 30, 80, 55)

col1, col2, col3 = st.columns(3)
with col1:
    use_ssa = st.checkbox("SSA (Boost Thymic Output)", value=True)
with col2:
    use_il7 = st.checkbox("IL-7 (Improve Naive T-cell Survival)", value=True)
with col3:
    use_rapa = st.checkbox("Rapamycin (Reduce Inflammaging)", value=True)

baseline = simulate(start_age=start_age, use_ssa=False, use_il7=False, use_rapa=False)
treated = simulate(start_age=start_age, use_ssa=use_ssa, use_il7=use_il7, use_rapa=use_rapa)

target_age = 80
ia = immune_age(baseline, treated, target_age=target_age)

st.metric("Immune Age (Proxy at Age 80)", f"{ia} years")

# Plot naive T-cell pool
fig1 = plt.figure()
plt.plot(baseline["age"], baseline["N"], label="Baseline")
plt.plot(treated["age"], treated["N"], label="With Interventions")
plt.xlabel("Age")
plt.ylabel("Naive T-cell Pool (Normalized)")
plt.legend()
st.pyplot(fig1)

# Plot inflammaging burden
fig2 = plt.figure()
plt.plot(baseline["age"], baseline["I"], label="Baseline")
plt.plot(treated["age"], treated["I"], label="With Interventions")
plt.xlabel("Age")
plt.ylabel("Inflammaging Burden (Arbitrary Units)")
plt.legend()
st.pyplot(fig2)

st.caption("Note: This is a simplified educational model. Parameters are illustrative, not clinically validated.")
