# Machine Learning Enhanced Reliability Evaluation - Code Guide

This guide accompanies the **MachineLearning.tex** lecture slides and provides hands-on Julia implementations of the key methodologies.

## 📚 Overview

The code demonstrates the integrated framework from Lin & Shang (IEEE TPWRS 2021), combining:
1. **DBBN** (Dynamic Bayesian Belief Network) for renewable forecasting
2. **PHM** (Proportional Hazard Model) for aging-aware failures
3. **RH-UC** (Rolling-Horizon Unit Commitment) for operational decisions
4. **Sequential MCS** (Monte Carlo Simulation) for reliability assessment

## 🎯 Learning Objectives

After working through this code, you should be able to:
- [ ] Implement inverse transform sampling for DBBN-based renewable generation
- [ ] Calculate hazard rates and failure probabilities using PHM
- [ ] Understand the coupling between UC dispatch and component aging
- [ ] Set up and solve a rolling-horizon optimization problem
- [ ] Run Sequential Monte Carlo with convergence monitoring
- [ ] Interpret reliability indices (LOLP, EENS) in context

## 🚀 Getting Started

### Prerequisites

Install required Julia packages:
```julia
using Pkg
Pkg.add(["Random", "Distributions", "Statistics", "JuMP", "HiGHS", "Plots"])
```

### Running the Code

**Option 1: Full demonstration**
```julia
include("ML_Reliability_Examples.jl")
```

**Option 2: Interactive exploration**
```julia
include("ML_Reliability_Examples.jl")

# Example: Sample from DBBN
dist = create_example_wind_distribution(14)  # 2 PM
samples = [sample_renewable(dist) for _ in 1:1000]
using Plots
histogram(samples, bins=20, xlabel="Wind Power (MW)", ylabel="Frequency", title="DBBN Samples at Hour 14")
```

## 📖 Code Structure and Lecture Connections

### Part 1: DBBN for Renewables (Lecture Section 1)

**Corresponds to slides:**
- Frame 1.1: "The Renewable Prediction Challenge"
- Frame 1.2: "Why DBBN?"
- Frame 3.2: "Sampling Renewable Generation from DBBN"

**Key functions:**
```julia
# Create time-dependent wind distribution (simulates trained DBBN output)
dist = create_example_wind_distribution(hour)

# Sample using inverse transform (Algorithm from slide 3.2)
wind_power = sample_renewable(dist)
```

**Mathematical connection:**
The code implements the inverse transform sampling from slides:
```
1. Build CDF: P(t,v) = Σᵢ₌₀ᵛ p(t,i)
2. Sample u ~ Uniform(0,1)
3. Find vᵢ: u ∈ [P(t,vᵢ), P(t,vᵢ₊₁)]
```

**Exercise:**
Modify `create_example_wind_distribution()` to include solar generation with different diurnal patterns (high at noon, zero at night).

---

### Part 2: Proportional Hazard Model (Lecture Section 2)

**Corresponds to slides:**
- Frame 2.1: "Beyond Constant Failure Rates"
- Frame 2.3: "PHM Hazard Function Formulation"
- Frame 2.4: "From Hazard to Outage Probability"

**Key functions:**
```julia
# Define PHM parameters (from slide 2.3)
phm = PHMParameters(
    a = 8760.0,   # Scale (hours)
    b = 1.2,      # Shape (aging effect)
    γ₁ = 0.5,     # Temperature sensitivity
    γ₂ = 2.0      # Loading sensitivity
)

# Compute hazard rate (slide 2.3 equation)
h = compute_hazard_rate(phm, t=1000.0, temp_category=1, loading_rate=0.8)

# Update failure probability recursively (slide 2.4)
F_new = update_outage_probability(F_prev, phm, t, temp, loading)
```

**Mathematical connection:**
Implements the hazard function from slide 2.3:
```
h(t, Z(t)) = a^(-b) * b * t^(b-1) * exp(γ₁*Z₁(t) + γ₂*Z₂(t)²)
```

And the recursive formula from slide 2.4:
```
F(t) = 1 - exp(-∫h(x,Z)dx) + F(t-1)*exp(-∫h(x,Z)dx)
```

**Exercise:**
Run `demonstrate_phm_effect()` and observe how loading and temperature affect failure rates. Create a heatmap showing hazard rate as a function of loading (x-axis) and temperature (y-axis).

---

### Part 3: Rolling-Horizon Unit Commitment (Lecture Section 4)

**Corresponds to slides:**
- Frame 4.1: "The Operational Decision Challenge"
- Frame 4.2: "Rolling-Horizon Concept"
- Frame 4.3: "RH-UC Control Loop Algorithm"
- Frame 4.4: "Outage-Aware UC Constraint Modifications"

**Key functions:**
```julia
generation, load_shed, obj = solve_simplified_uc(
    generators,
    load_window,      # Load over [t, t+m]
    renewable_window, # Renewable forecast
    outage_states,    # Which generators failed
    horizon_length
)
```

**Mathematical connection:**
Solves the optimization problem from slides:
```
Minimize: Σ(cost * generation) + penalty * load_shedding

Subject to:
- Power balance: Σ generation + renewables = load - shedding
- Generator limits (adjusted for outages, slide 4.4)
- Ramp rates (relaxed if failed, slide 4.4)
```

**Note:** The code uses a LINEAR relaxation for educational clarity. A full MILP implementation would include binary commitment variables u_{g,t}.

**Exercise:**
Modify `solve_simplified_uc()` to add binary on/off variables and minimum up/down time constraints. Compare results with the linear relaxation.

---

### Part 4: Sequential Monte Carlo Framework (Lecture Section 6)

**Corresponds to slides:**
- Frame 6.1: "Recall: Sequential MCS from Earlier Lectures"
- Frame 6.2: "Anatomy of One SMC Scenario"
- Frame 6.3: "Reliability Indices"
- Frame 6.4: "Convergence Criterion: Coefficient of Variation"

**Key functions:**
```julia
# Simulate one complete year (implements algorithm from slide 7.1)
indices = simulate_one_scenario(generators, load_profile, T, phm_params)

# Run multiple scenarios until convergence (slide 6.4)
results = run_sequential_monte_carlo(
    generators, load_profile, T, phm_params,
    cv_target = 0.05  # Stop when CV ≤ 5%
)
```

**Mathematical connection:**
Implements the convergence criterion from slide 6.4:
```
CV = σ / (μ√n) where σ = std(LOLP samples), μ = mean(LOLP)
```

Computes reliability indices from slide 6.3:
```
LOLP = (# hours with load shedding) / T
EENS = Σ load_shedding_t (MWh/year)
```

**Exercise:**
Vary the `cv_target` from 0.01 to 0.10 and observe the trade-off between accuracy and computational time. Plot CV vs. number of scenarios.

---

### Part 5: Integrated Example (Complete Framework)

**Corresponds to slide:**
- Frame 7.1: "Complete Algorithm (Integrated Framework)"
- Frame 7.2: "Key Methodological Couplings"

**Key function:**
```julia
results = run_example_case_study()
```

This demonstrates the **four critical coupling loops** from slide 7.2:
1. **DBBN → UC:** Renewable samples → UC boundary conditions
2. **UC → PHM:** Dispatch → loading → updated failure probability
3. **PHM → UC:** Failures → constrained UC problem
4. **Rolling Horizon → Reliability:** Adaptive decisions reduce forecast error

## 🔬 Experimental Exercises

### Exercise 1: Impact of Renewable Penetration
Modify the case study to include more wind capacity (add a 4th generator with zero cost but intermittent output). Observe how LOLP and EENS change.

```julia
# Add to run_example_case_study()
generators = [
    Generator(1, 50.0, 200.0, 30.0, 50.0, 50.0, 4, 4),
    Generator(2, 30.0, 150.0, 40.0, 40.0, 40.0, 2, 2),
    Generator(3, 20.0, 100.0, 50.0, 30.0, 30.0, 1, 1),
    Generator(4, 0.0, 150.0, 0.0, 150.0, 150.0, 0, 0),  # Wind (variable)
]
```

### Exercise 2: PHM Parameter Sensitivity
Create a sensitivity analysis showing how γ₂ (loading sensitivity) affects:
- Average failure frequency
- Correlation between high-load hours and failures

```julia
function sensitivity_analysis()
    γ₂_values = [0.5, 1.0, 1.5, 2.0, 2.5]
    results = []

    for γ₂ in γ₂_values
        phm_params = [PHMParameters(8760.0, 1.2, 0.5, γ₂) for _ in 1:3]
        result = run_sequential_monte_carlo(..., max_scenarios=100)
        push!(results, result)
    end

    # Plot EENS vs γ₂
end
```

### Exercise 3: Rolling-Horizon Length Study
Compare system performance with different horizon lengths m = {1, 3, 6, 12, 24}:
- m=1: Myopic (high EENS due to commitment inflexibility)
- m=24: Full day-ahead (may be overconfident in forecasts)
- m=6: Sweet spot (paper's choice)

### Exercise 4: Comparison with Traditional MCS
Implement a traditional Sequential MCS with:
- Fixed failure rates (no PHM)
- Deterministic renewable forecast (no DBBN)
- Day-ahead UC (no rolling horizon)

Compare LOLP and EENS with the ML-enhanced framework.

## 📊 Understanding the Output

### Convergence Plot
The code generates `reliability_convergence.png` showing:
- **Top panel:** LOLP convergence - should stabilize as scenarios increase
- **Middle panel:** EENS convergence - larger variance than LOLP
- **Bottom panel:** CV metric - should decrease below 5% target

### Key Observations to Note:
1. **DBBN effect:** Wind samples vary by hour (check histogram at different times)
2. **PHM effect:** Failures cluster during high-load periods (stress-dependent)
3. **RH-UC effect:** Load shedding is minimized by looking ahead m hours
4. **Convergence:** Typically needs 100-500 scenarios for CV < 5%

## 🔗 Connection to Other Course Materials

### Builds on Previous Lectures:
- **HL1 (Generating Adequacy):** Still assessing capacity adequacy, but with time-varying renewables
- **HL2 (Composite System):** Could extend to include transmission constraints in UC
- **Sequential MCS:** Extends basic framework with ML and dynamic failure rates
- **Markov Models:** PHM generalizes constant-rate exponential failures

### Key Differences from Traditional Methods:
| Aspect | Traditional | This Framework |
|--------|-------------|----------------|
| Renewables | Deterministic/Gaussian | DBBN (non-parametric) |
| Failures | λ = constant | PHM (stress-dependent) |
| Operations | Day-ahead UC | Rolling-horizon |
| Coupling | One-way | Bidirectional feedback |

## 🐛 Troubleshooting

**Issue:** "HiGHS solver not found"
```julia
Pkg.add("HiGHS")
```

**Issue:** UC optimization infeasible
- Check that total generation capacity > peak load
- Verify outage states don't remove too much capacity
- Increase generator ramp rates if needed

**Issue:** CV not converging
- Increase `max_scenarios` (try 500-1000)
- Relax `cv_target` to 0.10 for initial testing
- Check that EENS > 0 (system must have some risk)

**Issue:** Slow performance
- Reduce T (simulate 1 month instead of 1 year)
- Reduce horizon_length to 3 instead of 6
- Use fewer generators (3 is minimum for interesting behavior)

## 📚 Further Reading

1. **Original paper:** Lin & Shang, "Reliability Evaluation on a Joint Machine Learning and Optimization Framework," IEEE TPWRS, 2021
2. **DBBN:** Murphy, "Dynamic Bayesian Networks"
3. **PHM:** Cox, "Regression Models and Life-Tables"
4. **UC:** Padhy, "Unit Commitment - A Bibliographical Survey"
5. **MCS in Power Systems:** Billinton & Li, "Reliability Assessment of Electric Power Systems Using Monte Carlo Methods"

## 💡 Tips for Learning

1. **Start simple:** Run `demonstrate_phm_effect()` first to understand one component
2. **Visualize:** Create plots for each component before integrating
3. **Compare:** Run with/without each feature to see its impact
4. **Experiment:** Change parameters and observe sensitivity
5. **Debug:** Add print statements to `simulate_one_scenario()` to see hourly evolution
6. **Extend:** Implement the optional storage model from Section 5

## ✅ Validation Checklist

Your implementation is working correctly if:
- [ ] DBBN samples match the expected distribution (histogram test)
- [ ] PHM hazard rate increases with loading and temperature
- [ ] UC respects generation limits and outage constraints
- [ ] Load shedding only occurs when capacity is insufficient
- [ ] CV decreases as more scenarios are simulated
- [ ] LOLP is between 0.001 and 0.1 for a reasonable test system
- [ ] EENS increases with more frequent/longer failures

## 🎓 Assessment Questions

Test your understanding:

1. Why does DBBN output discrete distributions rather than parametric (Gaussian)?
2. How does the Z₂² term in PHM emphasize overloading effects?
3. What happens if the rolling horizon m=1 (pure real-time)?
4. Why is CV based on LOLP rather than EENS for convergence?
5. Identify the feedback loop: How does UC dispatch affect future failures?

Answers are in the lecture slides (Sections 1.1, 2.3, 4.1, 6.4, and 7.2 respectively).

---

**Happy learning!** 🚀

For questions or issues, refer to:
- Lecture slides: `MachineLearning.tex`
- Course materials: Other reliability assessment lectures
- Julia documentation: https://docs.julialang.org
