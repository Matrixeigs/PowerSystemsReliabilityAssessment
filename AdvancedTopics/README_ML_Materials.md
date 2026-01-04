# Machine Learning Enhanced Reliability Evaluation - Course Materials

This directory contains comprehensive materials for the **Machine Learning Enhanced Reliability Evaluation** lecture, integrating theory (LaTeX slides) with practice (Julia code).

## 📁 Files Overview

### 📊 Lecture Slides
- **`MachineLearning.tex`** - Main lecture presentation (compile with LaTeX)
  - Comprehensive theory and derivations
  - Integrated with course HL1, HL2, and Sequential MCS lectures
  - Based on Lin & Shang (IEEE TPWRS 2021)

### 💻 Julia Code Examples

1. **`ML_Tutorial_StepByStep.jl`** ⭐ **START HERE**
   - Beginner-friendly walkthrough
   - Each concept explained individually with simple examples
   - Generates 3 diagnostic plots
   - ~30 minutes to work through
   - **Run first to understand basics!**

2. **`ML_Reliability_Examples.jl`** ⭐ **MAIN IMPLEMENTATION**
   - Complete, production-ready implementation
   - Fully integrated DBBN + PHM + RH-UC + Sequential MCS
   - Demonstrates full case study
   - Can take 10-30 minutes to run (depending on scenarios)
   - **Run after understanding step-by-step tutorial**

3. **`ML_Reliability_Guide.md`** 📖 **COMPREHENSIVE GUIDE**
   - Detailed documentation
   - Maps code to lecture slides
   - Exercises and experiments
   - Troubleshooting tips
   - **Keep open as reference while coding**

4. **`README_ML_Materials.md`** (this file)
   - Navigation guide
   - Quick start instructions

## 🚀 Quick Start Guide

### For Complete Beginners

```bash
# Step 1: Understand the basics
julia ML_Tutorial_StepByStep.jl

# Step 2: Review the generated plots
# - step1_dbbn_samples.png
# - step2_phm_loading.png
# - step5_convergence.png

# Step 3: Read corresponding lecture slides
# Open MachineLearning.tex (or compiled PDF)
# Review Sections 1-3, 6

# Step 4: Run full implementation
julia ML_Reliability_Examples.jl

# Step 5: Experiment with exercises
# See ML_Reliability_Guide.md for ideas
```

### For Advanced Students

```bash
# Go directly to main implementation
julia ML_Reliability_Examples.jl

# Then experiment with modifications
# Refer to ML_Reliability_Guide.md for exercises
```

### For Instructors

```julia
# Use step-by-step tutorial in class demonstration
include("ML_Tutorial_StepByStep.jl")

# Assign exercises from guide for homework
# See "Experimental Exercises" section in ML_Reliability_Guide.md

# Use full implementation for project work
# Students can extend with storage, multi-area, etc.
```

## 📚 Learning Path

### Week 1: Foundations
- [ ] Read MachineLearning.tex slides (Sections 1-2: Introduction, DBBN, PHM)
- [ ] Run `ML_Tutorial_StepByStep.jl` - Steps 1-2
- [ ] Exercise: Modify DBBN distribution for solar (Guide, Exercise 1)

### Week 2: Integration
- [ ] Read MachineLearning.tex slides (Sections 3-4: Sampling, RH-UC)
- [ ] Run `ML_Tutorial_StepByStep.jl` - Steps 3-4
- [ ] Exercise: PHM sensitivity analysis (Guide, Exercise 2)

### Week 3: Full Framework
- [ ] Read MachineLearning.tex slides (Sections 5-7: SMC, Algorithm, Summary)
- [ ] Run `ML_Reliability_Examples.jl` - Full case study
- [ ] Exercise: Rolling-horizon length study (Guide, Exercise 3)

### Week 4: Advanced Topics
- [ ] Implement storage double-layer model (Slides Section 5)
- [ ] Compare with traditional MCS (Guide, Exercise 4)
- [ ] Final project: Extend framework with custom features

## 🔗 Connections to Course

### Prerequisite Lectures
You should understand these concepts first:
- **HL1 (Generation Adequacy):** COPT, LOLE, LOEE
- **HL2 (Composite System):** DC-OPF, power flow, network constraints
- **Sequential MCS:** Time-series simulation, TTF/TTR, convergence
- **Markov Models:** State space, transition rates, exponential distribution

### This Lecture Extends
- HL2 with **renewable uncertainty** (DBBN)
- Sequential MCS with **stress-dependent failures** (PHM)
- Standard UC with **rolling-horizon** (RH-UC)
- All three with **tight coupling** (integrated framework)

### Next Topics
This framework prepares you for:
- **Resilience Assessment:** Extreme events, cascading failures
- **Multi-area Reliability:** Inter-regional dependencies
- **Market-Reliability Co-optimization:** Coupling with electricity markets
- **Deep Learning Methods:** Neural networks for forecasting

## 🛠️ Setup Instructions

### Install Julia
```bash
# Download from https://julialang.org/downloads/
# Recommended: Julia 1.9 or later
```

### Install Required Packages
```julia
using Pkg
Pkg.add([
    "Random",
    "Distributions",
    "Statistics",
    "JuMP",
    "HiGHS",
    "Plots"
])
```

### Verify Installation
```julia
# Run this to test:
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@objective(model, Max, x)
@constraint(model, x <= 5)
optimize!(model)
println("Test passed: x = ", value(x))  # Should print 5.0
```

## 📊 Expected Outputs

### Tutorial Outputs
- **Console:** Step-by-step explanations with numerical examples
- **Plots (3 files):**
  - `step1_dbbn_samples.png` - Histogram of renewable samples
  - `step2_phm_loading.png` - Hazard rate vs loading
  - `step5_convergence.png` - CV convergence curve

### Main Example Outputs
- **Console:** Simulation progress, convergence status, final results
- **Plot (1 file):**
  - `reliability_convergence.png` - Three-panel convergence plot

### Typical Results
For the 3-generator test system:
- LOLP ≈ 0.005 - 0.015 (0.5% - 1.5%)
- EENS ≈ 50 - 200 MWh/year
- Convergence: 100-300 scenarios for CV < 5%

## 🐛 Common Issues & Solutions

### "Package not found"
```julia
# Install missing package
Pkg.add("PackageName")
```

### "UC optimization infeasible"
- Increase total generation capacity
- Check `generators` definition in code
- Reduce peak load in `load_profile`

### "Slow execution"
- Reduce simulation horizon: `T = 720` (1 month instead of 1 year)
- Reduce scenarios: `max_scenarios = 50`
- Reduce horizon length: `horizon_length = 3`

### "CV not converging"
- Increase `max_scenarios` to 500-1000
- Relax `cv_target` to 0.10
- Check that EENS > 0 (system must have some risk)

### "Plot not displaying"
```julia
# Save plot instead of displaying
using Plots
plot(...)
savefig("myplot.png")
```

## 📝 Assessment Rubric

Students should demonstrate understanding of:

### Knowledge (40%)
- [ ] Can explain DBBN advantages over parametric distributions
- [ ] Can derive PHM hazard function
- [ ] Can describe rolling-horizon mechanism
- [ ] Can calculate reliability indices (LOLP, EENS)

### Implementation (40%)
- [ ] Successfully runs tutorial and main example
- [ ] Completes at least 2 exercises from guide
- [ ] Code produces reasonable results
- [ ] Can debug and troubleshoot issues

### Analysis (20%)
- [ ] Interprets convergence plots correctly
- [ ] Explains coupling loop with evidence from results
- [ ] Compares traditional vs ML-enhanced approaches
- [ ] Proposes meaningful extensions

## 🎯 Key Takeaways

After completing these materials, students should:

1. **Understand integration:** ML (DBBN) + Optimization (UC) + Simulation (MCS) work together
2. **Recognize coupling:** Operational decisions affect reliability dynamically
3. **Appreciate realism:** Stress-dependent failures are more realistic than constant rates
4. **Apply methods:** Can implement and extend the framework for research
5. **Connect theory-practice:** Code implements mathematical equations from slides

## 📚 References

### Primary Source
- Lin, Y., & Shang, C. (2021). "Reliability Evaluation on a Joint Machine Learning and Optimization Framework." *IEEE Transactions on Power Systems*, 36(4), 3421-3434.

### Background Reading
- **DBBN:** Koller & Friedman, "Probabilistic Graphical Models"
- **PHM:** Cox, "Regression Models and Life-Tables"
- **Sequential MCS:** Billinton & Li, "Reliability Assessment Using Monte Carlo Methods"
- **Unit Commitment:** Wood, Wollenberg & Sheblé, "Power Generation, Operation, and Control"

### Course Materials
- Lecture slides: `MachineLearning.tex`
- Previous lectures: HL1, HL2, Sequential MCS, Markov Models
- Course textbook: Billinton & Allan, "Reliability Evaluation of Power Systems"

## 🤝 Contributing

Students are encouraged to:
- Report bugs or unclear explanations
- Suggest additional exercises
- Share interesting extensions
- Submit improved visualizations

Contact course instructor or create issues in course repository.

## 📜 License

Educational materials for "Electric Power System Reliability Assessment" course.
Code examples may be used for academic purposes with attribution.

---

**Happy Learning!** 🎓

*Last updated: January 2026*
*Course: Electric Power System Reliability Assessment*
*Topic: ML-Enhanced Reliability Evaluation*
