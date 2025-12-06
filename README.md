# Power Systems Reliability Assessment

A comprehensive MATLAB-based framework for assessing the reliability of bulk power systems using Monte Carlo Simulation (MCS). This package supports both **Hierarchical Level I (HL1)** and **Hierarchical Level II (HL2)** reliability studies, utilizing the IEEE RTS-24 bus test system.

## Overview

Reliability assessment is critical for planning and operating power systems. This repository implements two primary Monte Carlo simulation techniques:

1.  **Non-Sequential MCS (State Sampling)**: Evaluates independent system states to determine probability-based indices.
2.  **Sequential MCS (Chronological)**: Simulates the system evolution over time to capture frequency and duration indices.

## Hierarchical Levels

### HL1: Generation Adequacy
Focuses solely on the ability of generation resources to meet the total system load. Transmission constraints are ignored (Copper Sheet model).
*   **Key Question**: Is there enough generation capacity?
*   **Indices**: LOLE (Loss of Load Expectation), EENS (Expected Energy Not Supplied).

### HL2: Composite System Adequacy
Evaluates the combined generation and transmission system. It accounts for network topology, line limits, and voltage constraints using DC Optimal Power Flow (DC-OPF).
*   **Key Question**: Can the energy be delivered to the load points?
*   **Indices**: LOLE, EENS, plus Nodal Indices (Bus-specific reliability).

## Methodologies Implemented

### 1. Non-Sequential Monte Carlo (State Sampling)
Located in `Montecarlo_nsq_single/`.
*   **Approach**: Randomly samples component states (Up/Down) based on unavailability probabilities ($U$). Each sampled state is independent.
*   **Process**:
    1.  Sample state vector.
    2.  Run DC-OPF (minimizing load shedding).
    3.  Accumulate failure statistics.
*   **Best for**: Fast calculation of expected values (EENS, LOLE) when chronological correlation is not critical.

### 2. Sequential Monte Carlo (Chronological)
Located in `Montecarlo_seq/`.
*   **Approach**: Simulates the system hour-by-hour over many years. Component states transition based on Time-To-Failure (TTF) and Time-To-Repair (TTR) distributions.
*   **Features**:
    *   **Chronological Load**: Uses the IEEE RTS-79 hourly load profile (8736 hours/year).
    *   **Frequency & Duration**: Calculates how often failures occur and how long they last.
*   **Best for**: Detailed analysis requiring frequency (occ/yr) and duration (hr/occ) indices.

## Prerequisites

*   **MATLAB**: (Tested on R2020b+)
*   **MATPOWER**: Required for power flow and OPF calculations. [Download here](https://matpower.org/).
*   **Optimization Solver**:
    *   **MIPS**: Included with MATPOWER (Default).
    *   **CPLEX / Gurobi**: Highly recommended for faster performance, especially for the Sequential simulation which solves thousands of OPFs.

## Usage

### Running Non-Sequential Simulation
1.  Navigate to the folder:
    ```matlab
    cd Montecarlo_nsq_single
    ```
2.  Run the main script:
    ```matlab
    nsqMain
    ```
3.  **Outputs**: Convergence plots, `reliability_results.mat`, and `nodal_results.csv`.

### Running Sequential Simulation
1.  Navigate to the folder:
    ```matlab
    cd Montecarlo_seq
    ```
2.  Run the main script:
    ```matlab
    seqMain
    ```
3.  **Outputs**: Time-series plots, frequency/duration indices, and `seq_nodal_results.csv`.

## Key Files

*   `nsqMain.m` / `seqMain.m`: Master scripts for simulation control.
*   `mc_sampling.m` / `seq_mcsampling.m`: State generation logic.
*   `mc_simulation.m` / `seq_mcsimulation.m`: System evaluation (OPF) logic.
*   `case24_ieee_rts.m`: MATPOWER case file for the test system.
*   `anloducurve.m`: Generates the chronological load profile (Sequential only).

## License
This project is intended for educational and research purposes.
