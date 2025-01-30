# Neurovascular Unit (NVU) Model Simulation

**By:Elshin Mathias** 

## Description

This MATLAB code implements a comprehensive computational model of a Neurovascular Unit (NVU), integrating the dynamics of a neuron, an astrocyte, a smooth muscle cell (SMC) within the cerebral microvasculature, and the vessel wall mechanics. This model simulates the complex interactions between these cellular components, focusing on the transport of key ions (K+, Na+, Cl-, HCO3-, Ca2+), signaling molecules (IP3, EET, NO, cGMP), and oxygen. The model is designed to explore neurovascular coupling mechanisms, neural signaling, homeostasis, and the interplay between neuronal activity, astrocytic function, and vascular responses.

**Key Features of the Model Include:**

*   **Integrated Multi-Cellular Simulation:** Simulates the dynamics of the neuron, astrocyte, smooth muscle cell, and vessel wall as interconnected components.
*   **Ion Transport:** Models various ion channels and transporters, including Na+/K+ pump, K+ channels (BK, KDR, and KA), TRPV4 channels, NMDA channels, NKCC1, and KCC1.
*   **Volume and Pressure Dynamics:** Simulates changes in volume and surface area ratios in the astrocyte and synaptic cleft, as well as vessel radius.
*   **Calcium and Signaling Pathways:** Incorporates IP3 and calcium-dependent pathways, TRPV4 channel activation, NO signaling, and cGMP.
*   **Oxygen and Metabolism:** Models oxygen consumption in the neuron, ATP production, and its effects on membrane pump activity.
*   **Perivascular and Synaptic Compartments:** Accounts for ion concentrations and transport in the PVS and synaptic cleft.
*   **Vascular Dynamics:** Simulates the smooth muscle cell response, and vascular wall mechanics based on intracellular calcium and cGMP.
*   **BOLD Response:** Calculates the BOLD signal based on oxygen consumption and vessel radius.
*   **Flexible Parameterization:** Allows modification of parameters for specific scenarios and targeted experimentation.

## Features

*   **Modular and Object-Oriented Design:** The code is structured using MATLAB classes (`Astrocyte`, `Neuron`, `SmoothMuscleCell`, `VesselWall`, and `NVU`) to encapsulate model equations and parameters, facilitating extensibility and code organization.
*   **Detailed Multi-Compartment Modeling:** Models the fluxes and concentrations of several key ions across multiple cellular compartments, including the soma, dendrite, synaptic cleft, astrocyte, perivascular space and smooth muscle cell.
*   **Coupled Systems:** Implements interactions between cellular components using shared variables and functions.
*   **Output Analysis:** Provides a variety of outputs for assessing model behavior, such as membrane voltage, ion concentrations, fluxes, vessel radius and BOLD response.
*   **Flexible Time and State Variables:** Allows for simulation time specifications, initial conditions, and integration options.
*   **Parameter Customization:**  Includes multiple parameters allowing the user to fine-tune the desired outputs of each cell type.
*   **Neurovascular coupling:** Models the effect of neuronal activity on the blood oxygenation, vasodilation, and cerebral blood flow (CBF).

## How to Use

1.  **MATLAB Environment:** Ensure you have MATLAB installed on your system.
2.  **Copy the Code:** Copy the contents of the `Astrocyte.m`, `Neuron.m`, `SmoothMuscleCell.m`, `VesselWall.m`, and `NVU.m` files into separate MATLAB script files.
3.  **Run the Simulation:**
    *   Create instances of each class, `Astrocyte`, `Neuron`, `SmoothMuscleCell` and `VesselWall`.
    *   Instantiate the `NVU` class by passing these objects to the constructor.
    *   Call the `simulate()` method of the `NVU` object to run the simulation.
    *   Access simulation results using the `out()` method of the `NVU` object.

```matlab
% Example of how to run the simulation.
% Add the path to where the .m files are saved
addpath('path/to/your/files');

% Create instances of the model components
astro = Astrocyte();
neuron = Neuron();
smcec = SmoothMuscleCell();
wall = VesselWall();

% Create an instance of the NVU class
nvu = NVU(neuron, astro, wall, smcec, 'T', linspace(0, 1200, 1000));

% Run the simulation
nvu.simulate();

% Example of how to plot the membrane voltage of the astrocyte
time = nvu.T;
v_k = nvu.out('v_k');
plot(time, v_k)
xlabel('Time (s)');
ylabel('Membrane Voltage (V)');

% To access the outputs of each sub-model:

% Neuron outputs:
% example
% CBF = nvu.out('CBF');

% Astrocyte outputs:
% example
% J_N_BK_k = nvu.out('J_N_BK_k');

% SMC outputs:
% example
% R_smc = nvu.out('R_smc');

% Wall outputs:
% example
% R = nvu.out('R');
Use code with caution.
Markdown
*   **Optional Parameters:** You can optionally set model parameters when you call the `Astrocyte`, `Neuron`, `SmoothMuscleCell`, `VesselWall` or `NVU` classes. For example: `astro = Astrocyte('g_K_k', 30);` or `nvu = NVU(neuron, astro, wall, smcec, 'T', linspace(0, 100, 1000));`.

*   **Manual Initial Conditions:** It is possible to run the model with initial conditions of your choice by calling the `simulateManualICs` function, for this you must change the `u0` parameter of each of the models.
Use code with caution.
Analyze Results: The simulation returns a time vector, an array of concentrations, and arrays of different model outputs such as membrane voltage, ion concentrations, fluxes, and BOLD signal. You can use MATLAB plotting functions or use the varnames() function of each of the model classes, to extract the required parameters using the out() method from NVU class.

Structure
The project is organized as follows:

Astrocyte.m: Defines the Astrocyte class, which models the astrocyte cell dynamics and includes the dynamics of the perivascular space and synaptic cleft.

Neuron.m: Defines the Neuron class, which models the neuron dynamics.

SmoothMuscleCell.m: Defines the SmoothMuscleCell class, which models the dynamics of the smooth muscle cell in cerebral blood vessels, including key pathways, such as calcium and NO signalling.

WallMechanics.m: Defines the VesselWall class, which models vessel radius dynamics based on intracellular calcium levels.

NVU.m: Defines the NVU class, which integrates the other model components and manages the simulation.

NVUrunScript.m: Demonstrates how to set up and run a simulation using the NVU model, as well as how to access the results.

PLOTallVariables.m: Demonstrates how to plot all variables within each of the sub-models.

Assumptions and Limitations
The model makes simplifying assumptions about cell geometries, diffusion, and reaction rates.

The model does not explicitly account for all cell types and physiological processes present in a neurovascular unit.

The model is computationally intensive, and simulation times can be long depending on the complexity of the simulation.

The model uses ordinary differential equations, which assume a well-mixed environment within the compartments of the model.

Parameter values are based on existing literature and may not be universally applicable to all physiological conditions.

Future Improvements
Incorporate additional cell types, such as pericytes and microglia.

Model the spatial distribution of ions and signaling molecules.

Include more detailed models of cell metabolism and energy production.

Integrate experimental data to validate the model and refine parameters.

Implement sensitivity analysis to evaluate the impact of different parameters on model outputs.

Improve the computation speed by vectorising the sub-models and using mex files.

Add code to estimate the uncertainty in parameter values by fitting the model to experimental data

Technologies Used
MATLAB: Programming language and environment for numerical computation and simulation.

Contact Information
Feel free to reach out if you have any questions or suggestions regarding this model. Email:elshinjo.88@gmail.com