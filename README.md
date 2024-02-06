# 2D KSE and calmed KSE solver (Matlab)

This repository contains a 2D Kuramoto-Sivashinsky Equation (KSE) solver implemented in Matlab 
which simulates both solutions to KSE and a modified set of equations, called the calmed Kuramoto-Sivashinsky equations.
Calmed KSE is a modification of 2D KSE introduced in [this body of work](https://arxiv.org/abs/2304.10493). 

## Features

- Efficient 2D KSE solver using pseudo-spectral methods
- Uses IMEX time-stepping conditions
- Selects time-stepping to respect CFL conditions
- Customizable simulation settings
- Generates convergence plots to show that calmed KSE converges to KSE as calming parameter tends to 0

## Getting Started

### Prerequisites

To run the solver, you need to have Matlab installed on your system. You can download Matlab [here](https://www.mathworks.com/products/matlab.html).

### Installation

1. Clone the repository.

2. Navigate to the project folder in Matlab.


## Usage

### Running simulations
- To run a simulation of both KSE and calmed KSE, run 'KSE2D_IF_RK4_calm.m'. 
- You can configure the simulation parameters by modifying the parameters and default settings at the start of 'KSE2D_IF_RK4_calm.m'.
- To increase the number of linearly unstable Fourier modes, increase the value of lambda in default settings.
- The calming parameter is determined by p.epsilon. For calmed KSE to better approximate KSE, reduce the size of epsilon in default settings.
- To determine the type of calming function used, set type = 1, 2, or 3 in default settings.
- To see an updated figure display, set setting.display_on = 1 in the settings section. Otherwise the figure will not update.
- To set the figure to display every h units of time, set setting.update_time = h. On default, it will display the figure 100 times.

### Testing convergence
- To perform a convergence analysis, run 'calm_convergence_scriptCaller.m'. 
- For this script to run properly, set setting.convergence_on = 1 in script 'KSE2D_IF_RK4_calm.m'. This is the default setting. 
- On default, the calming paramater will be set within the range 1e-14 to 1e-4, but this can easily be adjusted for personal use. 


