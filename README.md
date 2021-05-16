# Lagrangian stochastic model for fire propagation in forests and the forest/urban interface

## Introduction
Large scale fires in forests and cities are a serious hazard for the environment and human population. In order to adopt proper measures to prevent fires or take actions to reduce the effects of a fire, it is very important to have tools able to compute fire propagation and ignition probability in short time. In this project, an ignition probability model developed and validated in the context of gas-turbine combustors will be extended to the prediction of ignition probability in large-scale fires. The model is based on the stochastic tracking of ”flame kernels” and we’ll improve the model by using data-driven method and machine learning. A specific fire scenario will be used as a case of study to extract data for data assimilation methods.

Data assimilation is used to optimally combine a surrogate model withs parse noisy data ie combine different sources of information to estimate possible states of a system as it evolves in time. Data assimilation is useful here since our model is sensitive to many parameters and can be improved by using the good parameters. The DA technique we’ll be using is the ensemble Kalman filter and there is a lot of documentation/articles that help our understanding of DA in this project (see references).

## 1. The Lagrangian stochastic model for fire propagation
A forest fire in general propagates due to various physical processes: **preheating** and **pyrolysis/devolatilisation**, **radiation**, **convection**, **embers** and **firebrand**, **effect of inclination**...

In order to model the fire propagation, we need to consider these processes and model them first. The model we will be using borrows ideas from the [cellular automata](https://en.wikipedia.org/wiki/Cellular_automaton) approach.
The first thing to do is the discretisation of space. We split the region (to be studied) in cells. Therefore, a cell is a **flammable** or **non-flammable material** and has at least 3 caracteristics:
- *ignition time delay*: time difference between initiation of the fire (e.g. by an ember) and self-sustaining burning.
- *maximum expected heat release*
- *burn duration* : total expected duration of the fire.


These parameters are pre-determined by examination of each cell.

The next step is to simulate the fire propagation process. The emission of the convection particles and radiation particles are modelled as follow:
- **convection**: A turbulent random walk of virtual convection particles, which are emitted by each individual cell and have 3 caracteristics: a velocity vector `Up`, a position vector `Xp` and a scalar `Yst,p` (between 0 and 1) denoting the burning status of the particle
- **radiation**: a radiation particle is emitted from a burning cell at a random angle, and a step `Lr`.

To learn more about the lagrangian stochastic model implemented here, read this article : [SPINTHIR](https://imperiallondon-my.sharepoint.com/personal/agiusti_ic_ac_uk/Documents/Microsoft%20Teams%20Chat%20Files/SPINTHIR_v1(1).pdf)

## 2. How to launch the model simulation ?

It's quite simple to launch the simulation. You can either use MATLAB or Python. 
- with MATLAB, you have to execute the file `./SPINTHIR/SPINTHIR_ToyProblems_v5.m`
- with python, I've implemented the same code in python version so you can execute `./SPINTHIR/SPINTHIR_ToyProblems_v5.py`.
However, a best implementation use python `class` and you can launch the simulation using `./main.py` and setting the necessary parameters.

We have 3 main classes (in the folder `LagrangianModel`) to represent our model objects:
- `Cell`: this class model a cell as described in the previous section.
- `Grid`: this class represent the area and is divided in cells.
- `Particle`: this class represent a single particle which can move from cell to cell.


## 3. What about the data assimilation ?
You can launch the data assimilation by executing the file `./data_assimilation.py`. 

