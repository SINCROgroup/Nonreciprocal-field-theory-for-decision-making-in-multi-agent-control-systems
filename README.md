# Nonreciprocal-field-theory-for-decision-making-in-multi-agent-control-systems

This folder contains the code used to generate the data for the results presented in [LAMA2025]. All code is written in MATLAB, version 2022A.

In [LAMA2025], the authors address the shepherding control problem, where a group of agents, the herders, must coordinate to steer the collective dynamics of another group of agents, the targets. Specifically, the scenario considered is one in which the herders must guide the targets to a prescribed goal region and contain them within it. It is assumed that herders repel the targets, which enables them to complete the task.

The herders' strategy involves two key decision-making components:
i) selecting one target to chase, based on its distance from the goal region, and
ii) steering the selected target toward the goal region by positioning themselves "behind" it.


The paper presents numerical results at two levels of description:

-the Agent-Based (AB) level, where the shepherding task is simulated in two dimensions (2D) within a periodic box, and

-the Continuum level, where the authors derive and simulate a "shepherding" field theory consisting of two coupled partial differential equations (PDEs), simulated on a 1D periodic domain.

All simulations can be launched from the main.m script. This file introduces and explains all parameters, and its variable naming is kept as consistent as possible with those used in [LAMA2025].
Below is a list of the functions in this folder, along with brief descriptions of their roles.

- main.m: Initializes and explains all parameters. From here, you can run both AB and PDE simulations. Simulation data is saved upon completion and also visualized during the simulation. The script also includes code for replaying saved simulations.

AB simulations*
- AB_radial.m:      Simulates the agent-based equations (Eqs. 7–8) to reproduce Figs. 3, 4(a–b), S1, S3. This case considers a circular goal region centered at the origin.
- AB_rectangular:   Simulates the agent-based equations (Eqs. 7–8, S2, S38, S41) to reproduce Figs. 5(i–j), S2. This case involves a rectangular goal region, where the objective is to either collect or expel the targets.
- AB_static_patterns: Simulates the agent-based equations (Eqs. 7–8, S43) to reproduce Fig. 5(k). In this case, the goal is to gather the targets into multiple equally spaced rectangular regions.
- AB_travelling_patterns: Simulates the agent-based equations (Eq. S45) to reproduce Fig. 5(l). The goal is to create a traveling wave state in which the herders push the targets along the periodic domain.

*PDEs simulations*
- PDE_simulation:  Simulates the PDEs describing the shepherding dynamics. From main.m, you can set the type of dynamics via the input variable type, which determines the forms of the coupling functions v1(x) and v2(x) (see [LAMA2025], Eqs (3-4)). The available types are:
  
- a) type ="main"; v1(x) and v2(x) are chosen as of Eqs. (15-16); needed to reproduce Figs (2, 4(c-d), S5(b)).
- b) type="containment"; v2 is a positive constant and v1 is a square wave with  2\pi/L period (Eq. S34); needed to reproduce Fig (5(e)).
- c) type="expulsion"; v2 is a positive constant and v1 is a square wave (with a minus sign) with  2\pi/L period (Eq. S35); needed to reproduce Fig (5(f)).
- d) type="static_patterns"; v2 is a positive constant and v1 is a square wave with  6\pi/L period (Eq. S36); needed to reproduce Fig (5(g)).
- e) type="travelling_patterns"; v1 and v2 are constants (with same sign) (Eq. S37); needed to reproduce Fig (5(h)).

The other functions in the folder are "auxiliary" functions needed to complete the AB simulations
-attraction.m: computes the attraction of herders to targets; depending on parameters, this function can include decision-making logic (e.g., target selection and approach strategy)
-cartesian_to_polar.m: converts an array of cartesian coordinates to an array of polar coordinates
-intial_pos.m: randomly generates uniformly distributed initial positions within a box.
-minimum_image_distance.m: compute relative positions and distances considering the periodicity of the domain.
-periodic.m: prepares an array of positions, defined in a periodic domain, for plotting
-repulsion.m: computes the repulsion between two (sets of) agents
