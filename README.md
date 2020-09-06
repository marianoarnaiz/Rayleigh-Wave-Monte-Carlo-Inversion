# Rayleigh-Wave-Monte-Carlo-Inversion
A code to jointly invert dispersion curves for R waves

MC_R_CyU_VariableThickness.m

This code performs Monte Carlo inversion (Markov Chain MC) of Rayleigh wave dispersion curves (Phase and Group jointly). It searches for a 1D high probability model that fits the observed data. The code uses surf96 from Herrmann’s Seismology codes (http://www.eas.slu.edu/eqc/eqccps.html). These set of codes should be installed before attempting to run the inversion. 

As any MC code, running this takes time! Speed depends on your computer and Matlab version. Inversion might take up to a day depending on the number of iterations, but experience suggest that 100000 iterations are OK to solve a problem. Be patient

Please read: Bosch (1999; 2001 and 2005) and Mosegaard & Sambridge (2002) for the method and background. 


The code begins from a model (the layers to modify and the number of layers to be considered should be set as explained) and it is perturbed to find the solution. This is done by choosing a layer or set of layers (randomly) and changing the velocity or thickness of the layers. Layers are not allowed to appear or disappear.  If a “perturbation” reduces de error it is accepted, otherwise it has a probability of been accepted or rejected. If rejected the previous model is repeated.

Parameters:

IT: Number of iterations
Model0: Initial model (a constant “uninformative mode” is fine). The file is composed of a set of constant thickness and velocities layers followed by the AK135.
R_Cobs: Observed phase velocity. Use the same period range as R_Uobs
R_Uobs: Observed group velocity. Use the same period range as R_Cobs
VpVs:Vp/Vs ratio.
maxVstop: Max Vs to accept in the candidate models on the top (crust?)
minVstop: Min Vs to accept in the candidate models on the top (crust?)
maxVsbottom: Max Vs to accept in the candidate models on the bottom (mantle?)
minVsbottom: Min Vs to accept in the candidate models on the bottom (mantle?)
fv: Initial Variance.
fv2: reduce variance at some point (this makes harder for new models to be accepted).
a = 2; This is the number of the first layer in the range to accept perturbations (there is a 0 km thick layers in the model, this prevents problems when the forward calculation is performed)
b = 11 This is the number of the last layer to in the range to accept perturbations (i.e. layer 12 won’t be chanced).
sl=5; Step layer. This layer defines where the top/bottom limit for the search is.
% Random amplitude of the assymetric perturbation
c = 0; Minimum number of layers to perturb. 1 means, perturb at least 1 layer. 2 means perturb a least 2. Could be 0 to allow for stronger variations.
d = 2; Half of the max number of layers to perturb.
smoothfactor=0.25 Smoothing parameter for the Vs profile. This avoids extreme Vs jumps. 0.25 works well. Recommended between 0.1 and 0.5.

In line 261 the code reads another file. This is the solution to the example data. You can comment this line and change the plot. 


Tips:

Make a test run before changing some parameters.
Change some parameters to see what happens!
Any questions, email me!.
