<img src="https://isrlab.github.io/images/julia%20robust%20control%20library.png" alt="Julia RCS Logo" width="350"/>

This is a robust control toolbox in Julia. The goal is to provide necessary functions for designing H2/Hinf controllers and estimators with model uncertainty.

The synthesis is formulated as a convex optimization problem, which is solved using JuMP/Convex. The default solver is SCS. If Mosek is available, it can be used. Future version will also interface with SLICOT.

Some functions have been implemented. See testCode.jl under src/ to see basic analysis usage. 

This toolbox is not ready for use yet. Please stay tuned. This is under heavy development.

# Installation 
``` julia

Pkg.add("https://github.com/isrlab/RobustControl");

```
# Usage

``` julia
# Code to test basic functionality and plotting
using RobustControl, PyPlot

A = [0 1;-1 -1];
B = reshape([0;1],2,1);
C = RobustControl.eye(2);
D = reshape([0;0],2,1);;
s = RobustControl.StateSpace(A,B,C);
rr,ii,om,d = RobustControl.rifd(s);
z = RobustControl.TransmissionZeros(s,2,1);
p = RobustControl.Poles(s);
d = RobustControl.Damping(s);
om = RobustControl.NaturalFrequencies(s);

# Frequency response with auto frequency grid.
G,om = RobustControl.FrequencyResponse(s,1,1);

# Frequency response with specified frequency grid.
om1 = RobustControl.logspace(1E-2,1E2,200);
G1 = RobustControl.FrequencyResponse(s,2,1,om1);

# Test plots
figure(1);clf(); pygui(true);
subplot(2,1,1);
RobustControl.BodeMagPlot(G,om); title("Bode Mag Plot");
RobustControl.BodeMagPlot(G1,om1);
grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
legend(["G1","G2"]);
tight_layout();

# Another way to get Bode plots -- just magnitude.
subplot(2,1,2);
RobustControl.BodeMagPlot(s,1,1);
RobustControl.BodeMagPlot(s,2,1);
grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
legend(["G1","G2"]);
tight_layout();

# 3rd way -- get both mag and phase in one plot
figure(2); clf();
RobustControl.BodePlot(s,1,1);
RobustControl.BodePlot(s,2,1);
subplot(211);grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
subplot(212);grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
subplot(211);legend(["G1","G2"]);
title("Bode Plot");
tight_layout();

# Nyquist Plot
figure(3);clf();
subplot(2,1,1);RobustControl.NyquistPlot(s,1,1); title("Nyquist Plot"); # With auto frequency grid
subplot(2,1,2);RobustControl.NyquistPlot(s,1,1,om1); # With specified frequency grid

# System Responses
figure(4); clf();
subplot(2,2,1);
t,y = RobustControl.ImpulseResponse(s,1,1);
plot(t,y); grid();
xlabel("Time (sec)");
#ylabel("Impulse Response");
title("Impulse Response");tight_layout();

subplot(2,2,2);
t,y = RobustControl.StepResponse(s,1,1);
plot(t,y); grid();
xlabel("Time (sec)");
#ylabel("Step Response");
title("Step Response");tight_layout();

subplot(2,2,3);
t = Vector(0:0.01:50);
u = ones(length(t)) + 0.1*randn(Float64,length(t)); # Step + Disturbance
x0 = [0.0;0.0];
y = RobustControl.ForcedResponse(s,1,1,t,u,x0);
plot(t,y); grid();
xlabel("Time (sec)");
#ylabel("Time Response to Step+Disturbance");
title("Forced Response");tight_layout();

subplot(2,2,4);
t = Vector(0:0.01:10);
x0 = [1.0;0.0];
t1,y1 = RobustControl.InitialResponse(s,x0); # Auto time grid
y2 = RobustControl.InitialResponse(s,x0,t); # Specified time grid
l1=plot(t1,y1,"b",linewidth=1);
l2=plot(t,y2,"r",linestyle=":",linewidth=2);
xlabel("Time (sec)");
#ylabel("Initial Response");
legend([l1[1],l2[1]],["(auto time grid)","(specified time grid)"]);
title("Initial Response");
grid(); tight_layout();
```
