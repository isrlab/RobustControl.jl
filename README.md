<img src="https://isrlab.github.io/assets/images/julia%20robust%20control%20library.png" alt="Julia RCS Logo" width="350"/>

This is a robust control toolbox in Julia. The goal is to provide necessary functions for designing H2/Hinf controllers and estimators with model uncertainty.

The synthesis is formulated as a convex optimization problem, which is solved using JuMP/Convex. Future version will also interface with SLICOT.

Some functions have been implemented. See examples/testCode.jl for basic analysis usage.

<!-- This toolbox is not ready for use yet. Please stay tuned. This is under heavy development. -->

# Installation

``` julia

Pkg.add("https://github.com/isrlab/RobustControl");

```

# Usage

``` julia
# Code to test basic functionality and plotting
using RobustControl, Plots, LaTeXStrings

default(titlefontsize = 10, legendfontsize = 8, labelfontsize=8, guidefont = (8, :black), 
tickfont = (6, :black), framestyle = :default, minorticks=:true , label=:none, 
yminorgrid = true, xminorgrid=true, xlabel = "", ylabel = "", linewidth=1, color=:dodgerblue2);


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
G1 = RobustControl.FrequencyResponse(s,1,1,om1);

# Test plots
p1 = plot(RobustControl.BodeMagPlot(G,om),title="Bode Magnitude Plot of G(1,1)"); 
display(p1);

# Nyquist Plot
p = plot(RobustControl.NyquistPlot(s,1,1),title="Nyquist Plot of G(1,1)");
display(p);

# System Responses
t,y = RobustControl.ImpulseResponse(s,1,1);
p1=plot(t,y,xlabel="Time (s)", title = "Impulse Response",label=:none); 

t,y = RobustControl.StepResponse(s,1,1);
p2=plot(t,y,xlabel="Time (s)", title = "Step Response",label=:none); 

t = Vector(0:0.01:50);
u = ones(length(t)) + 0.1*randn(Float64,length(t)); # Step + Disturbance
x0 = [0.0;0.0];
y = RobustControl.ForcedResponse(s,1,1,t,u,x0);
p3 = plot(t,y,xlabel="Time (s)", title = "Forced Response",label=:none); 

x0 = [1.0;0.0];
# t,y = RobustControl.InitialResponse(s,x0); # Auto time grid
t = Vector(0:0.01:10);
y = RobustControl.InitialResponse(s,x0,t); # Specified time grid
p4 = plot(t,y,xlabel="Time (s)", layout=(1,2), 
     label = [L"x_1" L"x_2"],title = "IC Response", yminorgrid=:false,
     legend=[:topright :bottomright]); 

p = plot(p1,p2,p3,p4,layout = (2,2));
display(p)
```

## Output of the above code

![](images/BodeMag.png)

![](images/Nyquist.png)

![](images/VariousResponses.png)
