# Code to test basic functionality and plotting

include("StateSpace.jl");
include("Utilities.jl");
include("Response.jl");

A = [0 1;-1 -1];
B = reshape([0;1],2,1);
C = eye(2);
D = reshape([0;0],2,1);;
s = StateSpace(A,B,C);
rr,ii,om,d = rifd(s);
z = TransmissionZeros(s,2,1);
p = Poles(s);
d = Damping(s);
om = NaturalFrequencies(s);

# Frequency response with auto frequency grid.
G,om = FrequencyResponse(s,1,1);

# Frequency response with specified frequency grid.
om1 = logspace(1E-2,1E2,200);
G1 = FrequencyResponse(s,2,1,om1);

# Test plots
using PyPlot # We will be using PyPlot for all our plotting needs.
include("Plotting.jl"); #close("all")
figure(1);clf();
subplot(2,1,1);
 BodeMagPlot(G,om); title("Bode Mag Plot");
 BodeMagPlot(G1,om1);
grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
legend(["G1","G2"]);
tight_layout();


# Another way to get Bode plots -- just magnitude.
subplot(2,1,2);
BodeMagPlot(s,1,1);
BodeMagPlot(s,2,1);
grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
legend(["G1","G2"]);
tight_layout();

# 3rd way -- get both mag and phase in one plot
figure(2); clf();
BodePlot(s,1,1);
BodePlot(s,2,1);
subplot(211);grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
subplot(212);grid(which="both",color=[1,1,1]*0.9,linewidth=0.35)
subplot(211);legend(["G1","G2"]);
title("Bode Plot");
tight_layout();

# Nyquist Plot
figure(3);clf();
subplot(2,1,1);NyquistPlot(s,1,1); title("Nyquist Plot"); # With auto frequency grid
subplot(2,1,2);NyquistPlot(s,1,1,om1); # With specified frequency grid

# System Responses
figure(4); clf();
subplot(2,2,1);
t,y = ImpulseResponse(s,1,1);
plot(t,y); grid();
xlabel("Time (sec)");
#ylabel("Impulse Response");
title("Impulse Response");tight_layout();

subplot(2,2,2);
t,y = StepResponse(s,1,1);
plot(t,y); grid();
xlabel("Time (sec)");
#ylabel("Step Response");
title("Step Response");tight_layout();

subplot(2,2,3);
t = Vector(0:0.01:50);
u = ones(length(t)) + 0.1*randn(Float64,length(t)); # Step + Disturbance
x0 = [0.0;0.0];
y = ForcedResponse(s,1,1,t,u,x0);
plot(t,y); grid();
xlabel("Time (sec)");
#ylabel("Time Response to Step+Disturbance");
title("Forced Response");tight_layout();

subplot(2,2,4);
t = Vector(0:0.01:10);
x0 = [1.0;0.0];
t1,y1 = InitialResponse(s,x0); # Auto time grid
y2 = InitialResponse(s,x0,t); # Specified time grid
l1=plot(t1,y1,"b",linewidth=1);
l2=plot(t,y2,"r",linestyle=":",linewidth=2);
xlabel("Time (sec)");
#ylabel("Initial Response");
legend([l1[1],l2[1]],["(auto time grid)","(specified time grid)"]);
title("Initial Response");
grid(); tight_layout();
