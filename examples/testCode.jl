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
p1 = plot(RobustControl.BodeMagPlot(G,om),title="G1"); 
p2 = plot(RobustControl.BodeMagPlot(G,om),title="G2"); 
display(plot(p1,p2));


# Nyquist Plot
p = RobustControl.NyquistPlot(s,1,1);
display(p);

# # System Responses
# figure(4); clf();
# subplot(2,2,1);
# t,y = RobustControl.ImpulseResponse(s,1,1);
# plot(t,y); grid();
# xlabel("Time (sec)");
# #ylabel("Impulse Response");
# title("Impulse Response");tight_layout();

# subplot(2,2,2);
# t,y = RobustControl.StepResponse(s,1,1);
# plot(t,y); grid();
# xlabel("Time (sec)");
# #ylabel("Step Response");
# title("Step Response");tight_layout();

# subplot(2,2,3);
# t = Vector(0:0.01:50);
# u = ones(length(t)) + 0.1*randn(Float64,length(t)); # Step + Disturbance
# x0 = [0.0;0.0];
# y = RobustControl.ForcedResponse(s,1,1,t,u,x0);
# plot(t,y); grid();
# xlabel("Time (sec)");
# #ylabel("Time Response to Step+Disturbance");
# title("Forced Response");tight_layout();

# subplot(2,2,4);
# t = Vector(0:0.01:10);
# x0 = [1.0;0.0];
# t1,y1 = RobustControl.InitialResponse(s,x0); # Auto time grid
# y2 = RobustControl.InitialResponse(s,x0,t); # Specified time grid
# l1=plot(t1,y1,"b",linewidth=1);
# l2=plot(t,y2,"r",linestyle=":",linewidth=2);
# xlabel("Time (sec)");
# #ylabel("Initial Response");
# legend([l1[1],l2[1]],["(auto time grid)","(specified time grid)"]);
# title("Initial Response");
# grid(); tight_layout();
