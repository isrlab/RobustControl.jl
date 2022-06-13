# Plotting routines -- uses Plots. Not using PyPlot Directly.

using Plots
# Plot defaults
default(titlefont = (10, "times"), legendfontsize = 10, guidefont = (8, :black), 
tickfont = (6, :black), framestyle = :default, minorticks=:false , label=:none, 
yminorgrid = false, xlabel = "", ylabel = "", linewidth=1, color=:dodgerblue2);

# Bodemag with frsp and om given
function BodeMagPlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    p = plot(om, abs.(freqResponse), xaxis=:log, yaxis=:log, xlabel="Frequency (rad/s)", ylabel = "Magnitude");
    display(p);
    return p
end

# Bodemag directly from system
function BodeMagPlot(sys::StateSpace,no::Integer, ni::Integer, om::Vector{Float64}=Array{Float64}(undef,0))
    if (ni>sys.nu) error("Input index exceeds number of system inputs") end
    if (no>sys.ny) error("Output index exceeds number of system output") end

    if isempty(om)
        freqResponse,om = FrequencyResponse(sys,no,ni);
    else
        freqResponse = FrequencyResponse(sys,no,ni,om);
    end
    return BodeMagPlot(freqResponse,om);
end

# Bode phase with frsp and om given -- this method should be sufficient. Other methods are not necessary.
function BodePhasePlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    phase = 180*atan.(imag.(freqResponse),real.(freqResponse))/pi; # Convert it to [-180,180]
    p = plot(om, phase, xaxis=:log, yaxis=:log, xlabel="Frequency (rad/s)", ylabel = "Phase (deg)");
    display(p);
    return p;
end

function BodePlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    p1 = BodeMagPlot(freqResponse,om);
    p2 = BodePhasePlot(freqResponse,om);
    p = plot(p1,p2); # Need to check dimension and layout.
    display(p);
    return p;
end

# Bodemag directly from system
function BodePlot(sys::StateSpace,no::Integer, ni::Integer, om::Vector{Float64}=Array{Float64}(undef,0))
    if (ni>sys.nu) error("Input index exceeds number of system inputs") end
    if (no>sys.ny) error("Output index exceeds number of system output") end

    if isempty(om)
        freqResponse,om = FrequencyResponse(sys,no,ni);
    else
        freqResponse = FrequencyResponse(sys,no,ni,om);
    end
    BodePlot(freqResponse,om);
end

function NyquistPlot(sys::StateSpace,no::Integer, ni::Integer, om::Vector{Float64}=Array{Float64}(undef,0))
    if (ni>sys.nu) error("Input index exceeds number of system inputs") end
    if (no>sys.ny) error("Output index exceeds number of system output") end

    if isempty(om)
        f1,om = FrequencyResponse(sys,no,ni);
        f2 = FrequencyResponse(sys,no,ni,-om);
    else
        f1 = FrequencyResponse(sys,no,ni,om);
        f2 = FrequencyResponse(sys,no,ni,-om);
    end
    p = plot(real.(f1),imag.(f1),color=:red); # Positive ω
    plot!(real.(f2),imag.(f2)); # Negative ω
    plot!(grid=:true,grid_lw=0.5, xlabel="Real", "Imaginary");
    display(p);
    return p
end

function PoleZeroMap()
    # Todo ...
end
