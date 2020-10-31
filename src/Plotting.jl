# Plotting routines -- uses PyPlot

using PyPlot

# include("Utilities.jl");
# include("Response.jl");

# Bodemag with frsp and om given
function BodeMagPlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    plot(om,abs.(freqResponse),linewidth=.75);
    xscale("log");
    yscale("log");
    xlabel("Frequency (rad/s)");
    ylabel("Magnitude");
    tight_layout();
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
    BodeMagPlot(freqResponse,om);
end

# Bode phase with frsp and om given -- this method should be sufficient. Other methods are not necessary.
function BodePhasePlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    phase = 180*atan.(imag.(freqResponse),real.(freqResponse))/pi; # Convert it to [-180,180]
    plot(om,phase,linewidth=.75);
    xscale("log");
    #yscale("log");
    xlabel("Frequency (rad/s)");
    ylabel("Phase");
    tight_layout();
end

function BodePlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    subplot(211); BodeMagPlot(freqResponse,om);
    subplot(212); BodePhasePlot(freqResponse,om);
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
    plot(real.(f1),imag.(f1),linewidth=0.5); # Positive ω
    plot(real.(f2),imag.(f2),linewidth=0.5); # Negative ω
    grid(which="both",color=[1,1,1]*0.9,linewidth=0.25)
    xlabel("Real");
    ylabel("Imaginary");
    tight_layout();
end

function PoleZeroMap()
    # Todo ...
end
