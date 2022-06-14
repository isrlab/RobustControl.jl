# Plotting routines -- uses Plots. Not using PyPlot Directly.
# Bodemag with frsp and om given

function BodeMagPlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    p = plot(om, abs.(freqResponse), xaxis=:log, yaxis=:log, label=:none, xlabel="Frequency (rad/s)", ylabel = "Magnitude");
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
    p = plot(om, phase, xaxis=:log, yaxis=:log, label=:none, xlabel="Frequency (rad/s)", ylabel = "Phase (deg)");
    return p;
end

function BodePlot(freqResponse::Vector{ComplexF64},om::Vector{Float64})
    p1 = BodeMagPlot(freqResponse,om);
    p2 = BodePhasePlot(freqResponse,om);
    return p1,p2;
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
    return BodePlot(freqResponse,om);
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
    p = plot(real.(f1),imag.(f1),color=:red,label="Positive"); # Positive ω
    plot!(real.(f2),imag.(f2),label="Negative"); # Negative ω
    plot!(xlabel="Real", ylabel = "Imaginary");
    scatter!([real(f1[1])],[imag(f1[1])],color=:white,label=L"\omega = 0")
    scatter!([real(f1[end])],[imag(f1[end])],color=:black,label=L"\omega = \pm \infty")        
    return p
end

function PoleZeroMap()
    # Todo ...
end
