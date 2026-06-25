
include("Aux.jl")
import Plots
import CairoMakie


# TODO: Make sure that the documentation for structs work?
"""
    @kwdef mutable struct SequenceBlock
    
    Creates an event block. Events can be an `RF` event, `Delay` event, (gradient) `Spoiling` event, or `Ignore`, in which case the event will be skipped.
    
    # Fields
    `RF`: Complex RF signal [μT]
    `T`: Time array with the same size as RF [s]
    `Δω`: Frequency offset of RF [Hz]
    `ϕ`: Phase offset of RF [rad]
    `EventType`: `:RF`, `:Delay`, `:Spoiling`, `:ADC`, or `:Ignore`
    `EventName`: Name of the event e.g., "Sat RF"

    If an event is coded as `:Spoiling`, Gradient spoiling is assumed; For RF spoiling, please create an `RF` event
"""

@with_kw mutable struct SequenceBlock
    # Not to be conflicted with KomaMRI's Sequence() mutable struct
    # 3 types of EVentName accepted: RF, Delay, and Spoiling
    RF::Union{Array{ComplexF64,1}, ComplexF64, Nothing}=nothing
    T::Union{Array{Float64,1}, Nothing}=[0.0]
    Δω::Union{Float64, Nothing}=nothing
    ϕ::Union{Float64, Nothing}=nothing
    EventType::Symbol=:Ignore
    EventName::String=""
end

@with_kw mutable struct CESTSequence
    System::Union{NamedTuple, CESTSys}=defaultSystem()
    Events::Array{SequenceBlock,1}=[SequenceBlock(EventName="Seq Init")]
    Author::String=""
    Name::String=""
end

"""
    CoreBlocks

    Mutable Struct containing core `SequenceBlock`(s) to be used for simulations. The usual event flows can be, but are not limited to, the following:
        Sat -> (SatSpoil) -> SatTd -> SatImgSpoil -> (SatImgDelay) -> Img -> (ImgSpoil) -> ImgTE -> (ImgTR) -> Recovery
        
    # Fields 
        `Sat`: Individual Saturation RF event, also known as pulse duration, i.e., Sat RF On
        `SatSpoil`: Spoiling event between each Saturation RF; Setting Mxy = 0
        `SatTd`: Delay between saturation RF events, i.e., Sat RF OFF

        `SatImgSpoil`: Spoiling between last Sat RF and first Img RF events
        `SatImgDelay`: Delay between last SatRF (or SatImgSpoil) and first Img RF events
        
        `Img`: Imaging RF event; i.e., Img RF On
        `ImgSpoil`: Spoiling event after each Img RF event
        `ImgTE`: Delay between Img RF event and start of Readout
        `ImgTR`: Delay between each Img RF event
        `ImgTRO`: ADC on

        `Recovery`: Delay between each Sat-Img block

"""

@with_kw mutable struct CoreBlocks
    Sat::SequenceBlock # Individual saturation RF event duration
    SatSpoil::SequenceBlock=SequenceBlock(EventType=:Ignore, EventName="") 
    SatTd::SequenceBlock # Delay between two consecutive saturation RF event
    
    SatImgSpoil::SequenceBlock
    SatImgDelay::SequenceBlock=SequenceBlock(EventType=:Ignore, EventName="") # delay between last saturation RF event and first imaging RF event
    
    Img::SequenceBlock=SequenceBlock(EventType=:Ignore, EventName="") # Individual imaging RF event duration
    ImgSpoil::SequenceBlock=SequenceBlock(EventType=:Ignore, EventName="")
    ImgTE::SequenceBlock=SequenceBlock(EventType=:Ignore, EventName="") # Time to Echo
    ImgTR::SequenceBlock=SequenceBlock(EventType=:Ignore, EventName="") # Repetition Time
    ImgTRO::SequenceBlock # Readout Time
    
    Recovery::SequenceBlock # recovery time after each full block 
end

"""
    function addBlock!(seq::CESTSequence, block::SequenceBlock; roundingdigits::Union{Nothing, Int64}=nothing)

    In-place append `block` to existing `seq`

    # Arguments
        - seq - of type `CESTSequence`
        - block - of type `SequenceBlock`

    # Returns
        - seq - modified sequence struct of type `SequenceBlock`(s)
"""
function addBlock!(seq::CESTSequence, block::SequenceBlock; roundingdigits::Union{Nothing, Int64}=nothing)
    # isnothing(seq.Events[end].T) ? τend = 0.0 : τend = seq.Events[end].T[end] # redundant with the new Sequence definition
    newblock = deepcopy(block)
    if !isnothing(roundingdigits)
        newblock.T = round.(block.T .+ seq.Events[end].T[end], digits=roundingdigits)
    else
        newblock.T = block.T .+ seq.Events[end].T[end]
    end
    seq.Events = [seq.Events; newblock]
    return seq
end


"""
    function removeBlock!(seq::CESTSequence, blocknum::Int64; doublecheck::Bool=true)

    In-place remove `block` from `seq`

    # Arguments
        - seq - of type `CESTSequence`
        - blocknum - the index of the Event block to be removed, of type `Int64`
        - block - of type `SequenceBlock`
        - (optional) doublecheck - whether to check the name and type of Event block before removal; of type `Bool`

    # Returns
        - seq - modified sequence struct of type `SequenceBlock`(s)
"""
function removeBlock!(seq::CESTSequence, blocknum::Int64; doublecheck::Bool=true)
    if blocknum > 1
        eventname = seq.Events[blocknum].EventName
        eventtype = seq.Events[blocknum].EventType
        if doublecheck
            println("=> Block number `" * string(blocknum) * "` has name `" * eventname * "`, of type `" * string(eventtype) * "`.")
            println("==> Are you sure you want to delete this block?")
            println("===> Type `Y` for Yes; any other key for No")
            print("====> ")
            checker = readline()
            if checker === "Y"
                # Recalculate time vectors after this event block
                if blocknum === 2
                    δτ = seq.Events[blocknum].T[end]
                else
                    δτ = seq.Events[blocknum].T[end] - seq.Events[blocknum-1].T[end]
                end
                deleteat!(seq.Events, blocknum)
                for n = blocknum : length(seq.Events)
                    seq.Events[n].T = seq.Events[n].T .- δτ
                end
                println("=====> Block number " * string(blocknum) * " deleted.")
            else
                println("====> No changes to sequence")
            end
        else
            println("=> Deleted block number `" * string(blocknum) * "`, name `" * eventname * "`, of type `" * string(eventtype) * "`.")
            
            # Recalculate time vectors after this event block
            if blocknum === 2
                δτ = seq.Events[blocknum].T[end] # this is now redundant due to the new seq definition, oops!
            else
                δτ = seq.Events[blocknum].T[end] - seq.Events[blocknum-1].T[end]
            end
            deleteat!(seq.Events, blocknum)
            for n = blocknum : length(seq.Events)
                seq.Events[n].T = seq.Events[n].T .- δτ
            end

        end
    else
        throw(ArgumentError("Cannot remove Event `Block 1` (Initialisation block)"))
    end
    return seq
end

"""
    function insertBlock!(seq::CESTSequence, block::SequenceBlock, blocknum::Int64)

    In-place insert `block` into `seq` at `blocknum`

    # Arguments
        - seq - of type `CESTSequence`
        - block - of type `SequenceBlock`
        - blocknum - the index of the Event block to be inserted at, of type `Int64`
        
    # Returns
        - seq - modified sequence struct of type `SequenceBlock`(s)
"""
function insertBlock!(seq::CESTSequence, block::SequenceBlock, blocknum::Int64)
    if blocknum > 1
        # δτ = seq.Events[blocknum-1].T[end] + block.T[end] - seq.Events[blocknum].T[end]
        newblock = deepcopy(block)
        newblock.T = block.T .+ seq.Events[blocknum-1].T[end]

        for n = blocknum:length(seq.Events)
            seq.Events[n].T = seq.Events[n].T .+ block.T[end]
        end

        insert!(seq.Events, blocknum, newblock)

    else
        throw(ArgumentError("Cannot insert to Event `Block 1` (Initialisation block)"))
    end
end

"""
    function replaceBlock!(seq::CESTSequence, block::SequenceBlock, blocknum::Int64; doublecheck::Bool=true)

    In-place replace Event at `blocknum` in `seq` by `block` 

    # Arguments
        - seq - of type `CESTSequence`
        - block - of type `SequenceBlock`
        - blocknum - the index of the Event block to be replaced, of type `Int64`
        - (optional) doublecheck - whether to check the name and type of Event block before removal; of type `Bool`
        
    # Returns
        - seq - modified sequence struct of type `SequenceBlock`(s)
"""
function replaceBlock!(seq::CESTSequence, block::SequenceBlock, blocknum::Int64; doublecheck::Bool=true)
    if blocknum > 1
        removeBlock!(seq, blocknum; doublecheck=doublecheck)
        insertBlock!(seq, block, blocknum)
    else
        throw(ArgumentError("Cannot replace Event `Block 1` (Initialisation block)"))
    end
end


function plot_CESTseq(seq::CESTSequence, fast::Bool, kwargs...)

    f = Plots.plot(size=(500,400), title=seq.Name, y_foreground_color_border=:red)
    f2 = Plots.twinx(f)

    for n = eachindex(seq.Events)
        if seq.Events[n].EventType === :RF
            if isa(seq.Events[n].RF, ComplexF64)
                rf = convert(Array{ComplexF64,1}, repeat([seq.Events[n].RF], length(seq.Events[n].T)))
                Plots.plot!(f, seq.Events[n].T , abs.(rf), label="",color=:red)
                Plots.plot!(f2, seq.Events[n].T, angle.(rf),label="", color=:orange, linestyle=:dot, y_foreground_color_border=:orange)
            else 
                Plots.plot!(f, seq.Events[n].T, abs.(seq.Events[n].RF), label="",color=:red)
                Plots.plot!(f2, seq.Events[n].T, angle.(seq.Events[n].RF),label="", color=:orange, linestyle=:dot, y_foreground_color_border=:orange)
            end
        elseif seq.Events[n].EventType === :Delay
            Plots.plot!(f, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", color=:gray)
        end

    end
    Plots.ylabel!(f, "|B₁| [μT]")
    Plots.ylabel!(f2, "ϕ [rad]")
    Plots.xlabel!(f, "Time [s]")
    Plots.xlabel!(f2, "")
    Plots.ylims!(f2, -π, π)

    return f
end  


function plot_CESTseq(seq::CESTSequence)

    f = CairoMakie.Figure(size=(750,500))
    axL = CairoMakie.Axis(f[1,1], yaxisposition=:left, limits=(nothing, nothing, 0.0, nothing), ylabel="|B₁| [μT]", xlabel="Time [s]", leftspinecolor=:red, ylabelcolor=:red, yticklabelcolor=:red, ytickcolor=:red, title=seq.Name)
    axR = CairoMakie.Axis(f[1,1], yaxisposition=:right, limits=(nothing, nothing, -π, π), ylabel="ϕ [rad]", rightspinecolor=:orange, ylabelcolor=:orange, yticklabelcolor=:orange, ytickcolor=:orange, leftspinecolor=:red)
    CairoMakie.linkxaxes!(axL, axR)

    for n = eachindex(seq.Events)
        if seq.Events[n].EventType === :RF
            CairoMakie.lines!(axL, seq.Events[n].T, abs.(seq.Events[n].RF), label="",color=:red)
            CairoMakie.lines!(axR, seq.Events[n].T, angle.(seq.Events[n].RF),label="", color=:orange, linestyle=:dash)
        elseif seq.Events[n].EventType === :Delay
            CairoMakie.lines!(axL, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", color=:red)
        elseif seq.Events[n].EventType === :Spoiling
            CairoMakie.lines!(axR, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", linestyle=:dot, color=:blue)
        elseif seq.Events[n].EventType === :ADC
            CairoMakie.lines!(axR, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", linestyle=:solid, color=:green)
        end
    end

    return f
end  

@with_kw struct CESTParamsCW
    N::Int64
    R::Array{Float64,2}
    C::Array{Float64,1}
    Ω::Array{Float64,1}
    Kx::Array{Float64,2}
    Ky::Array{Float64,2}
    Kz::Array{Float64,2}
    B1::Float64
    ϕ::Float64
end

@with_kw struct CESTParamsPulsed # Usually the case for modulated amplitude and phase
    N::Int64
    R::Array{Float64,2}
    C::Array{Float64,1}
    Ω::Array{Float64,1}
    Kx::Array{Float64,2}
    Ky::Array{Float64,2}
    Kz::Array{Float64,2}
    B1::ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{LinRange{Float64, Int64}}} 
    ϕ::ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{LinRange{Float64, Int64}}}
end

"""
    function interp_b1(block::SequenceBlock)
    
        Interpolates RF given the Event block.
    
        Returns a function of B1(t) and ϕ(t), where T = block.T
"""
function interp_b1(block::SequenceBlock)
    fb1 = scale(interpolate(abs.(block.RF), BSpline(Constant(Previous))), LinRange(block.T[1], block.T[end], length(block.T)))
    fphi = scale(interpolate(angle.(block.RF) .+ event.ϕ, BSpline(Constant(Previous))), LinRange(block.T[1], block.T[end], length(block.T)))    
    return fb1, fphi
end

"""
    Method 1

    function create_params(model::CESTModel, sys::CESTSys, Δω::Float64)
    Create a NamedTuple with backbone information for simulation.

"""
function create_params(model::CESTModel, sys::CESTSys, Δω::Float64)
    # Δω = B1 in Hz

    N = model.nSpecies
    f = model.fⁿ # fractional concentration vector [au]
    k = model.Cⁿ # exchange rate vector [Hz]

    w = model.wⁿppm .* (sys.γ * sys.B₀)
    Mz0 = 1.0

    (Kx, Ky, Kz) = (zeros(N, N), zeros(N, N), zeros(N, N))
    R = zeros(3N, 3N)
    C = zeros(3N)
    Ω = zeros(N)

    for n in 1:N
        Ω[n] = 2π * (w[n] + w[1] - Δω)
        idx = 3 * (n - 1) 
        R[idx+1, idx+1] = -model.r2ⁿ[n]
        R[idx+2, idx+2] = -model.r2ⁿ[n]
        R[idx+3, idx+3] = -model.r1ⁿ[n]
        C[idx+3] = Mz0 * model.r1ⁿ[n] * f[n]
    end

    # TODO: deal with the below
    # A bit silly as this is inherited from the original script and 
    # isn't in the same order as M0 and K but i don't want to change too much for now
    # I should prob change it if considering more than 100 pools - but thats unlikely (FOR NOW!)

    Kx[1, :] .= [-sum(k .* f[2:end]); k]
    Ky[1, :] .= [-sum(k .* f[2:end]); k]
    Kz[1, :] .= [-sum(k .* f[2:end]); k]

    # other pools
    for n in 2:N
        Kx[n, 1] = k[n-1] * f[n]
        Ky[n, 1] = k[n-1] * f[n]
        Kz[n, 1] = k[n-1] * f[n]

        Kx[n, n] = -k[n-1]
        Ky[n, n] = -k[n-1]
        Kz[n, n] = -k[n-1]
        
    end

    p = (; N=N, R=R, C=C, Ω=Ω, Kx=Kx, Ky=Ky, Kz=Kz)
    return p
end


"""
    Method 2

    function create_params(p::NamedTuple, event::SequenceBlock, type::Symbol)

    Create a parameter Struct for simulation

    Arguments
        - p: Initialisation NamedTuple parameter (see Method 1)
        - event: the Event block with RF and T information
        - type: one of :Relax, :Pulsed, or :CW 

    Returns: One of
        - ::CESTParamsFP, or
        - ::CESTParamsPulsed, or
        - ::CESTParamsCW

"""
function create_params(p::NamedTuple, event::SequenceBlock, type::Symbol)
    if type === :Relax
        p = CESTParamsFP(N=p.N, R=p.R, C=p.C, Ω=p.Ω, Kx=p.Kx, Ky=p.Ky, Kz=p.Kz)
    elseif type === :Pulsed
        fB₁, fϕ = interp_b1(event)
        p = CESTParamsPulsed(N=p.N, R=p.R, C=p.C, Ω=p.Ω, Kx=p.Kx, Ky=p.Ky, Kz=p.Kz, B1=fB₁, ϕ=fϕ)
    elseif type === :CW
        p = CESTParamsCW(N=p.N, R=p.R, C=p.C, Ω=p.Ω, Kx=p.Kx, Ky=p.Ky, Kz=p.Kz, B1=abs(event.RF), ϕ=angle(event.RF) + event.ϕ)
    end
    return p
end



## WIP
module simple
    function makeBlockPulse(A::Float64, 
                            T::Array{Float64, 2},
                            Δt::Float64=1e-6)
        N = floor((T[2] - T[1]) / Δt) |> Int
        T[2] = (T[1]:Δt:T[2])[end] # To deal with suboptimal rounding
        A = convert(ComplexF64, A)
        pulse = (; A = repeat([A], N), T = LinRange(T[1], T[2], N))
        return pulse
    end


    # Pulse with pre-defined A and times
    function makeArbitraryPulse(A, T)
        if size(A) != size(T)
            throw(ArgumentError("Amplitude and Time arrays have different sizes. Please check"))
        end

        if istype(A[1]) != Union{ComplexF64, ComplexF32, ComplexF64}
            throw(ArgumentError("Amplitude must be a Complex Array"))
        end

        Δt = T[2] - T[1]
        pulse = (; A = A, T = T, Δt = Δt)
        return pulse
    end
end

module pulseq
    # Ported from pyPulseq
    using DSP
    using LinearAlgebra
    using Trapz
    function make_sinc_pulse(
        flip_angle::Float64;
        apodization::Float64=0.0,
        delay::Float64=0.0,
        duration::Float64=4e-3,
        dwell::Float64=0.0,
        system,
        center_pos::Float64=0.5,
        freq_offset::Float64=0.0,
        phase_offset::Float64=0.0,
        time_bw_product::Float64=4.0,
        use::String=""
    )
        """
        Creates a radio-frequency sinc pulse event and optionally accompanying slice select and slice select rephasing
        trapezoidal gradient events.

        Parameters
        ----------
        flip_angle : Float64
            Flip angle in radians.
        apodization : Float64, default=0.0
            Apodization.
        center_pos : Float64, default=0.5
            Position of peak.5 (midway).
        delay : Float64, default=0.0
            Delay in seconds (s).
        duration : Float64, default=4e-3
            Duration in seconds (s).
        dwell : Float64, default=0.0
        freq_offset : Float64, default=0.0
            Frequency offset in Hertz (Hz).
        max_grad : Float64, default=0.0
            Maximum gradient strength of accompanying slice select trapezoidal event.
        max_slew : Float64, default=0.0
            Maximum slew rate of accompanying slice select trapezoidal event.
        phase_offset : Float64, default=0.0
            Phase offset in Hertz (Hz).
        return_gz : Bool, default=false
            Boolean flag to indicate if slice-selective gradient has to be returned.
        slice_thickness : Float64, default=0.0
            Slice thickness of accompanying slice select trapezoidal event. The slice thickness determines the area of the
            slice select event.
        system : Opts, default=Opts()
            System limits. Default is a system limits object initialized to default values.
        time_bw_product : Float64, default=4.0
            Time-bandwidth product.
        use : String, default=""
            Use of radio-frequency sinc pulse. Must be one of "excitation", "refocusing" or "inversion".

        Returns
        -------
        rf : NamedTuple
            Radio-frequency sinc pulse event.
        gz : NamedTuple, optional
            Accompanying slice select trapezoidal gradient event. Returned only if `slice_thickness` is provided.
        gzr : NamedTuple, optional
            Accompanying slice select rephasing trapezoidal gradient event. Returned only if `slice_thickness` is provided.

        Raises
        ------
        ArgumentError
            If invalid `use` parameter was passed. Must be one of "excitation", "refocusing" or "inversion".
            If `return_gz=true` and `slice_thickness` was not provided.
        """

        if isempty(system)
            system = defaultsystem()
        end

        valid_pulse_uses = ["excitation", "refocusing", "inversion", "saturation"]
        if !isempty(use) && !(use in valid_pulse_uses)
            throw(ArgumentError("Invalid use parameter. Must be one of $valid_pulse_uses. Passed: $use"))
        end

        if dwell == 0.0
            dwell = system.rf_raster_time
        end

        if duration <= 0.0
            throw(ArgumentError("RF pulse duration must be positive."))
        end

        bandwidth = time_bw_product / duration
        alpha = apodization
        n_samples = round(Int, duration / dwell)
        t = (collect(1:n_samples) .- 0.5) .* dwell
        tt = t .- (duration * center_pos)
        window = 1 .- alpha .+ alpha .* cos.(2 * π * tt / duration)
        signal = window .* sinc.(bandwidth .* tt)
        flip = sum(signal) * dwell * 2 * π
        signal = signal .* (flip_angle / flip)

        rf = (
            type="rf",
            signal=signal,
            t=t,
            shape_dur=n_samples * dwell,
            freq_offset=freq_offset,
            phase_offset=phase_offset,
            dead_time=system.rf_dead_time,
            ringdown_time=system.rf_ringdown_time,
            delay=delay,
            use=use,
            name="Sinc"
        )

        if rf.dead_time > rf.delay
            @warn "Specified RF delay $(rf.delay*1e6) µs is less than the dead time $(rf.dead_time*1e6) µs. Delay was increased to the dead time."
            rf = merge(rf, (delay=rf.dead_time,))
        end

        negative_zero_indices = findall(x -> x == -0.0, rf.signal)
        rf.signal[negative_zero_indices] .= 0.0

        return rf
    end

    function make_gauss_pulse(
        flip_angle::Float64;
        apodization::Float64=0.0,
        delay::Float64=0.0,
        duration::Float64=0.0,
        dwell::Float64=0.0,
        system,
        center_pos::Float64=0.5,
        freq_offset::Float64=0.0,
        phase_offset::Float64=0.0,
        time_bw_product::Float64=3.0,
        use::String=""
    )

        valid_pulse_uses = ["excitation", "refocusing", "inversion", "saturation"]
        if !isempty(use) && !(use in valid_pulse_uses)
            throw(ArgumentError("Invalid use parameter. Must be one of $valid_pulse_uses. Passed: $use"))
        end

        if dwell == 0.0
            dwell = system.rf_raster_time
        end

        if duration <= 0.0
            throw(ArgumentError("RF pulse duration must be positive."))
        end

        bandwidth = time_bw_product / duration
        alpha = apodization
        n_samples = round(Int, duration / dwell)
        t = (collect(1:n_samples) .- 0.5) .* dwell
        tt = t .- (duration * center_pos)
        window = 1 .- alpha .+ alpha .* cos.(2 * π * tt / duration)
        # signal = window .* gauss.(bandwidth .* tt)
        signal = window .* exp.(-π.*(bandwidth .* tt).^2);
        flip = sum(signal) * dwell * 2 * π
        signal = signal .* (flip_angle / flip)

        rf = (
            type="rf",
            signal=signal,
            t=t,
            shape_dur=n_samples * dwell,
            freq_offset=freq_offset,
            phase_offset=phase_offset,
            dead_time=system.rf_dead_time,
            ringdown_time=system.rf_ringdown_time,
            delay=delay,
            use=use,
            name="Gauss"
        )

        if rf.dead_time > rf.delay
            @warn "Specified RF delay $(rf.delay*1e6) µs is less than the dead time $(rf.dead_time*1e6) µs. Delay was increased to the dead time."
            rf = merge(rf, (delay=rf.dead_time,))
        end

        negative_zero_indices = findall(x -> x == -0.0, rf.signal)
        rf.signal[negative_zero_indices] .= 0.0

        return rf

        function gauss(x)
            #gauss Calculate the Gaussian function:
            #   gauss(x) = exp(-pi*x^2)
            
            # This is a useful helper function for those without the signal 
            # processing toolbox 
            
            y = exp(-pi*x.^2);
            return y
        end
    end


    function make_gauss_pulse_hanning(
        flip_angle::Float64;
        apodization::Float64=0.0,
        delay::Float64=0.0,
        duration::Float64=0.0,
        dwell::Float64=0.0,
        system,
        center_pos::Float64=0.5,
        freq_offset::Float64=0.0,
        phase_offset::Float64=0.0,
        time_bw_product::Float64=3.0,
        use::String=""
    )

        valid_pulse_uses = ["excitation", "refocusing", "inversion", "saturation"]
        if !isempty(use) && !(use in valid_pulse_uses)
            throw(ArgumentError("Invalid use parameter. Must be one of $valid_pulse_uses. Passed: $use"))
        end

        if dwell == 0.0
            dwell = system.rf_raster_time
        end

        if duration <= 0.0
            throw(ArgumentError("RF pulse duration must be positive."))
        end

        bandwidth = time_bw_product / duration
        alpha = apodization
        n_samples = round(Int, duration / dwell)
        t = (collect(1:n_samples) .- 0.5) .* dwell
        tt = t .- (duration * center_pos)
        window = 1 .- alpha .+ alpha .* cos.(2 * π * tt / duration)
        
        signal = window .* exp.(-pi.*(bandwidth .* tt).^2);
        # signal = window .* gauss.(bandwidth .* tt)
        flip = sum(signal) * dwell * 2 * π
        signal = signal .* (flip_angle / flip)
        
        # z = diff(t;dims=1)' * (signal[1:end-1,:] + signal[2:end,:])/2;

        hanning_shape = hanning(length(signal));
        signal = hanning_shape ./ trapz(t, hanning_shape) * (flip_angle ./ (2π));


        rf = (
            type="rf",
            signal=signal,
            t=t,
            shape_dur=n_samples * dwell,
            freq_offset=freq_offset,
            phase_offset=phase_offset,
            dead_time=system.rf_dead_time,
            ringdown_time=system.rf_ringdown_time,
            delay=delay,
            use=use,
            name="Hann. Gauss"
        )

        if rf.dead_time > rf.delay
            @warn "Specified RF delay $(rf.delay*1e6) µs is less than the dead time $(rf.dead_time*1e6) µs. Delay was increased to the dead time."
            rf = merge(rf, (delay=rf.dead_time,))
        end

        negative_zero_indices = findall(x -> x == -0.0, rf.signal)
        rf.signal[negative_zero_indices] .= 0.0

        return rf
    end

    function make_rect_pulse(
        flip_angle::Float64;
        delay::Float64=0.0,
        duration::Float64=0.0,
        dwell::Float64=0.0,
        system,
        center_pos::Float64=0.5,
        freq_offset::Float64=0.0,
        phase_offset::Float64=0.0,
        time_bw_product::Float64=3.0,
        use::String="",
    )

        valid_pulse_uses = ["excitation", "refocusing", "inversion", "saturation"]
        if !isempty(use) && !(use in valid_pulse_uses)
            throw(ArgumentError("Invalid use parameter. Must be one of $valid_pulse_uses. Passed: $use"))
        end

        if dwell == 0.0
            dwell = system.rf_raster_time
        end

        if duration <= 0.0
            throw(ArgumentError("RF pulse duration must be positive."))
        end

        bandwidth = 1/(4*duration);
        n_samples = round(Int, duration / dwell)
        t = [0; n_samples] .* dwell
        
        signal = flip_angle/(2π)/duration*ones(size(t));


        rf = (
            type="rf",
            signal=signal,
            t=t,
            shape_dur=n_samples * dwell,
            freq_offset=freq_offset,
            phase_offset=phase_offset,
            dead_time=system.rf_dead_time,
            ringdown_time=system.rf_ringdown_time,
            delay=delay,
            use=use,
            name="Rect"
        )

        if rf.dead_time > rf.delay
            @warn "Specified RF delay $(rf.delay*1e6) µs is less than the dead time $(rf.dead_time*1e6) µs. Delay was increased to the dead time."
            rf = merge(rf, (delay=rf.dead_time,))
        end

        negative_zero_indices = findall(x -> x == -0.0, rf.signal)
        rf.signal[negative_zero_indices] .= 0.0

        return rf

    end
end


