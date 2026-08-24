
include("Aux.jl")
import Plots
import Makie

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
    # 3 types of EVentName accepted: RF, Delay, ADC, and Spoiling
    RF::Union{Array{ComplexF64,1}, ComplexF64, Nothing}=nothing
    T::Union{Array{Float64,1}, Nothing}=[0.0]
    Δω::Union{Float64, Nothing}=0.0
    ϕ::Union{Float64, Nothing}=0.0
    EventType::Symbol # removed default being :Ignore to reduce human errors
    EventName::String=""
end

@with_kw mutable struct CESTSequence
    System::Union{NamedTuple, CESTSys}=defaultSystem()
    Events::Array{SequenceBlock,1}=[SequenceBlock(EventName="Seq Init", EventType=:Ignore)]
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
# TODO: Replace this with exported pre-defined `SequenceBlock`s instead?
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
                # rf = convert(Array{ComplexF64,1}, repeat([seq.Events[n].RF], length(seq.Events[n].T)))
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

    f = Makie.Figure(size=(750,500))
    axL = Makie.Axis(f[1,1], yaxisposition=:left, limits=(nothing, nothing, 0.0, nothing), ylabel="|B₁| [μT]", xlabel="Time [s]", leftspinecolor=:red, ylabelcolor=:red, yticklabelcolor=:red, ytickcolor=:red, title=seq.Name)
    axR = Makie.Axis(f[1,1], yaxisposition=:right, limits=(nothing, nothing, -π, π), ylabel="ϕ [rad]", rightspinecolor=:orange, ylabelcolor=:orange, yticklabelcolor=:orange, ytickcolor=:orange, leftspinecolor=:red)
    Makie.linkxaxes!(axL, axR)

    for n = eachindex(seq.Events)
        if seq.Events[n].EventType === :RF
            mag =  seq.Events[n].RF
            if isa(mag, ComplexF64)
                mag = mag .+ zeros(length(seq.Events[n].T))
            end
            Makie.lines!(axL, seq.Events[n].T, abs.(mag), label="",color=:red)
            Makie.lines!(axR, seq.Events[n].T, angle.(mag),label="", color=:orange, linestyle=:dash)
        elseif seq.Events[n].EventType === :Delay
            Makie.lines!(axL, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", color=:red)
        elseif seq.Events[n].EventType === :Spoiling
            Makie.lines!(axR, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", linestyle=:dot, color=:blue)
        elseif seq.Events[n].EventType === :ADC
            Makie.lines!(axR, seq.Events[n].T, zeros(size(seq.Events[n].T)), label="", linestyle=:solid, color=:green)
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

## Some default blocks
blockADC = SequenceBlock(T=0.0:1e-6:2e-6, EventType=:ADC, EventName="Default readout, 2 μs")
blockdummyRF = SequenceBlock(RF=1 + 0.0im, T=0.0:1e-6:1e-3, EventType=:RF, EventName="Dummy RF, 1 ms")
dummySeq = CESTSequence()
addBlock!(dummySeq, blockdummyRF)
addBlock!(dummySeq, blockADC)

##

# @with_kw struct CESTParamsPulsed # Usually the case for modulated amplitude and phase
#     N::Int64
#     R::Array{Float64,2}
#     C::Array{Float64,1}
#     Ω::Array{Float64,1}
#     Kx::Array{Float64,2}
#     Ky::Array{Float64,2}
#     Kz::Array{Float64,2}
#     B1::ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{LinRange{Float64, Int64}}} 
#     ϕ::ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Constant{Previous, Throw{OnGrid}}}, Tuple{LinRange{Float64, Int64}}}
# end

@with_kw struct CESTParamsPulsed # Usually the case for modulated amplitude and phase
    N::Int64
    R::Array{Float64,2}
    C::Array{Float64,1}
    Ω::Array{Float64,1}
    Kx::Array{Float64,2}
    Ky::Array{Float64,2}
    Kz::Array{Float64,2}
    B1::FastInterpolations.AbstractInterpolant{Float64, Float64}
    ϕ::FastInterpolations.AbstractInterpolant{Float64, Float64}
end

"""
    function interp_b1(block::SequenceBlock)
    
        Interpolates RF given the Event block.
    
        Returns a function of B1(t) and ϕ(t), where T = block.T
"""
# function interp_b1(block::SequenceBlock, ver::String) # Deprecated
#     fB₁ = scale(interpolate(abs.(block.RF), BSpline(Constant(Previous))), LinRange(block.T[1], block.T[end], length(block.T)))
#     fϕ = scale(interpolate(angle.(block.RF) .+ block.ϕ, BSpline(Constant(Previous))), LinRange(block.T[1], block.T[end], length(block.T)))    
#     return fB₁, fϕ
# end

function interp_b1(event::SequenceBlock, order::Symbol=:Linear)
    if order == :Linear
        fB₁ = linear_interp(event.T, abs.(event.RF)             ; extrap=ClampExtrap())
        fϕ  = linear_interp(event.T, angle.(event.RF) .+ event.ϕ; extrap=ClampExtrap())
    elseif order == :Constant
        fB₁ = constant_interp(event.T, abs.(event.RF)             ; side=LeftSide(), extrap=ClampExtrap())
        fϕ  = constant_interp(event.T, angle.(event.RF) .+ event.ϕ; side=LeftSide(), extrap=ClampExtrap())
    elseif order == :Quadratic 
        fB₁ = quadratic_interp(event.T, abs.(event.RF)             ; extrap=ClampExtrap())
        fϕ  = quadratic_interp(event.T, angle.(event.RF) .+ event.ϕ; extrap=ClampExtrap())
    elseif order == :Cubic 
        fB₁ = cubic_interp(event.T, abs.(event.RF)             ; extrap=ClampExtrap())
        fϕ  = cubic_interp(event.T, angle.(event.RF) .+ event.ϕ; extrap=ClampExtrap())
    end
    return fB₁, fϕ
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
    δω = model.δω
    Mz0 = 1.0

    (Kx, Ky, Kz) = (zeros(N, N), zeros(N, N), zeros(N, N))
    R = zeros(3N, 3N)
    C = zeros(3N)
    Ω = zeros(N)

    for n = 1:N
        # Ω[n] = 2π * (w[n] + w[1] - Δω) # moved away from this
        Ω[n] = 2π * (w[n] - Δω + δω) # Made this consistent across all species;
        # Previously, we only needed to change the water offset to simulate B0 inhomo, and all other species offsets are adjusted accordingly. 
        # This is practically very convenient but not very proper.
        idx = 3 * (n - 1) 
        R[idx+1, idx+1] = -model.r2ⁿ[n]
        R[idx+2, idx+2] = -model.r2ⁿ[n]
        R[idx+3, idx+3] = -model.r1ⁿ[n]
        C[idx+3] = Mz0 * model.r1ⁿ[n] * f[n]
    end
    # Ω[1] = 2π * (w[1] - Δω) # moved away from this

    # TODO: deal with the below
    # A bit silly as this is inherited from the original script and 
    # isn't in the same order as M0 and K but i don't want to change too much for now
    # I should prob change it if considering more than 100 pools - but thats unlikely (FOR NOW!)


    # other pools
    if N > 1
        Kx[1, :] .= [-sum(k[2:end] .* f[2:end]); k[2:end]]
        Ky[1, :] .= [-sum(k[2:end] .* f[2:end]); k[2:end]]
        Kz[1, :] .= [-sum(k[2:end] .* f[2:end]); k[2:end]]
        for n = 2:N
            Kx[n, 1] = k[n] * f[n]
            Ky[n, 1] = k[n] * f[n]
            Kz[n, 1] = k[n] * f[n]

            Kx[n, n] = -k[n]
            Ky[n, n] = -k[n]
            Kz[n, n] = -k[n]
            # Kx[n, 1] = k[n-1] * f[n]
            # Ky[n, 1] = k[n-1] * f[n]
            # Kz[n, 1] = k[n-1] * f[n]
# 
            # Kx[n, n] = -k[n-1]
            # Ky[n, n] = -k[n-1]
            # Kz[n, n] = -k[n-1]
            
        end
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
function create_params(p::NamedTuple, event::SequenceBlock, type::Symbol; order::Symbol=:Linear)
    if type === :Relax
        p = CESTParamsFP(N=p.N, R=p.R, C=p.C, Ω=p.Ω, Kx=p.Kx, Ky=p.Ky, Kz=p.Kz)
    elseif type === :Pulsed
        fB₁, fϕ = interp_b1(event, order)
        p = CESTParamsPulsed(N=p.N, R=p.R, C=p.C, Ω=p.Ω, Kx=p.Kx, Ky=p.Ky, Kz=p.Kz, B1=fB₁, ϕ=fϕ)
    elseif type === :CW
        p = CESTParamsCW(N=p.N, R=p.R, C=p.C, Ω=p.Ω, Kx=p.Kx, Ky=p.Ky, Kz=p.Kz, B1=abs(event.RF), ϕ=angle(event.RF) + event.ϕ)
    end
    return p
end


function plot_CESTtraj(RF::SequenceBlock, ADC::SequenceBlock, Δω::Float64, sys::CESTSys, model::CESTModel; species::Int=1)
    # If no adiabatic
    Makie.set_theme!()
    f = Makie.Figure(size=(1500,1000));
    ax = [Makie.Axis(f[1,n]) for n = 1:3]
    clrs = [:firebrick, :blue]
    ws = [-Δω Δω] .* sys.B₀MHz
    rs = [0.0 0.0]
    
    if isa(RF.RF, ComplexF64)
        # pRF = repeat([RF.RF], length(RF.T))
        pRF = RF.RF .+ zeros(ComplexF64, length(RF.T))
    else
        pRF = RF.RF
    end
    Makie.lines!(ax[1], RF.T, pRF .|> abs)
    Makie.lines!(ax[2], RF.T, pRF .|> angle)
    ax[1].xlabel="Time [s]"; ax[1].ylabel="B₁ [μT]"; ax[1].title = "|RF|"
    ax[2].xlabel="Time [s]"; ax[2].ylabel="ϕ [rad]"; ax[2].title = "∠RF"
    RF2 = deepcopy(RF)
    axs = [Makie.Axis(f[2,n], aspect=1) for n = 1:3]
    ax3 = Makie.Axis3(f[:,4:5], title="Trajectory", aspect=:equal, limits=(-1,1,-1,1,-1,1), azimuth=1.6π)
    for nw = eachindex(ws)
        Dw     = ws[nw]
        RF2.Δω = Dw
        seq    = CESTSequence(System=sys)
        addBlock!(seq, RF2)
        addBlock!(seq, ADC)
        simp = CESTSimParams(Δω = [Δω], model = model, solver="Rodas5P(autodiff=AutoFiniteDiff())", sys=sys, reltol=1e-6, abstol=1e-6)
        r = simulate_CEST(seq, simp, true, true; species=species)
        Mx = r.Mx; My = r.My; Mz = r.Mz
        Makie.lines!(ax[3], r.t, Mx, label="Mx")
        Makie.lines!(ax[3], r.t, My, label="My")
        Makie.lines!(ax[3], r.t, Mz, label="Mz")
        ax[3].limits=(nothing,nothing, 0.0,1)


        Makie.lines!(axs[1], Mx, My, color=clrs[nw]); Makie.arc!(axs[1], 0.0, 1, -π, π, linestyle=:dash, color=:black, alpha=.5); axs[1].xlabel="X"; axs[1].ylabel="Y"
        Makie.lines!(axs[2], Mx, Mz, color=clrs[nw]); Makie.arc!(axs[2], 0.0, 1, -π, π, linestyle=:dash, color=:black, alpha=.5); axs[2].xlabel="X"; axs[2].ylabel="Z"
        Makie.lines!(axs[3], My, Mz, color=clrs[nw]); Makie.arc!(axs[3], 0.0, 1, -π, π, linestyle=:dash, color=:black, alpha=.5); axs[3].xlabel="Y"; axs[3].ylabel="Z"

        Makie.arc!(ax3, 0.0, 1, -π, π, linestyle=:dash, color=:black, alpha=.5)
        Makie.lines!(ax3, zeros(63), sin.(-π:.1:π), cos.(-π:.1:π), linestyle=:dash, color=:black, alpha=.5)
        Makie.lines!(ax3, sin.(-π:.1:π), zeros(63), cos.(-π:.1:π), linestyle=:dash, color=:black, alpha=.5)
        Makie.lines!(ax3, [0 0], [0 0], [-1.2 1.2], color=:black, alpha=.6, linecap=:square)
        Makie.lines!(ax3, [0 0], [-1.2 1.2], [0 0], color=:black, alpha=.6, linecap=:square)
        Makie.lines!(ax3, [-1.2 1.2], [0 0], [0 0], color=:black, alpha=.6, linecap=:square)
        Makie.lines!(ax3, Mx, My, Mz,color=clrs[nw])
        Makie.arrows3d!(ax3, [0.0], [0.0], [0.0], [Mx[end]], [My[end]], [Mz[end]], shaftradius=0.01, tipradius=0.04, tiplength=0.1, color=clrs[nw])
        Makie.text!(ax3, Mx[end], My[end], Mz[end]; text="  (" * string(round(Mx[end];sigdigits=2)) * "," * string(round(My[end];sigdigits=2)) * "," * string(round(Mz[end];sigdigits=2)) * ")")

        println(string("Mx = ", Mx[end], " My = ", My[end]," Mz = ", Mz[end]))
        # Makie.Label(f[4+nw,:]; text=string("Mx = ", Mx[end], " My = ", My[end]," Mz = ", Mz[end]), fontsize=15, color=clrs[nw])
        rs[nw] = Mz[end]
    end
    Makie.Label(f[0,:]; text="Δω = ±"*string(Δω) * " ppm, MTRAsym = " * string(round(((rs[1] - rs[2])/rs[1])*100; sigdigits=4)) * " %", fontsize=20)
    # Makie.Label(f[4,:]; text=, fontsize=15)
    println(((rs[1] - rs[2])/rs[1])*100)
    return f
end

function plot_CESTtraj(seq::CESTSequence, Δω::Float64, sys::CESTSys, model::CESTModel; both::Bool=true, species::Int=1)

    scale = model.fⁿ[species]
    # If no adiabatic
    Makie.set_theme!()
    f = Makie.Figure(size=(1500,1000));
    ax = [Makie.Axis(f[1,n]) for n = 1:3]
    clrs = [:firebrick, :blue]
    if both
        ws = [-Δω Δω] .* sys.B₀MHz
        rs = [0.0 0.0]
    else
        ws = [Δω] .* sys.B₀MHz
        rs = [0.0]
    end
    
    for n = eachindex(seq.Events)
        if seq.Events[n].EventType == :RF
            if isa(seq.Events[n].RF, ComplexF64)
                RF = seq.Events[n].RF .+ zeros(ComplexF64, length(seq.Events[n].T))
            else
                RF = seq.Events[n].RF
            end
            Makie.lines!(ax[1], seq.Events[n].T, RF .|> abs)
            Makie.lines!(ax[2], seq.Events[n].T, RF .|> angle)
        end
    end
    ax[1].xlabel="Time [s]"; ax[1].ylabel="B₁ [μT]"; ax[1].title = "|RF|"
    ax[2].xlabel="Time [s]"; ax[2].ylabel="ϕ [rad]"; ax[2].title = "∠RF"
    # RF2 = deepcopy(RF)
    seq2 = deepcopy(seq)
    axs = [Makie.Axis(f[2,n], aspect=1) for n = 1:3]
    ax3 = Makie.Axis3(f[:,4:5], title="Trajectory", aspect=:equal, limits=scale.*(-1,1,-1,1,-1,1), azimuth=1.6π)
    for nw = eachindex(ws)
        Dw = ws[nw]
        for n = eachindex(seq.Events)
            seq2.Events[n].Δω = Dw
        end
        simp = CESTSimParams(Δω = [Δω], model = model, solver="Rodas5P(autodiff=AutoFiniteDiff())", sys=sys, reltol=1e-6, abstol=1e-6)
        r = simulate_CEST(seq2, simp, true, true; species=species)
        Mx = r.Mx; My = r.My; Mz = r.Mz
        Makie.lines!(ax[3], r.t, Mx, label="Mx")
        Makie.lines!(ax[3], r.t, My, label="My")
        Makie.lines!(ax[3], r.t, Mz, label="Mz")
        ax[3].limits=(nothing,nothing, 0.0,scale*1)


        Makie.lines!(axs[1], Mx, My, color=clrs[nw]); Makie.arc!(axs[1], 0.0, scale*1, -scale*π, scale*π, linestyle=:dash, color=:black, alpha=.5); axs[1].xlabel="X"; axs[1].ylabel="Y"
        Makie.lines!(axs[2], Mx, Mz, color=clrs[nw]); Makie.arc!(axs[2], 0.0, scale*1, -scale*π, scale*π, linestyle=:dash, color=:black, alpha=.5); axs[2].xlabel="X"; axs[2].ylabel="Z"
        Makie.lines!(axs[3], My, Mz, color=clrs[nw]); Makie.arc!(axs[3], 0.0, scale*1, -scale*π, scale*π, linestyle=:dash, color=:black, alpha=.5); axs[3].xlabel="Y"; axs[3].ylabel="Z"


        Makie.arc!(ax3, 0.0, 1*scale, -π, π, linestyle=:dash, color=:black, alpha=.5)
        Makie.lines!(ax3, zeros(63), scale.*sin.(-π:.1:π), scale.*cos.(-π:.1:π), linestyle=:dash, color=:black, alpha=.5)
        Makie.lines!(ax3, scale.*sin.(-π:.1:π), zeros(63), scale.*cos.(-π:.1:π), linestyle=:dash, color=:black, alpha=.5)
        Makie.lines!(ax3, [0 0], [0 0], scale.*[-1.2 1.2], color=:black, alpha=.6, linecap=:square)
        Makie.lines!(ax3, [0 0], scale.*[-1.2 1.2], [0 0], color=:black, alpha=.6, linecap=:square)
        Makie.lines!(ax3, scale.*[-1.2 1.2], [0 0], [0 0], color=:black, alpha=.6, linecap=:square)
        Makie.lines!(ax3, Mx, My, Mz,color=clrs[nw])
        Makie.arrows3d!(ax3, [0.0], [0.0], [0.0], [Mx[end]], [My[end]], [Mz[end]], shaftradius=0.01, tipradius=0.04, tiplength=0.1, color=clrs[nw])
        Makie.text!(ax3, Mx[end], My[end], Mz[end]; text="  (" * string(round(Mx[end];sigdigits=2)) * "," * string(round(My[end];sigdigits=2)) * "," * string(round(Mz[end];sigdigits=2)) * ")")

        println(string("Mx = ", Mx[end], " My = ", My[end]," Mz = ", Mz[end]))
        # Makie.Label(f[4+nw,:]; text=string("Mx = ", Mx[end], " My = ", My[end]," Mz = ", Mz[end]), fontsize=15, color=clrs[nw])
        rs[nw] = Mz[end]
        println(string("θ = ", atan(Mx[end]/Mz[end])))
    end
    if both
        Makie.Label(f[0,:]; text="Δω = ±"*string(Δω) * " ppm, MTRAsym = " * string(round(((rs[1] - rs[2])/rs[1])*100; sigdigits=4)) * " %", fontsize=20)
        println(((rs[1] - rs[2])/rs[1])*100)
    else
        Makie.Label(f[0,:]; text="Δω = ±"*string(Δω) * " ppm", fontsize=20)
    end
        # Makie.Label(f[4,:]; text=, fontsize=15)
    return f
end

function plot_CESTtraj(seq::CESTSequence, sys::CESTSys, model::CESTModel; species::Int=1)

    scale = model.fⁿ[species]
    # If no adiabatic
    Makie.set_theme!()
    f = Makie.Figure(size=(1500,1000));
    ax = [Makie.Axis(f[1,n]) for n = 1:3]
    clrs = [:firebrick, :blue]
    
    for n = eachindex(seq.Events)
        if seq.Events[n].EventType == :RF
            if isa(seq.Events[n].RF, ComplexF64)
                RF = seq.Events[n].RF .+ zeros(ComplexF64, length(seq.Events[n].T))
            else
                RF = seq.Events[n].RF
            end
            Makie.lines!(ax[1], seq.Events[n].T, RF .|> abs)
            Makie.lines!(ax[2], seq.Events[n].T, RF .|> angle)
        end
    end
    ax[1].xlabel="Time [s]"; ax[1].ylabel="B₁ [μT]"; ax[1].title = "|RF|"
    ax[2].xlabel="Time [s]"; ax[2].ylabel="ϕ [rad]"; ax[2].title = "∠RF"
    # RF2 = deepcopy(RF)
    seq2 = deepcopy(seq)
    axs = [Makie.Axis(f[2,n], aspect=1) for n = 1:3]
    ax3 = Makie.Axis3(f[:,4:5], title="Trajectory", aspect=:equal, limits=scale.*(-1,1,-1,1,-1,1), azimuth=1.6π)
    nw = 1

    simp = CESTSimParams(Δω = [0.0], model = model, solver="Rodas5P(autodiff=AutoFiniteDiff())", sys=sys, reltol=1e-6, abstol=1e-6)
    r = simulate_CEST(seq2, simp, true, true; species=species)
    Mx = r.Mx; My = r.My; Mz = r.Mz
    Makie.lines!(ax[3], r.t, Mx, label="Mx")
    Makie.lines!(ax[3], r.t, My, label="My")
    Makie.lines!(ax[3], r.t, Mz, label="Mz")
    ax[3].limits=(nothing,nothing, 0.0,scale*1)


    Makie.lines!(axs[1], Mx, My, color=clrs[nw]); Makie.arc!(axs[1], 0.0, scale*1, -scale*π, scale*π, linestyle=:dash, color=:black, alpha=.5); axs[1].xlabel="X"; axs[1].ylabel="Y"
    Makie.lines!(axs[2], Mx, Mz, color=clrs[nw]); Makie.arc!(axs[2], 0.0, scale*1, -scale*π, scale*π, linestyle=:dash, color=:black, alpha=.5); axs[2].xlabel="X"; axs[2].ylabel="Z"
    Makie.lines!(axs[3], My, Mz, color=clrs[nw]); Makie.arc!(axs[3], 0.0, scale*1, -scale*π, scale*π, linestyle=:dash, color=:black, alpha=.5); axs[3].xlabel="Y"; axs[3].ylabel="Z"


    Makie.arc!(ax3, 0.0, 1*scale, -π, π, linestyle=:dash, color=:black, alpha=.5)
    Makie.lines!(ax3, zeros(63), scale.*sin.(-π:.1:π), scale.*cos.(-π:.1:π), linestyle=:dash, color=:black, alpha=.5)
    Makie.lines!(ax3, scale.*sin.(-π:.1:π), zeros(63), scale.*cos.(-π:.1:π), linestyle=:dash, color=:black, alpha=.5)
    Makie.lines!(ax3, [0 0], [0 0], scale.*[-1.2 1.2], color=:black, alpha=.6, linecap=:square)
    Makie.lines!(ax3, [0 0], scale.*[-1.2 1.2], [0 0], color=:black, alpha=.6, linecap=:square)
    Makie.lines!(ax3, scale.*[-1.2 1.2], [0 0], [0 0], color=:black, alpha=.6, linecap=:square)
    Makie.lines!(ax3, Mx, My, Mz,color=clrs[nw])
    Makie.arrows3d!(ax3, [0.0], [0.0], [0.0], [Mx[end]], [My[end]], [Mz[end]], shaftradius=0.01, tipradius=0.04, tiplength=0.1, color=clrs[nw])
    Makie.text!(ax3, Mx[end], My[end], Mz[end]; text="  (" * string(round(Mx[end];sigdigits=2)) * "," * string(round(My[end];sigdigits=2)) * "," * string(round(Mz[end];sigdigits=2)) * ")")

    println(string("Mx = ", Mx[end], " My = ", My[end]," Mz = ", Mz[end]))

    return f
end

