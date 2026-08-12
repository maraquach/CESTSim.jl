using Parameters
@kwdef mutable struct CESTSys
    B₀::Float64
    # maxSlew::Float64
    # maxGrad::Float64
    # rfRasterTime::Float64
    # rfDeadTime::Float64
    # rfRingDownTime::Float64
    γ::Float64
    γ̄::Float64
    B₀MHz::Float64
end


function defaultSystem(verbose::Bool=true) 
    sys = CESTSys(B₀ = 6.98,
    # maxSlew = 200.0,
    # maxGrad = 70.0,
    # rfRasterTime = 1.67e-3,
    # rfDeadTime = 1e-5,
    # rfRingDownTime = 1e-5,
    γ = γ * 1e-6,
    γ̄ = γ̄ * 1e-6,
    B₀MHz = 6.98 * γ * 1e-6)
    if verbose
        println("Default System Limits for Siemens 7T Magnetom Plus")
    end
    return sys
end

@with_kw mutable struct CESTSimParams
    solver   ::String="Rodas5P(autodiff=AutoFiniteDiff())"
    reltol   ::Float64=1e-4
    abstol   ::Float64=1e-4
    interp   ::Symbol=:Linear
    sys      ::CESTSys=defaultSystem()
    Δω       ::Array{Float64,1}
    model    ::CESTModel
    plotZspec::Bool=true
    opts     ::String=""
end

@with_kw struct CESTParamsFP
    N ::Int64
    R ::Array{Float64,2}
    C ::Array{Float64,1}
    Ω ::Array{Float64,1}
    Kx::Array{Float64,2}
    Ky::Array{Float64,2}
    Kz::Array{Float64,2}
end


"""
    function create_M0(model::CESTModel, M0ᴬ::Float64=1.0)
    
    Returns:
    - M0::Array{Float64, 1}, containing [MxiA, MxiB, MyiA, MyiB, MziA, MziB].
        For 2 pools, this would be [0,0,0,0,MziA,MziB]
        
"""
function create_M0(model::CESTModel, M0ᴬ::Float64=1.0)
    # M0a: Initial Z magnetisation of pool A/water (normalised) 
    N = model.nSpecies
    M0 = zeros(3N)
    M0[2N+1:end] = M0ᴬ .* model.fⁿ .|> Float64
    return M0
end


## Pretty things -- GUI in the making
import ColorSchemes
import GLMakie

const CESTSimColours = get(ColorSchemes.rainbow, range(0.0,1.0, length=10)) |> ColorSchemes.ColorScheme

CESTSimTheme = GLMakie.Theme(
    backgroundcolor = :black,
    size=(700,700),
    textcolor = :white,
    linecolor = :white,
    colormap = :lightrainbow,
    palette = (
        color = CESTSimColours,
        marker = [:circle, :cross, :rect, :xcross, :pentagon],
        linestyle = [:solid, :dash, :dot, :dashdot, :dashdotdot],
        side = [:left, :right]
    ),
    Axis = (
        backgroundcolor = :transparent,
        bottomspinecolor = :white,
        topspinecolor = :white,
        leftspinecolor = :white,
        rightspinecolor = :white,
        xgridcolor = :gray,
        ygridcolor = :gray,
        xtickcolor = :white,
        ytickcolor = :white,
    ),
    Legend = (
        framecolor = :white,
        backgroundcolor = :black,
        fontsize = 20
    ),
    Axis3 = (
        xgridcolor = :gray,
        ygridcolor = :gray,
        zgridcolor = :gray,
        xspinecolor_1 = :white,
        yspinecolor_1 = :white,
        zspinecolor_1 = :white,
        xspinecolor_2 = :white,
        yspinecolor_2 = :white,
        zspinecolor_2 = :white,
        xspinecolor_3 = :white,
        yspinecolor_3 = :white,
        zspinecolor_3 = :white,
        xspinecolor_4 = :white,
        yspinecolor_4 = :white,
        zspinecolor_4 = :white,
        xticklabelpad = 3,
        yticklabelpad = 3,
        zticklabelpad = 6,
        xtickcolor = :white,
        ytickcolor = :white,
        ztickcolor = :white,
    ),
    Colorbar = (
        tickcolor = :white,
        spinecolor = :white,
        topspinecolor = :white,
        bottomspinecolor = :white,
        leftspinecolor = :white,
        rightspinecolor = :white,
    ),
    tooltip = (
        backgroundcolor=:black
    ),
    DataInspector = (
        backgroundcolor=:black,
        textcolor=:white
    )
)


function plot_CESTresult(r, params::CESTSimParams, hasref::Bool=true)
    f = GLMakie.Figure(size=(800,400))
    axZ = GLMakie.Axis(f[1,1], title="Z-spectrum", ylabel="Z(Δω)/Z(M0) [a.u]", xlabel="Δω [Hz]")
    if hasref
        GLMakie.lines!(axZ, params.Δω[2:end], r[2:end]./r[1])
    else
        GLMakie.lines!(axZ, params.Δω, r)
    end
    
    axM = GLMakie.Axis(f[1,2], title="MTRAsym", ylabel="MTRAsym [%]", xlabel="Δω [Hz]")
    
    
    if hasref
        r = r[2:end]./r[1]
        w = params.Δω[2:end]
        w0idx = findfirst(w .== 0.0)
        GLMakie.lines!(axM, w[w0idx+1:end], 100 .*(reverse(r[1:w0idx-1]) .- r[w0idx+1:end]) ./ reverse(r[1:w0idx-1]))
    else
        w = params.Δω
        w0idx = findfirst(w .== 0.0)
        GLMakie.lines!(axM, w[w0idx+1:end], 100 .*(reverse(r[1:w0idx-1]) .- r[w0idx+1:end]) ./ reverse(r[1:w0idx-1]))
    end
    return f, axZ, axM

end
## TODO: Move this to CESTTools.jl
function calcMTR(dw, Z, type::Symbol=:ratio)
        zloc = findfirst(dw .== 0.0)
        neg = reverse(Z[1:zloc-1])
        pos = Z[zloc+1:end]
        if type == :ratio
            mtr = (; vals = 100 .* (neg .- pos) ./ neg, dw = dw[zloc+1:end])
        elseif type == :S0
            mtr = (; vals = 100 .* (neg .- pos), dw = dw[zloc+1:end])
        end
    return mtr

end