using Parameters
@kwdef mutable struct CESTSys
    B₀::Float64
    maxSlew::Float64
    maxGrad::Float64
    rf_raster_time::Float64
    rf_dead_time::Float64
    rf_ringdown_time::Float64
    γ::Float64
    γ̄::Float64
    B₀MHz::Float64
end


function defaultSystem() # Could also be a NamedTuple? 
    sys = CESTSys(B₀ = 6.98,
    maxSlew = 200.0,
    maxGrad = 70.0,
    rf_raster_time = 1.67e-3,
    rf_dead_time = 1e-5,
    rf_ringdown_time = 1e-5,
    γ = γ * 1e-6,
    γ̄ = γ̄ * 1e-6,
    B₀MHz = 6.98 * γ * 1e-6)
    println("Default System Limits for Siemens 7T Magnetom Plus")
    return sys
end

@with_kw mutable struct CESTSimParams
    solver::String="Rodas5P(), reltol = 0.00001, abstol = 0.00001"
    sys::CESTSys=defaultSystem()
    Δω::Array{Float64,1}
    model::CESTModel
    plotZspec::Bool=true
    opts::String=""
end

@with_kw struct CESTParamsFP
    N::Int64
    R::Array{Float64,2}
    C::Array{Float64,1}
    Ω::Array{Float64,1}
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

function plot_CESTresult(r, params::CESTSimParams, hasref::Bool=true)
    f = CairoMakie.Figure(size=(800,400))
    axZ = CairoMakie.Axis(f[1,1], title="Z-spectrum", ylabel="Z(Δω)/Z(M0) [a.u]", xlabel="Δω [Hz]")
    if hasref
        CairoMakie.lines!(axZ, params.Δω[2:end], r[2:end]./r[1])
    else
        CairoMakie.lines!(axZ, params.Δω, r)
    end
    
    axM = CairoMakie.Axis(f[1,2], title="MTRAsym", ylabel="MTRAsym [%]", xlabel="Δω [Hz]")
    
    
    if hasref
        r = r[2:end]./r[1]
        w = params.Δω[2:end]
        w0idx = findfirst(w .== 0.0)
        CairoMakie.lines!(axM, w[w0idx+1:end], 100 .*(reverse(r[1:w0idx-1]) .- r[w0idx+1:end]) ./ reverse(r[1:w0idx-1]))
    else
        w = params.Δω
        w0idx = findfirst(w .== 0.0)
        CairoMakie.lines!(axM, w[w0idx+1:end], 100 .*(reverse(r[1:w0idx-1]) .- r[w0idx+1:end]) ./ reverse(r[1:w0idx-1]))
    end
    return f, axZ, axM

end