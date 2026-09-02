"""Simulation support"""

# include("Aux.jl")
# using DifferentialEquations
# using FastInterpolations

## TODO: Consider if I should remove the first CESTSequence `SequenceBlock` (currently a dummy init block)
# TODO: Alert when no ADC `SequenceBlock`s 
function simulate_CEST(seq::CESTSequence, params::CESTSimParams, resetM0::Bool=false; species::Int=1)
    @assert sum(getfield.(seq.Events, :EventType) .== :ADC) > 0 "💥 No ADC event found in CESTSequence. Aborting."
    # @assert params.model.wⁿppm .* sys.B₀MHz == params.model.wⁿHz "💥 CESTSys and CESTModel not consistent. Check CESTSys.B₀MHz or γ"
    Zspec  = similar(params.Δω)
    XYspec = similar(params.Δω)
    θ      = similar(params.Δω)
    M0_0   = create_M0(params.model)
    p0     = create_params(params.model, params.sys, 0.0)
    nadc   = 0
    for n = eachindex(seq.Events)
        global M0
        if n == 1
            continue
        elseif n == 2
            M0 = M0_0
        end
        # event = deepcopy(seq.Events[n])
        if seq.Events[n].EventType === :RF      
            p1 = create_params(params.model, params.sys, seq.Events[n].Δω)
            if isa(seq.Events[n].RF, ComplexF64)
                p = create_params(p1, seq.Events[n], :CW) # This will produce p::CESTParamsCW
            elseif isa(seq.Events[n].RF, Array{ComplexF64,1})
                p = create_params(p1, seq.Events[n], :Pulsed; order=params.interp) # This will produce p::CESTParamsPulsed
            else
                println("uhOh")
            end
            # println(typeof(p))
            prob = ODEProblem(solve_BME!, M0, (seq.Events[n].T[1], seq.Events[n].T[end]), p) 
            sol = solve(prob, eval(Meta.parse(params.solver)), reltol=params.reltol, abstol=params.abstol, dense=false, save_everystep=false)
            M0 = sol.u[end];
        elseif seq.Events[n].EventType === :Delay || seq.Events[n].EventType === :ADC || seq.Events[n].EventType === :Spoiling
            p = create_params(p0, seq.Events[n], :Relax) # This will produce p::CESTParamsFP
            prob = ODEProblem(solve_BME!, M0, (seq.Events[n].T[1], seq.Events[n].T[end]), p)
            sol = solve(prob, eval(Meta.parse(params.solver)), reltol=1e-6, abstol=1e-6, dense=false, save_everystep=false)
            M0 = sol.u[end];
            if seq.Events[n].EventType === :Spoiling
                M0[params.model.nSpecies * 0 + species] = 0.0 
                M0[params.model.nSpecies * 1 + species] = 0.0 
            elseif seq.Events[n].EventType === :ADC && resetM0
                M0 = M0_0
            end
        else # if :Ignore
            continue
        end
        if seq.Events[n].EventType === :ADC
            nadc += 1 
            Zspec[nadc]  = getindex.(sol.u, params.model.nSpecies * 2 + species)[end]
            Mx           = getindex.(sol.u, species)[end]
            My           = getindex.(sol.u, params.model.nSpecies * 1 + species)[end]
            XYspec[nadc] = sqrt.(Mx.^2 .+ My.^2)
            θ[nadc]      = atan(Mx/Zspec[nadc])
        end
    end
    # return Zspec, XYspec, θ
    return (; Z = Zspec, XY = XYspec, θ = θ, Δω = params.Δω)
end


function simulate_CEST(seq::CESTSequence, params::CESTSimParams, detailed::Bool, resetM0::Bool; species::Int=1)
    @assert sum(getfield.(seq.Events, :EventType) .== :ADC) > 0 "💥 No ADC event found in CESTSequence. Aborting."
    println("Running detailed version")
    Mz = [];
    Mx = [];
    My = [];
    ts = [];

    Zspec = similar(params.Δω)
    M0_0  = create_M0(params.model)
    p0    = create_params(params.model, params.sys, 0.0)
    nadc  = 0
    for n = eachindex(seq.Events)
        global M0

        if n == 1
            continue
        elseif n == 2
            M0 = M0_0
        end

        # event = deepcopy(seq.Events[n])
        if seq.Events[n].EventType === :RF      
            p1 = create_params(params.model, params.sys, seq.Events[n].Δω)
            if isa(seq.Events[n].RF, ComplexF64)
                p = create_params(p1, seq.Events[n], :CW) # This will produce p::CESTParamsCW
            elseif isa(seq.Events[n].RF, Array{ComplexF64,1})
                p = create_params(p1, seq.Events[n], :Pulsed) # This will produce p::CESTParamsPulsed
            else
                println("uhOh")
            end
            prob = ODEProblem(solve_BME!, M0, (seq.Events[n].T[1], seq.Events[n].T[end]), p)
            sol  = solve(prob, eval(Meta.parse(params.solver)), reltol=params.reltol, abstol=params.abstol)
            M0   = sol.u[end];
            append!(Mz, getindex.(sol.u, params.model.nSpecies * 2 + species))
            append!(Mx, getindex.(sol.u, params.model.nSpecies * 0 + species))
            append!(My, getindex.(sol.u, params.model.nSpecies * 1 + species))
            append!(ts, sol.t)
        elseif seq.Events[n].EventType === :Delay || seq.Events[n].EventType === :ADC || seq.Events[n].EventType === :Spoiling
            p = create_params(p0, seq.Events[n], :Relax) # This will produce p::CESTParamsFP
            prob = ODEProblem(solve_BME!, M0, (seq.Events[n].T[1], seq.Events[n].T[end]), p)
            sol = solve(prob, eval(Meta.parse(params.solver)), reltol=1e-6, abstol=1e-6)
            M0 = sol.u[end];
            if seq.Events[n].EventType === :Spoiling
                M0[params.model.nSpecies * 0 + 1] = 0.0 
                M0[params.model.nSpecies * 1 + 1] = 0.0 
            elseif seq.Events[n].EventType === :ADC && resetM0
                M0 = M0_0
            end
            append!(Mz, getindex.(sol.u, params.model.nSpecies * 2 + species))
            append!(Mx, getindex.(sol.u, params.model.nSpecies * 0 + species))
            append!(My, getindex.(sol.u, params.model.nSpecies * 1 + species))
    
            append!(ts, sol.t)
        else # if :Ignore
            continue
        end
        if seq.Events[n].EventType === :ADC
            nadc += 1 # append!(dw, Δω)
            # println(nadc)
            Zspec[nadc] = getindex.(sol.u, params.model.nSpecies * 2 + 1)[end]
        end
    end
    # return Zspec, Mz, ts, Mx, My 
    return (; Zspec=Zspec, t=ts, Mx=Mx, My=My, Mz=Mz) 
end

# TODO: Create another method for simulate_CEST that uses EssembleProblem


# Method 1: Continuous wave
function solve_BME!(dM, M, p::CESTParamsCW, t)
    
    N = p.N                                                             
    Mx = @view M[1:N]
    My = @view M[N+1:2N]
    Mz = @view M[2N+1:3N]

    dMx = @view dM[1:N]
    dMy = @view dM[N+1:2N]
    dMz = @view dM[2N+1:3N]

    ωx = γ̄M * p.B1 * cos(p.ϕ) # make sure B₁ is [μT] and ϕ is [rad]
    ωy = γ̄M * p.B1 * sin(p.ϕ)
    
    R = p.R # semi-alloc Relaxation 
    C = p.C # semi-alloc Exchange
    @inbounds @simd for n = 1:N::Int64
        idx = 3(n-1) # * self-reminder: dont change idx to m or k as tempting as it is lol
        Ωⁿ = p.Ω[n]
        R[idx+1, idx+2] = Ωⁿ
        R[idx+2, idx+1] = -Ωⁿ
        R[idx+1, idx+3] = -ωy
        R[idx+2, idx+3] = ωx
        R[idx+3, idx+1] = ωy
        R[idx+3, idx+2] = -ωx

        dMx[n] = R[idx+1, idx+1] * Mx[n] + R[idx+1, idx+2] * My[n] + R[idx+1, idx+3] * Mz[n] + C[idx+1]
        dMy[n] = R[idx+2, idx+1] * Mx[n] + R[idx+2, idx+2] * My[n] + R[idx+2, idx+3] * Mz[n] + C[idx+2]
        dMz[n] = R[idx+3, idx+1] * Mx[n] + R[idx+3, idx+2] * My[n] + R[idx+3, idx+3] * Mz[n] + C[idx+3]
    end

    @inbounds @simd for n = 1:N::Int64
        for idx = 1:N::Int64
            dMx[n] += p.Kx[n, idx] * Mx[idx]
            dMy[n] += p.Ky[n, idx] * My[idx]
            dMz[n] += p.Kz[n, idx] * Mz[idx]
        end
    end

    return nothing
end

# Method 2: Free precession
function solve_BME!(dM, M, p::CESTParamsFP, t)
    
    N = p.N
    Mx = @view M[1:N]
    My = @view M[N+1:2N]
    Mz = @view M[2N+1:3N]

    dMx = @view dM[1:N]
    dMy = @view dM[N+1:2N]
    dMz = @view dM[2N+1:3N]
    
    R = p.R # semi-alloc Relaxation 
    C = p.C # semi-alloc Exchange
    @inbounds @simd for n = 1:N::Int64
        idx = 3(n-1) # * self-reminder: dont change idx to m or k as tempting as it is lol
        Ωⁿ = p.Ω[n]
        R[idx+1, idx+2] = Ωⁿ
        R[idx+2, idx+1] = -Ωⁿ
        
        dMx[n] = R[idx+1, idx+1] * Mx[n] + R[idx+1, idx+2] * My[n] + C[idx+1]
        dMy[n] = R[idx+2, idx+1] * Mx[n] + R[idx+2, idx+2] * My[n] + C[idx+2]
        dMz[n] = R[idx+3, idx+3] * Mz[n] + C[idx+3]
    end

    @inbounds @simd for n = 1:N::Int64
        for idx = 1:N::Int64
            dMx[n] += p.Kx[n, idx] * Mx[idx]
            dMy[n] += p.Ky[n, idx] * My[idx]
            dMz[n] += p.Kz[n, idx] * Mz[idx]
        end
    end

    return nothing
end

# Method 3: Varying B1 amplitude+phase
function solve_BME!(dM, M, p::CESTParamsPulsed, t)
    
    N = p.N::Int64
    Mx = @view M[1:N]
    My = @view M[N+1:2N]
    Mz = @view M[2N+1:3N]

    dMx = @view dM[1:N]
    dMy = @view dM[N+1:2N]
    dMz = @view dM[2N+1:3N]

    ωx = γ̄M * p.B1(t) * cos(p.ϕ(t))::Float64 # make sure B₁ is [μT] and ϕ is [rad]
    ωy = γ̄M * p.B1(t) * sin(p.ϕ(t))::Float64
    
    R = p.R # semi-alloc Relaxation 
    C = p.C # semi-alloc Exchange
    @inbounds @simd for n = 1:N::Int64
        idx = 3(n-1)::Int64 # * self-reminder: dont change idx to m or k as tempting as it is lol
        Ωⁿ = p.Ω[n]
        R[idx+1, idx+2] = Ωⁿ
        R[idx+2, idx+1] = -Ωⁿ
        R[idx+1, idx+3] = -ωy
        R[idx+2, idx+3] = ωx
        R[idx+3, idx+1] = ωy
        R[idx+3, idx+2] = -ωx

        dMx[n] = R[idx+1, idx+1] * Mx[n] + R[idx+1, idx+2] * My[n] + R[idx+1, idx+3] * Mz[n] + C[idx+1]
        dMy[n] = R[idx+2, idx+1] * Mx[n] + R[idx+2, idx+2] * My[n] + R[idx+2, idx+3] * Mz[n] + C[idx+2]
        dMz[n] = R[idx+3, idx+1] * Mx[n] + R[idx+3, idx+2] * My[n] + R[idx+3, idx+3] * Mz[n] + C[idx+3]
    end

    @inbounds @simd for n = 1:N::Int64
        for idx = 1:N::Int64
            dMx[n] += p.Kx[n, idx] * Mx[idx]
            dMy[n] += p.Ky[n, idx] * My[idx]
            dMz[n] += p.Kz[n, idx] * Mz[idx]
        end
    end

    return nothing
end

# Method 4: NamedTuple

# function solve_BME!(dM, M, p, t)
    
#     N = p.N
#     Mx = @view M[1:N]
#     My = @view M[N+1:2N]
#     Mz = @view M[2N+1:3N]

#     dMx = @view dM[1:N]
#     dMy = @view dM[N+1:2N]
#     dMz = @view dM[2N+1:3N]

#     isa(p.ϕ, Number) ? ϕ = p.ϕ : ϕ = p.ϕ(t)
#     isa(p.b1, Number) ? B₁ = p.b1 : B₁ = p.b1(t)

#     ωx = γ̄M * B₁ * cos(ϕ) # make sure B₁ is [μT] and ϕ is [rad]
#     ωy = γ̄M * B₁ * sin(ϕ)
    
#     R = p.R # semi-alloc Relaxation 
#     C = p.C # semi-alloc Exchange
#     @inbounds @simd for n in 1:N
#         idx = 3 * (n - 1) # * self-reminder: dont change idx to m or k as tempting as it is lol
#         Ωⁿ = p.Ω[n]
#         R[idx+1, idx+2] = Ωⁿ
#         R[idx+2, idx+1] = -Ωⁿ
#         R[idx+1, idx+3] = -ωy
#         R[idx+2, idx+3] = ωx
#         R[idx+3, idx+1] = ωy
#         R[idx+3, idx+2] = -ωx

#         dMx[n] = R[idx+1, idx+1] * Mx[n] + R[idx+1, idx+2] * My[n] + R[idx+1, idx+3] * Mz[n] + C[idx+1]
#         dMy[n] = R[idx+2, idx+1] * Mx[n] + R[idx+2, idx+2] * My[n] + R[idx+2, idx+3] * Mz[n] + C[idx+2]
#         dMz[n] = R[idx+3, idx+1] * Mx[n] + R[idx+3, idx+2] * My[n] + R[idx+3, idx+3] * Mz[n] + C[idx+3]
#     end

#     @inbounds @simd for n in 1:N
#         for idx in 1:N
#             dMx[n] += p.Kx[n, idx] * Mx[idx]
#             dMy[n] += p.Ky[n, idx] * My[idx]
#             dMz[n] += p.Kz[n, idx] * Mz[idx]
#         end
#     end

#     return nothing
# end

