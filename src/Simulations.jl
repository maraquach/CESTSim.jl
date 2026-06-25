"""Simulation support"""

# include("Aux.jl")
using DifferentialEquations
using Interpolations


function simulate_CEST(seq::CESTSequence, params::CESTSimParams, resetM0::Bool=true)
    Zspec = copy(params.Δω)
    M0_0 = create_M0(params.model)
    p0 = create_params(params.model, params.sys, 0.0)
    global nadc
    global p
    global event
    nadc=0
    for n = eachindex(seq.Events)
        global M0
        if n == 1
            continue
        elseif n == 2
            M0 = M0_0
        end
        event = deepcopy(seq.Events[n])
        if event.EventType === :RF      
            p1 = create_params(params.model, params.sys, event.Δω)
            if isa(event.RF, ComplexF64)
                p = create_params(p1, event, :CW) # This will produce p::CESTParamsCW
            elseif isa(event.RF, Array{ComplexF64,1})
                p = create_params(p1, event, :Pulsed) # This will produce p::CESTParamsPulsed
            else
                println("uhOh")
            end
            prob = ODEProblem(solve_BME!, M0, (event.T[1], event.T[end]), p) 
            sol = solve(prob, eval(Meta.parse(params.solver)), dense=false, save_everystep=false)
            M0 = sol.u[end];
        elseif event.EventType === :Delay || event.EventType === :ADC || event.EventType === :Spoiling
            p = create_params(p0, event, :Relax) # This will produce p::CESTParamsFP
            prob = ODEProblem(solve_BME!, M0, (event.T[1], event.T[end]), p)
            sol = solve(prob, eval(Meta.parse(params.solver)), dense=false, save_everystep=false)
            M0 = sol.u[end];
            if event.EventType === :Spoiling
                M0[params.model.nSpecies * 0 + 1] = 0.0 
                M0[params.model.nSpecies * 1 + 1] = 0.0 
            elseif event.EventType === :ADC && resetM0
                M0 = M0_0
            end
        else # if :Ignore
            continue
        end
        if event.EventType === :ADC
            nadc += 1 # append!(dw, Δω)
            # println(nadc)
            Zspec[nadc] = getindex.(sol.u, params.model.nSpecies * 2 + 1)[end]
        end
    end
    return Zspec
end


function simulate_CEST(seq::CESTSequence, params::CESTSimParams; detailed::Bool=true)
    Mz = [];
    Mx = [];
    My = [];
    ts = [];

    Zspec = copy(params.Δω)
    M0_0 = create_M0(params.model)
    p0 = create_params(params.model, params.sys, 0.0)
    global nadc
    global p
    global event
    nadc=0
    for n = eachindex(seq.Events)
        global M0
        if n == 1
            continue
        elseif n == 2
            M0 = M0_0
        end
        event = deepcopy(seq.Events[n])
        if event.EventType === :RF      
            p1 = create_params(params.model, params.sys, event.Δω)
            if isa(event.RF, ComplexF64)
                p = create_params(p1, event, :CW) # This will produce p::CESTParamsCW
            elseif isa(event.RF, Array{ComplexF64,1})
                p = create_params(p1, event, :Pulsed) # This will produce p::CESTParamsPulsed
            else
                println("uhOh")
            end
            prob = ODEProblem(solve_BME!, M0, (event.T[1], event.T[end]), p) 
            sol = solve(prob, eval(Meta.parse(params.solver)))
            M0 = sol.u[end];
        elseif event.EventType === :Delay || event.EventType === :ADC || event.EventType === :Spoiling
            p = create_params(p0, event, :Relax) # This will produce p::CESTParamsFP
            prob = ODEProblem(solve_BME!, M0, (event.T[1], event.T[end]), p)
            sol = solve(prob, eval(Meta.parse(params.solver)))
            M0 = sol.u[end];
            if event.EventType === :Spoiling
                M0[params.model.nSpecies * 0 + 1] = 0.0 
                M0[params.model.nSpecies * 1 + 1] = 0.0 
            end
            append!(Mz, getindex.(sol.u, params.model.nSpecies * 2 + 1))
            append!(Mx, getindex.(sol.u, params.model.nSpecies * 0 + 1))
            append!(My, getindex.(sol.u, params.model.nSpecies * 1 + 1))
            append!(ts, sol.t)
        else # if :Ignore
            continue
        end
        if event.EventType === :ADC
            nadc += 1 # append!(dw, Δω)
            # println(nadc)
            Zspec[nadc] = getindex.(sol.u, params.model.nSpecies * 2 + 1)[end]
        end
    end
    return Zspec, Mz, ts, Mx, My 
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

function solve_BME!(dM, M, p, t)
    
    N = p.N
    Mx = @view M[1:N]
    My = @view M[N+1:2N]
    Mz = @view M[2N+1:3N]

    dMx = @view dM[1:N]
    dMy = @view dM[N+1:2N]
    dMz = @view dM[2N+1:3N]

    isa(p.ϕ, Number) ? ϕ = p.ϕ : ϕ = p.ϕ(t)
    isa(p.b1, Number) ? B₁ = p.b1 : B₁ = p.b1(t)

    ωx = γ̄M * B₁ * cos(ϕ) # make sure B₁ is [μT] and ϕ is [rad]
    ωy = γ̄M * B₁ * sin(ϕ)
    
    R = p.R # semi-alloc Relaxation 
    C = p.C # semi-alloc Exchange
    @inbounds @simd for n in 1:N
        idx = 3 * (n - 1) # * self-reminder: dont change idx to m or k as tempting as it is lol
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

    @inbounds @simd for n in 1:N
        for idx in 1:N
            dMx[n] += p.Kx[n, idx] * Mx[idx]
            dMy[n] += p.Ky[n, idx] * My[idx]
            dMz[n] += p.Kz[n, idx] * Mz[idx]
        end
    end

    return nothing
end


# function simulate_CEST(seq::CESTSequence, params::CESTSimParams)
#     Mz = [];
#     Mx = [];
#     My = [];
#     ts = [];
#     Zspec = []
#     dw = []
#     M0_0 = create_M0(params.model)
#     for n = eachindex(seq.Events)
#         global M0
#         global Δω
#         if n == 1
#             continue
#         elseif n == 2
#             M0 = M0_0
#         end
#         event = deepcopy(seq.Events[n])
#         if event.EventType === :RF
#             p = create_params(params.model, params.sys, event.Δω)
#             fB₁, fϕ = interp_b1(event)
#             prob = ODEProblem(solve_BME!, M0, (event.T[1], event.T[end]), merge(p, (; b1=fB₁, ϕ=fϕ)))
#             sol = solve(prob, eval(Meta.parse(params.solver)))
#             M0 = sol.u[end];
#             Δω = event.Δω
#         elseif event.EventType === :Delay || event.EventType === :ADC || event.EventType === :Spoiling
#             p = create_params(params.model, params.sys, 0.0)
#             prob = ODEProblem(solve_BME!, M0, (event.T[1], event.T[end]), merge(p, (; b1=0.0, ϕ=0.0)))
#             sol = solve(prob, eval(Meta.parse(params.solver)))
#             M0 = sol.u[end];
#             if event.EventType === :Spoiling
#                 M0[params.model.nSpecies * 0 + 1] = 0.0 
#                 M0[params.model.nSpecies * 1 + 1] = 0.0 
#             end
#         else # if :Ignore
#             continue
#         end
#         append!(Mz, getindex.(sol.u, params.model.nSpecies * 2 + 1))
#         append!(Mx, getindex.(sol.u, params.model.nSpecies * 0 + 1))
#         append!(My, getindex.(sol.u, params.model.nSpecies * 1 + 1))
#         append!(ts, sol.t)
        
#         if event.EventType === :ADC
#             append!(dw, Δω)
#             append!(Zspec, Mz[end])
#         end
#     end
#     return Zspec, dw, Mz, ts, Mx, My 
# end