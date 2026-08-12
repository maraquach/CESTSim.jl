
module pulses
    using DSP

    function calcRFpower(rf, T, Δt::Float64=1e-6) # Translated from Pulseq -- https://github.com/pulseq/pulseq
        # adapted from Pulseq
        rfs_sq = rf .* conj(rf)
        total_energy = sum(rfs_sq) * Δt
        peak_pwr = maximum(rfs_sq)
        rf_rms = sqrt(total_energy/(T[end]-T[1]))
        output = (; totalEnergy = total_energy, peakPower = peak_pwr, rfRMS=rf_rms)   
    return output

    end
    function rect(A::Float64, T::Union{Array{Float64,1},AbstractRange{Float64}}, ϕ::Float64=0.0)
        Λ = A .+ zeros(length(T))
        pulse = Λ.*exp(ϕ*im)
        return pulse
    end

    function rect(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, ϕ::Float64=0.0)
        N = floor((tend - t0)/Δt)|> Int
        newt = LinRange(t0, tend, N)
        # Λ = repeat([A], N)
        Λ = A .+ zeros(N)
        pulse = Λ.*exp(ϕ*im)
        return pulse, newt
    end

    # * Gaussian

    function gaussian(A::Float64, T::Array{Float64,1}, ϕ::Float64=0.0; σ::Float64=7.434, windowed::Symbol=:none)
        N = length(T)
        σ = (T[end]-T[1])/σ
        t0 = T[1]; tend = T[end];
        Λ = A .* exp.(-(T .- (tend-t0)/2 .- t0).^2/(2σ^2)) 
        if windowed == :hamming || windowed == :Hamming
            🪟 = DSP.hamming(N)
            Λ = Λ .* 🪟
        elseif windowed == :hanning || windowed == :Hanning
            🪟 = DSP.hanning(N)
            Λ = Λ .* 🪟
        end
        pulse = Λ.*exp(ϕ*im)
        return pulse
    end

    function gaussian(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, ϕ::Float64=0.0; σ::Float64=7.434, windowed::Symbol=:none)
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        σ = (tend-t0)/σ
        Λ = A .* exp.(-(T .- (tend-t0)/2 .- t0).^2/(2σ^2)) 
        if windowed == :hamming || windowed == :Hamming
            🪟 = DSP.hamming(N)
            Λ = Λ .* 🪟
        elseif windowed == :hanning || windowed == :Hanning
            🪟 = DSP.hanning(N)
            Λ = Λ .* 🪟
        end
        pulse =  Λ.*exp(ϕ*im)
        return pulse, T
    end

    # * Hanning 

    function hanning(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, ϕ::Float64=0.0)
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        Λ = A .* DSP.hanning(N)
        pulse = Λ.*exp(ϕ*1im)
        return pulse, T
    end
    function hanning(A::Float64, T::Array{Float64,1}, ϕ::Float64=0.0)
        N = length(T)
        Λ = A .* DSP.hanning(N)
        pulse = Λ.*exp(ϕ*1im)
        return pulse, T
    end

    # * Hamming 

    function hamming(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, ϕ::Float64=0.0)
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        Λ = A .* DSP.hamming(N)
        pulse = Λ.*exp(ϕ*1im)
        return pulse, T
    end
    function hamming(A::Float64, T::Array{Float64,1}, ϕ::Float64=0.0)
        N = length(T)
        Λ = A .* DSP.hamming(N)
        pulse = Λ.*exp(ϕ*1im)
        return pulse, T
    end

    # * Sinc
    function sinc(A::Float64, T::Array{Float64,1}, ϕ::Float64=0.0; TBP::Float64=4.0, apod::Union{Symbol, Float64}=:none)

        if      apod == :Hamming;   α = 0.46;
        elseif  apod == :Hanning;   α = 0.5;
        elseif  apod == :none;      α = 0.0;
        else;                       α = apod;                
        end 

        τ = T .- (T[end] - T[1])/2 .- T[1]
        t0 = (τ[end] - τ[1])/TBP
        Λ = ((t0.*sin.(π.*τ ./ t0)) ./ (π.*τ)) .* ((1-α)+α.*cos.((π.*τ)./(TBP*t0)))
        replace!(Λ, NaN => maximum(Λ))
        pulse = A.*Λ .* exp(ϕ*im)
        return pulse
    end

    function sinc(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, ϕ::Float64=0.0; TBP::Float64=4.0, apod::Union{Symbol, Float64}=:none)

        if      apod == :Hamming;   α = 0.46;
        elseif  apod == :Hanning;   α = 0.5;
        elseif  apod == :none;      α = 0.0;
        else;                       α = apod;                
        end 

        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)

        τ = T .- (tend - t0)/2 .- t0
        t0 = (τ[end] - τ[1])/TBP
        Λ = ((t0.*sin.(π.*τ ./ t0)) ./ (π.*τ)) .* ((1-α).+α.*cos.((π.*τ)./(TBP*t0)))
        replace!(Λ, NaN => maximum(Λ))
        # Λ[1] = 0.0; Λ[end] = 0.0
        pulse = A.*Λ .* exp(ϕ*im)
        return pulse, T
    end

    # * Fermi
    function fermi(A::Float64, T::Array{Float64,1}, ϕ::Float64=0.0; PW::Float64=.25*(T[end]-T[1]))

        tw = ((T[end] - T[1]) - 2 * PW)/13.81
        Λ = A ./ (1 .+ exp.((abs.(T .- (T[end] - T[1])/2 .- T[1]) .- PW) ./ tw))
        pulse = Λ .* exp(ϕ*im)
        return pulse
    end

    function fermi(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, ϕ::Float64=0.0; PW::Float64=.25*(tend-t0))

        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        
        tw = ((tend - t0) - 2 * PW)/13.81
        Λ = A ./ (1 .+ exp.((abs.(T .- (tend - t0)/2 .- t0) .- PW) ./ tw))
        pulse = Λ .* exp(ϕ*im)

        return pulse, T
    end

    # & Adiabatic pulses
    # * HS - RF spectral BW = μβ ; b1max >> (β√μ)/γ threshold for adiabacity; this is no longer similar to MRIlab
    
    function HS(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=100.0, μ::Float64=25.5, ef::Float64=3.5, twindow::Float64=2.5e-3)
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
 
        tp = tend - t0;

        if type == :full
            Λ = @. A * sech(π * β * (T - tp/2))
            ω = @. (β*μ*π) * tanh((β*π) * (T - tp/2))
        elseif type == :half
            Λ = @. A * sech(π * β * (T - tp))
            ω = @. (β*μ*π) * tanh((β*π) * (T - tp))
            ω = ω .- ω[findmin(abs.(ω))[2]]
        end
        ϕ = cumsum(ω) .* (tp / length(T))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end
        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T, ω
    end    
    
    function HSn(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; n::Int=8, β::Float64=15.0, μ::Float64=25.1, SLAC::Bool=false, ξ::Float64=1.0)
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        tp = tend - t0;

        if type == :full
            Λ = @. sech((π * β * (T-tp/2))^n)
            ω = μ*β*π*(cumsum(Λ.^2) ./ sum(Λ.^2) .- .5)
            ω = ω .- ω[findmin(abs.(ω))[2]]
            Λ = @. A .* Λ
        elseif type == :half
            Λ = @. sech((π * β * (T-tp))^n)
            ω = μ*β*π*(cumsum(Λ.^2) ./ sum(Λ.^2) .- 1)
            ω = ω .- ω[findmin(abs.(ω))[2]]
            Λ = @. A .* Λ
        end

        ϕ = cumsum(ω) .* Δt

        if SLAC    
            γ = 42.577478461
            α = @. atan(2π*γ*Λ/ω)
            α[end] = α[end-1]

            dαdt = [0.0; diff(α)] ./ Δt
            dαdt[1] = dαdt[2]
            dαdt[end] = dαdt[end-1]

            ΛSL = @. -1/(2π*γ*ξ) * dαdt

            xAdia = @. Λ * cos(ϕ)
            yAdia = @. Λ * sin(ϕ)
            xSL   = @. ΛSL * cos(ϕ-π/2)
            ySL   = @. ΛSL * sin(ϕ-π/2)

            x = @. xAdia + xSL
            y = @. yAdia + ySL

            Λ = @. sqrt(x^2 + y^2)
            ϕ = @. atan(y,x)
            ϕ = unwrap(ϕ) .- maximum(unwrap(ϕ))
        end
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end
        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T, ω
    end    
    
    
    function HS_AM(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=5.3, b::Float64=157.0, Pe::Float64=.5)
        # Based on Warnking and Pike, 2004: https://doi.org/10.1002/mrm.20262 
        γ = 42.577478461 * 2π
        a = sqrt(b^2 - 4*(acosh(cosh(b/2)*sqrt(1-Pe)))^2)
        μ = b/π

        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        
        tp = tend - t0;
        T2 = T .- tp/2
        
        Λ = @. (a/(2π)) * (4*β/tp) * sech(4*β*T2/tp/2) / γ
        ω = @. (b/(2π)) * (4*β/tp) * tanh(4*β*T2/tp/2)
        ϕ = cumsum(ω) .* (tp / length(T))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        
        if type == :half && tip == :down
            Λ = Λ[1:Int(N/2)]
            ϕ = ϕ[1:Int(N/2)]
            T = T[1:Int(N/2)]
        elseif type == :half && tip == :up
            Λ = Λ[Int(N/2+1):end]
            ϕ = reverse(ϕ[Int(N/2+1):end])
            T = T[1:Int(N/2)]
        end
        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T, ω
    end    


    function HS_SLAC(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=100.0, μ::Float64=25.5, ξ::Float64=1.0)
        # Based on Green et al., 2019: https://doi.org/10.1016/j.jmr.2019.106595
        γ = 42.577478461
        
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        tp = tend - t0
        Λ = @. ξ * A * sech(π * β * (T - tp))
        ω = @. (β*μ*π) * tanh((β*π) * (T - tp))
        ω = ω .- ω[findmin(abs.(ω))[2]]

        ϕ = cumsum(ω) .* Δt
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        α = @. atan(2π*γ*Λ/ω)
        α[end] = α[end-1]

        dαdt = [0.0; diff(α)] ./ Δt
        dαdt[1] = dαdt[2]
        dαdt[end] = dαdt[end-1]

        ΛSL = @. -1/(2π*γ*ξ) * dαdt

        xAdia = @. Λ * cos(ϕ)
        yAdia = @. Λ * sin(ϕ)
        xSL   = @. ΛSL * cos(ϕ-π/2)
        ySL   = @. ΛSL * sin(ϕ-π/2)

        x = @. xAdia + xSL
        y = @. yAdia + ySL

        Λ = @. sqrt(x^2 + y^2)
        ϕ = @. atan(y,x)
        ϕ = unwrap(ϕ) .- maximum(unwrap(ϕ))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π

        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end

        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T

    end


    function HSExp(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=100.0, μ::Float64=25.5, ef::Float64=3.5, twindow::Float64=2.5e-3)
        # Based on: Herz et al., 2018: https://doi.org/10.1002/mrm.27380
        
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        
        H = @. .42 - .5*cos(π * T[T .< twindow]/twindow) + .08*cos(2*π*T[T .< twindow]/twindow)
        H[1] = 0.0
        tp = tend - t0;

        if type == :full
            Λ = @. A * sech(π * β * (T - tp/2))
            Λ[T .< twindow] = Λ[T .< twindow] .* H
            ω = @. -(β*μ*π) * exp(-T/(tp/2) * ef)
            # ω = @. ω - ω[findmin(abs.(ω))[2]]
        elseif type == :half
            Λ = @. A * sech(π * β * (T - tp))
            Λ[T .< twindow] = Λ[T .< twindow] .* H
            ω = @. -(β*μ*π) * exp(-T/tp * ef)
            ω = ω .- ω[findmin(abs.(ω))[2]]
        end
        ϕ = cumsum(ω) .* (tp / length(T))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end
        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T
    end    
    # * BIR

    function BIR(A::Float64, T::Array{Float64,1}, type::Symbol=:BIR4; β::Float64=1.4711, λ::Float64=10.0, F::Float64=10e3)

        tend = T[end]-T[1]; t0 = 0.0
        T = LinRange(t0, tend, length(T))
        if type == :BIR1
            ζ = π / (tend-t0);
            τ1 = T[T.<(tend-t0)/2]
            ϕ1 = F .* sin.(ζ .* τ1)
            
            τ2 = T[T.>=(tend-t0)/2]
            ϕ2 = -F .* sin.(ζ .* τ2) .+ (π/2)

            Λ = A .* cos.(ζ .* T)
            ϕ = [ϕ1; ϕ2]

        elseif type == :BIR2
            ζ = 2π / (tend - t0)
            τ1 = T[T.<(tend-t0)/4]
            τ2 = T[T.<(tend-t0)/2 .&& T.>=(tend-t0)/4]
            τ3 = T[T.>=(tend-t0)/2]
            
            ϕ1 = abs.(F .* sin.(ζ.*τ1)) 
            ϕ2 = abs.(F .* sin.(ζ.*τ2)) .+ π/2
            ϕ3 = abs.(F .* sin.(ζ.*τ3)) .+ (π/2+π)

            Λ = abs.(A .* cos.(ζ .* T))
            ϕ = [ϕ1; ϕ2; ϕ3]
        elseif type == :BIR4
            τ1 = T[T.<(tend-t0)/4];
            Λ1 = A .* tanh.(λ.*(1 .- 4 .* τ1./(tend-t0)));
            ϕ1 = tan.(4 .* β .* τ1./(tend-t0))./(tan(β)*2π);
            τ2 = T[T.>=(tend-t0)/4 .&& T.<(tend-t0)/2];
            Λ2 = A .* tanh.(λ .* (4 .* τ2 ./(tend-t0) .- 1));
            ϕ2 = tan.(β .* (4 .* τ2./(tend-t0).-2)) ./ (tan(β)*2π);
            τ3 = T[T.>=(tend-t0)/2 .&& T.<3*(tend-t0)/4];
            Λ3 = A .* tanh.(λ .* (3 .- 4 .* τ3 ./ (tend-t0)));
            ϕ3 = tan.(β .* (4 .* τ3 ./ (tend-t0) .- 2)) ./ (tan(β)*2π);
            τ4 = T[T .>= 3*(tend-t0)/4];
            Λ4 = A .*tanh.(λ .* (4 .* τ4./(tend-t0) .- 3));
            ϕ4 = tan.(β .* (4 .* τ4 ./ (tend-t0) .- 4)) ./ (tan(β) * 2π);
            
            Λ = [Λ1; Λ2; Λ3; Λ4]
            ϕ = [ϕ1; ϕ2; ϕ3; ϕ4]
        end
        pulse = Λ .* exp.(ϕ.*im)
        return pulse
    end

    function BIR(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:BIR4; β::Float64=1.4711, λ::Float64=10.0, F::Float64=10e3)

        N = floor((tend - t0)/Δt)|> Int
        tend = tend - t0; t0 = 0.0
        T = LinRange(t0, tend, N)

        if type == :BIR1
            ζ = π / (tend-t0);
            τ1 = T[T.<(tend-t0)/2]
            ϕ1 = F .* sin.(ζ .* τ1)
            
            τ2 = T[T.>=(tend-t0)/2]
            ϕ2 = -F .* sin.(ζ .* τ2) .+ (π/2)

            Λ = A .* cos.(ζ .* T)
            ϕ = [ϕ1; ϕ2]

        elseif type == :BIR2
            ζ = 2π / (tend - t0)
            τ1 = T[T.<(tend-t0)/4]
            τ2 = T[T.<(tend-t0)/2 .&& T.>=(tend-t0)/4]
            τ3 = T[T.>=(tend-t0)/2]
            
            ϕ1 = abs.(F .* sin.(ζ.*τ1))
            ϕ2 = abs.(F .* sin.(ζ.*τ2)) .+ π/2
            ϕ3 = abs.(F .* sin.(ζ.*τ3)) .+ (π/2+π)

            Λ = abs.(A .* cos.(ζ .* T))
            ϕ = [ϕ1; ϕ2; ϕ3]
        elseif type == :BIR4
            τ1 = T[T.<(tend-t0)/4];
            Λ1 = A .* tanh.(λ.*(1 .- 4 .* τ1./(tend-t0)));
            ϕ1 = tan.(4 .* β .* τ1./(tend-t0))./(tan(β)*2π);
            τ2 = T[T.>=(tend-t0)/4 .&& T.<(tend-t0)/2];
            Λ2 = A .* tanh.(λ .* (4 .* τ2 ./(tend-t0) .- 1));
            ϕ2 = tan.(β .* (4 .* τ2./(tend-t0).-2)) ./ (tan(β)*2π);
            τ3 = T[T.>=(tend-t0)/2 .&& T.<3*(tend-t0)/4];
            Λ3 = A .* tanh.(λ .* (3 .- 4 .* τ3 ./ (tend-t0)));
            ϕ3 = tan.(β .* (4 .* τ3 ./ (tend-t0) .- 2)) ./ (tan(β)*2π);
            τ4 = T[T .>= 3*(tend-t0)/4];
            Λ4 = A .*tanh.(λ .* (4 .* τ4./(tend-t0) .- 3));
            ϕ4 = tan.(β .* (4 .* τ4 ./ (tend-t0) .- 4)) ./ (tan(β) * 2π);
            
            Λ = [Λ1; Λ2; Λ3; Λ4]
            ϕ = [ϕ1; ϕ2; ϕ3; ϕ4]
        end
        pulse = Λ .* exp.(ϕ.*im)
        T = T .+ t0
        return pulse, T
    end


function tanhtan(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, tip::Symbol=:down, Δ⁻::Symbol=:pos; Δf::Float64=7000.0, ξ::Float64=10.0, κ::Float64=1.537, SLAC::Bool=false, SLACξ=1.3)
    tp = tend - t0
    N = floor((tend - t0)/Δt)|> Int
    T = LinRange(t0, tend, N)
    Λ = @. A * tanh(ξ * T / tp)
    ω = @. Δf * (tan(κ*(T/tp - 1))/tan(κ))
    ϕ = cumsum(ω) .* Δt

    if SLAC    
        γ = 42.577478461
        α = @. atan(2π*γ*Λ/ω)
        α[end] = α[end-1]

        dαdt = [0.0; diff(α)] ./ Δt
        dαdt[1] = dαdt[2]
        dαdt[end] = dαdt[end-1]

        ΛSL = @. -1/(2π*γ*SLACξ) * dαdt

        xAdia = @. Λ * cos(ϕ)
        yAdia = @. Λ * sin(ϕ)
        xSL   = @. ΛSL * cos(ϕ-π/2)
        ySL   = @. ΛSL * sin(ϕ-π/2)

        x = @. xAdia + xSL
        y = @. yAdia + ySL

        Λ = @. sqrt(x^2 + y^2)
        ϕ = @. atan(y,x)
        ϕ = unwrap(ϕ) .- maximum(unwrap(ϕ))
    end
    ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
    if Δ⁻ == :neg
        ϕ = 0.0 .- ϕ
    end
    if tip == :down
        Λ = Λ
        ϕ = ϕ
    elseif tip == :up            
        Λ = reverse(Λ)
        ϕ = reverse(0.0 .- ϕ)
    end
    pulse = Λ .* exp.(ϕ .* 1im)
    return pulse, T
end

function unwrap(ϕ)
    unwrapped = copy(ϕ)
    for n = 2:length(ϕ)
        while unwrapped[n] - unwrapped[n-1] >= π
            unwrapped[n] -= 2π
        end
        while unwrapped[n] - unwrapped[n-1] <= -π
            unwrapped[n] += 2π
        end
    end
    return unwrapped
end

end

module pulsesWIP


    function HSmod(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=100.0, μ::Float64=25.5, ef::Float64=3.5, twindow::Float64=2.5e-3)
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        
        tp = tend - t0;

        if type == :full
            Λ = @. A * sech(π * β * (T - tp/2))
            ω = @. (β*μ*π) * tanh((β*π) * (T - tp/2))
            ω = @. ω - ω[findmin(abs.(ω))[2]]
        elseif type == :half
            Λ = @. A * sech(π * β * (T - tp))
            # ω = @. (β*μ*π) * tanh((β*π) * (T - tp))
            ω = @. -(β*μ*π) * exp(-T/tp * ef)
            ω = ω .- ω[findmin(abs.(ω))[2]]
        end
        ϕ = cumsum(ω) .* (tp / length(T))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end
        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T
    end    
    function HS_SLAC_mod(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=100.0, μ::Float64=25.5, ξ::Float64=1.0, ef = 3.5)
        γ = 42.577478461
        
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)
        tp = tend - t0
        Λ = @. ξ * A * sech(π * β * (T - tp))
        ω = @. -(β*μ*π) * exp(-T/(tp/2) * ef)
        # ω = ω .- ω[findmin(abs.(ω))[2]]

        ϕ = cumsum(ω) .* Δt
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π
        
        α = @. atan(2π*γ*Λ/ω)
        α[end] = α[end-1]

        dαdt = [0.0; diff(α)] ./ Δt
        dαdt[1] = dαdt[2]
        dαdt[end] = dαdt[end-1]

        ΛSL = @. -1/(2π*γ*ξ) * dαdt

        xAdia = @. Λ * cos(ϕ)
        yAdia = @. Λ * sin(ϕ)
        xSL   = @. ΛSL * cos(ϕ-π/2)
        ySL   = @. ΛSL * sin(ϕ-π/2)

        x = @. xAdia + xSL
        y = @. yAdia + ySL

        Λ = @. sqrt(x^2 + y^2)
        ϕ = @. atan(y,x)
        ϕ = unwrap(ϕ) .- maximum(unwrap(ϕ))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π

        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end

        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T

    end

    function HSExp_SLAC(A::Float64, t0::Float64, tend::Float64, Δt::Float64=1e-6, type::Symbol=:half, tip::Symbol=:down, Δ⁻::Symbol=:pos; β::Float64=100.0, μ::Float64=25.5, ξ::Float64=1.0, twindow=2.5e-3, ef = 3.5)
        γ = 42.577478461
        
        N = floor((tend - t0)/Δt)|> Int
        T = LinRange(t0, tend, N)

        H = @. .42 - .5*cos(π * T[T .< twindow]/twindow) + .08*cos(2*π*T[T .< twindow]/twindow)
        H[1] = 0.0

        tp = tend - t0
        Λ = @. ξ * A * sech(π * β * (T - tp)) 
        Λ[T .< twindow] = Λ[T .< twindow] .* H
        ω = @. -(β*μ*π) * exp(-T/(tp/2) * ef)
        # ω = @. ω - ω[findmin(abs.(ω))[2]]
        # ω = @. (β*μ*π) * tanh((β*π) * (T - tp))
        ϕ = cumsum(ω) .* Δt
        
        α = @. atan(2π*γ*Λ/ω)
        α[end] = α[end-1]

        dαdt = [0.0; diff(α)] ./ Δt
        dαdt[1] = dαdt[2]
        dαdt[end] = dαdt[end-1]

        ΛSL = @. -1/(2π*γ*ξ) * dαdt

        xAdia = @. Λ * cos(ϕ)
        yAdia = @. Λ * sin(ϕ)
        xSL   = @. ΛSL * cos(ϕ-π/2)
        ySL   = @. ΛSL * sin(ϕ-π/2)

        x = @. xAdia + xSL
        y = @. yAdia + ySL

        Λ = @. sqrt(x^2 + y^2)
        ϕ = @. atan(y,x)
        ϕ = unwrap(ϕ) .- maximum(unwrap(ϕ))
        ϕ = mod.(ϕ .- ϕ[end] .+ π, 2π) .- π

        if Δ⁻ == :neg
            ϕ = 0.0 .- ϕ
        end
        if tip == :down
            Λ = Λ
            ϕ = ϕ
        elseif tip == :up            
            Λ = reverse(Λ)
            ϕ = reverse(0.0 .- ϕ)
        end
        
        pulse = Λ .* exp.(ϕ .* 1im)
        return pulse, T

    end

end