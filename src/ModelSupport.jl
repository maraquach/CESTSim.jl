
@kwdef mutable struct CESTModel
    path::Union{Nothing, String}
    names::Array{String}
    t1ⁿ::Array{Float64,1}
    t2ⁿ::Array{Float64,1}
    r1ⁿ::Array{Float64,1}
    r2ⁿ::Array{Float64,1}
    fⁿ::Array{Float64,1}
    Cⁿ::Union{Nothing, Array{Float64,1}}
    wⁿppm::Array{Float64,1}
    wⁿHz::Array{Float64,1}
    nSpecies::Int64
    maxRFsamples::Int64
end

""" 
    function load_model(modelpath::String)
    Load model from a YAML file as per PulseqCEST standards. Returns a struct of type CESTModel.
    Part of CESTSim.jl
        
    Arguments
        - modelpath::String: relative or absolute path to YAML model file
    Returns
        - ::CESTModel
"""
function load_model(modelpath::String) 
    # The returned parameter "model" is mutable struct -- be careful when making changes; otherwise feeding in a NamedTuple to solve_BME! also works.
    modeldict = YAML.load_file(modelpath)

    t1ⁿ     = [modeldict["water_pool"]["t1"]]
    t2ⁿ     = [modeldict["water_pool"]["t2"]]
    fⁿ      = [Float64(modeldict["water_pool"]["f"])]
    wⁿppm   = [Float64(modeldict["b0_inhom"])]
    Cⁿ      = []
    names   = ["water"]

    for (key, value) in modeldict["cest_pool"]
        # println(value["t1"])
        if !isnothing(value)
            # println(value["t1"])
            append!(t1ⁿ     , value["t1"])
            append!(t2ⁿ     , value["t2"])
            append!(fⁿ      , value["f"])
            append!(wⁿppm   , value["dw"])
            append!(Cⁿ      , value["k"])
            append!(names   , [key])
        end
    end

    fⁿ = fⁿ / fⁿ[1] # put in fraction relative to water concentration
    wⁿHz = wⁿppm .* (modeldict["b0"]*γ*1e-6) 

    model = CESTModel(
        path        = modelpath,
        names       = names,
        t1ⁿ         = t1ⁿ,
        t2ⁿ         = t2ⁿ,
        r1ⁿ         = 1.0 ./ t1ⁿ,
        r2ⁿ         = 1.0 ./ t2ⁿ,
        fⁿ          = fⁿ,
        Cⁿ          = Cⁿ,
        wⁿppm       = wⁿppm,
        wⁿHz       = wⁿHz,
        nSpecies    = size(t2ⁿ, 1) |> Int,
        maxRFsamples= modeldict["max_pulse_samples"] |> Int
    )

    return model
end

""" Method 1

    function create_model(template::String)
    Creates a ::CESTModel struct from a template. Current templates offered are `2PoolsGlu7T` or `2PoolsAmide7T`; More to be added
    Part of CESTSim.jl
    
    Arguments:
        - template::String: name of the template you would like.
    Returns:
        - ::CESTModel
"""
function create_model(template::String)
    if template === "2PoolsGlu7T"
        model = CESTModel(
            path    = "Private",
            t1ⁿ     = [1.8  ; 1.0],
            t2ⁿ     = [55e-3; 6.9e-3],
            r1ⁿ     = 1.0 ./ [1.8  ; 1.0],
            r2ⁿ     = 1.0 ./ [55e-3; 6.9e-3],
            fⁿ      = [111; 34.5e-3] ./ 111,
            Cⁿ      = [7400],
            wⁿppm   = [0; 3.2],
            wⁿHz   = [0; 3.2] .* (6.98 * γ * 1e-6),
            names   = ["water"; "glu"],
            nSpecies = 2,
            maxRFsamples = 200
        )    
    elseif template === "2PoolsAmide7T"
        model = CESTModel(
            path    = "Private",
            t1ⁿ     = [1.8  ; 1.0],
            t2ⁿ     = [55e-3; 10e-3],
            r1ⁿ     = 1.0 ./ [1.8  ; 1.0],
            r2ⁿ     = 1.0 ./ [55e-3; 10e-3],
            fⁿ      = [111; 72e-3] ./ 111,
            Cⁿ      = [20],
            wⁿppm   = [0; 3.5],
            wⁿHz   = [0; 3.5] .* (6.98 * γ * 1e-6),
            names   = ["water"; "amide"],
            nSpecies = 2,
            maxRFsamples = 200
        ) 
    else
        throw(ArgumentError("Template name not recognised. Try again with either `2PoolsGlu7T`, `2PoolsAmide7T`, or pass arguments to create your own model struct"))
    end
end

""" Method 2

    function create_model(names::Array{String,1},
                        T1::Array{Float64,1}, 
                        T2::Array{Float64,1}, 
                        f::Array{Float64,1}, 
                        k::Union{Nothing, Array{Float64,1}}, 
                        dw::Array{Float64,1},
                        B0::Float64, 
                        maxRFsamples::Int64=200)

    Explicitly create a ::CESTModel struct from input arguments. 
    Part of CESTSim.jl
    
    NOTE: please ensure the inputs are of the following units:
        - T1: [s] 
        - T2: [s] 
        - f: [M] or [a.u.] (≤ 1.0) 
        - k: [Hz] 
        - dw: [ppm]
        - B0: [T]
    Returns:
        - ::CESTModel
"""
function create_model(names::Array{String,1},
                        T1::Array{Float64,1}, 
                        T2::Array{Float64,1}, 
                        f::Array{Float64,1}, 
                        k::Union{Nothing, Array{Float64,1}}, 
                        dw::Array{Float64,1},
                        B0::Float64, 
                        maxRFsamples::Int64=200)
    # note that the input arguments are different than the CESTModule struct to follow more common notations e.g., k for exchange rates not C, dw not Δω for quick typing
    model = CESTModel(
        path    = "Private",
        names   = names,
        t1ⁿ     = T1,
        t2ⁿ     = T2,
        r1ⁿ     = 1.0 ./ T1,
        r2ⁿ     = 1.0 ./ T2,
        fⁿ      = f,
        Cⁿ      = k,
        wⁿppm   = dw,
        wⁿHz    = dw .* (B0 * γ * 1e-6),
        nSpecies = size(T1, 1) |> Int,
        maxRFsamples = maxRFsamples
    )
    return model
end
