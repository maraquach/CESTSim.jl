# A simple example of a continuous wave saturation scheme and a pseudo-readout event


## 1. Set up scanner and protocols
using CESTSim
import DifferentialEquations

sys = CESTSys(B₀ = 3.0, γ=267.5154109126009/2/π, γ̄=267.5154109126009, B₀MHz = 3.0 * 42.57639999999999)
model = load_model("examples/models/case3.yaml")

# * Alternatively, create a model struct with another function method: create_model(T1::Vector{…}, T2::Vector{…}, f::Vector{…}, k::Vector{…}, dw::Vector{…}, B0::Float64, maxRFsamples::Int64)
# model = create_model(["water", "glu"], [1.8; 1.0], [55e-3; 6.9e-3], [111; 36e-3]./111, [7450.0], [0.0; 3.2])

# CW parameters -- You can set this up however you'd like, this is just an example :)
p = (; 
    S = (; Δt = 1e-6, tp = 2.0, td = 0.0e-3, tspoil = 6.5e-3, np = 1), # CW Saturation 
    I = (; Δt = 1e-6, tp = 0e-3, te = 0.0e-3, tr = 0.0e-3, tspoil = 0, np = 0), # Imaging 
    D = (; Δt = 1e-6, spoil = 6.5e-3, td = 0.0e-3, tro = 0e-3, trec = 5) # Other delays
) # all units in seconds

rf = 2.0 + 0.0im

blocks = CoreBlocks( # This is good if we have the same types of RF for sat and imaging -- no fancy things
    Sat         = SequenceBlock(RF=rf,      T=p.S.Δt:p.S.Δt:p.S.tp    , Δω=nothing, ϕ=0.0    , EventType=:RF      , EventName="Pulse duration"), 
    SatTd       = SequenceBlock(RF=nothing, T=p.S.Δt:p.S.Δt:p.S.td    , Δω=nothing, ϕ=nothing, EventType=:Ignore  , EventName="Inter-pulse delay"),
    SatImgSpoil = SequenceBlock(RF=nothing, T=p.D.Δt:p.D.Δt:p.D.spoil , Δω=nothing, ϕ=nothing, EventType=:Spoiling, EventName="Pre-img spoiling"), 
    ImgTRO      = SequenceBlock(RF=nothing, T=p.D.Δt:p.D.Δt:p.D.spoil   , Δω=nothing, ϕ=nothing, EventType=:ADC     , EventName="Readout"),
    Recovery    = SequenceBlock(RF=nothing, T=p.D.Δt:p.D.Δt:p.D.trec  , Δω=nothing, ϕ=nothing, EventType=:Delay   , EventName="Recovery time")
)

## 2. Set up "sequence" 

w = [-300 ; LinRange(-15:0.1:15)] .* sys.B₀MHz

seq = CESTSequence(Name="Case 3 CW Seq", Author="Mara Quach", System=sys)

for Δω = w
    # satBlock = deepcopy(blocks.Sat) # copy instead of using the same allocation!
    blocks.Sat.Δω = Δω
    addBlock!(seq, blocks.Sat)
    addBlock!(seq, blocks.ImgTRO)
    addBlock!(seq, blocks.Recovery)
end

f = plot_CESTseq(seq, true)


## 3. Simulate
simParams = CESTSimParams(Δω = w, sys = sys, model = model, solver="Rodas5P(autodiff=AutoFiniteDiff())",reltol=1e-6,abstol=1e-6, interp=:Constant)

# If resetting M0 after each offset
r = simulate_CEST(seq, simParams,true)
# otherwise
# r = simulate_CEST(seq, simParams,false) 

fZ = CESTSim.plot_CESTresult(r.Z, simParams)
fZ[1]
