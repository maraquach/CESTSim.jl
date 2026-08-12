using CESTSim
using Test

@testset "CESTSim.jl" begin
    p = CESTSimParams(sys=CESTSim.defaultSystem(false),
                    model=create_model("2PoolsAmide7T"),
                    Δω=[0.0])
    r = simulate_CEST(CESTSim.dummySeq, p)
    
    testrf = pulses.rect(5.0, 0.0, 1e-3, 1e-6)
end
