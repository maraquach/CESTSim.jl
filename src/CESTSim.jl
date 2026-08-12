## ^ This is a simple, (hopefully) non-convoluted Julia module to run CEST simulations. 

module CESTSim

using DifferentialEquations
using Parameters
using Plots
using LinearAlgebra
using LoopVectorization
using FastInterpolations
using YAML
using DSP

const γ = 42.577478461e6
const γ̄ = 42.577478461e6 * 2π
const γ̄M = 42.577478461 * 2π

include("ModelSupport.jl")
include("SeqSupport.jl")
include("Simulations.jl")
include("PulseSupport.jl")

export CESTModel, CESTSys, CESTSequence, CoreBlocks, SequenceBlock, CESTSimParams
export addBlock!, removeBlock!, insertBlock!, replaceBlock!
export create_M0, create_params, create_model, load_model
export plot_CESTtraj, plot_CESTseq
export solve_BME!, simulate_CEST
export blockADC, blockdummyRF 
export CESTSimTheme
export pulses

end
