# CESTSim

Small Julia package for CEST/Bloch-McConnell simulations.

The package is in development and therefore not (yet) registered.

To install and use, type in your Julia REPL

```
]add git@github.com:maraquach/CESTSim.jl.git

using CESTSim
```

Two simple examples are provided in the examples directory. 

### Main features & relevant functions
* Import YAML model files or create your own struct (::CESTModel) - `load_model()` or `create_model()` 
* Create and modify Sequence objects/structs (::CESTSequence) from event blocks (::SequenceBlock) and core event definitions (::CoreBlocks) - `addBlock!()`, `removeBlock!()`, `insertBlock!()`, `replaceBlock!()`. Sequences can be constructed similarly (albeit not identically) to Pulseq.
* Simulate Z-spectra from ::CESTModel, ::CESTSys, ::CESTSequence, and ::CESTSimParams - `solve_BME!()`, `simulate_CEST()` 

More information can be found by looking at each individual function's documentation (e.g., `?addBlock!`).

### Simple workflow

Please see `examples/simple_CW.jl` for a sample workflow. 

‼️Important: In order for a Z-spectrum to be returned from `simulate_CEST()`, the ::CESTSequence argument must contain at least one ::SequenceBlock of EventType = :ADC

### Other helpful packages

It is also possible to load a Pulseq .seq file using **KomaMRI.jl**'s `read_seq()`. This can then be adapted to a ::CESTSequence and fed into `simulate_CEST()`

https://github.com/JuliaHealth/KomaMRI.jl

### Acknowledgment

- The ::CESTSequence structure is influenced by Pulseq (https://github.com/pulseq-admin/pulseq)
- Some WIP pulse constructions are adapted from MRiLab (https://github.com/leoliuf/MRiLab) and PyPulseq (https://github.com/imr-framework/pypulseq)
- Some results have been compared through the BMSim Challenge initiative (https://github.com/pulseq-cest/BMsim_challenge)

### Immediate next steps

* Pulse constructors
* Upload more examples + results
* More intuitive integration with KomaMRI

[![Build Status](https://github.com/maraquach/CESTNumSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maraquach/CESTNumSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/maraquach/CESTNumSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/maraquach/CESTNumSim.jl)
