module CoupledODETools
    using OrdinaryDiffEq
    using DiffEqGPU
    using LabelledArrays
    using MLStyle: @λ
    using MacroTools: @capture

    export Component, @Component, Network, scan

    include("./utils.jl")
    include("./Component.jl")
    include("./Network.jl")


end # module
