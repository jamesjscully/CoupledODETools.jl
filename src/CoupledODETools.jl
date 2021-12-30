module CoupledODETools
    using LabelledArrays
    using MLStyle
    using MacroTools: @capture
    using RecursiveArrayTools
    using CUDA

    export Component, @Component, Network, SharedPar
    export generate, generate_cuarray, generate_ensemble

    include("./utils.jl")
    include("./Component.jl")
    include("./Network.jl")


end # module
