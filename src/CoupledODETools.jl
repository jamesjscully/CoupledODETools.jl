module CoupledODETools
    using LabelledArrays
    using MLStyle: @λ
    using MacroTools: @capture
    using RecursiveArrayTools, CuArrays

    export Component, @Component, Network

    include("./utils.jl")
    include("./Component.jl")
    include("./Network.jl")


end # module
