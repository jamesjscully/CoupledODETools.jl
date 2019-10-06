module CoupledODETools

using LabelledArrays
using MLStyle: @λ
using MacroTools: @capture
using IterTools: product

export Component, @Component, Network

include("./utils.jl")
include("./Component.jl")
include("./Network.jl")


end # module
