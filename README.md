# CoupledODEToolkit
```julia
@Component FitzughCell begin
    Dv = v - v^3 - w + Iext + gsyn, .0
    Dw = v + a -b*w, .1
    to_alpha = vpre -> ~v
    elec = gsyn -> gelec*(~v - v)
end begin
    Iext = .5
    a = :free
    b = .7
    gelec = .02
end

@Component AlphaSynapse begin
    DS = α*(1-S)/(1+exp(♉-vpre)/k) - β*S, 0
    synout = gsyn -> g*~S*(v-E)
end begin
    ♉ = -20
    k = 20
    E = -5
    g = .1
end

c1 = FitzughCell(:c1, elec = :c2)
c2 = FitzughCell(:c2, v0 = .1, w0 = .2, a = :free, elec = :c1)
s1 = AlphaSynapse(:s1)
net = Network([c1,c2])()
```
State variables are : D(sym) = (expr, ic)
couplings are : name = flag -> expr
The ~ Symbol is used to indicate that the variable belongs to the source, not target component.


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jamesjscully.github.io/CoupledODETools.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jamesjscully.github.io/CoupledODETools.jl/dev)
[![Build Status](https://travis-ci.com/jamesjscully/CoupledODETools.jl.svg?branch=master)](https://travis-ci.com/jamesjscully/CoupledODETools.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jamesjscully/CoupledODETools.jl?svg=true)](https://ci.appveyor.com/project/jamesjscully/CoupledODETools-jl)
[![Codecov](https://codecov.io/gh/jamesjscully/CoupledODETools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jamesjscully/CoupledODETools.jl)
