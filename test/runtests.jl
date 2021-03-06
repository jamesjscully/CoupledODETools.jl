using CoupledODETools
using Test

import CoupledODETools
import CoupledODETools: flaginsert
import CoupledODETools: flagreplace

@testset "CoupledODETools.jl" begin
    #test metaprogramming functions
    @test flaginsert(:f, :(f+a+t), :t) == :(t+f+a+t)
    @test flaginsert(:f, :(f-a-t), :t) == :(t+f-a-t)
    @test flaginsert(:f, :(√(f-a*t)), :q) == :(√(q+f-a*t))
    @test flaginsert(:f, :(a*f*b), :t) == :(t*a*f*b)
    @test flaginsert(:f, :(x^f), :a) == :(x^(a+f))
    @test flaginsert(:f, :(t/f), :a) == :(t/(a*f))
    @test flaginsert(:f, :(t/(1+f)), :a) == :(t/(a+1+f))

end
flagreplace(:s, :((f+a+t)), :a) == :(f+a+t)
flagreplace(:a, :((f+a+t)), :s) == :(f+s+t)


@testset "free parameters -> SLArray" begin

    @Component FitzughCell begin
        Dv = v - v^3 - w + Iext + gsyn, .0
        Dw = v + a -b*w, .1
        to_alpha = vpre -> ~v
        elec = gsyn -> gelec*(~v - v)
    end begin
        Iext = .5
        a = :free
        b = .7
        gelec = 3
    end
    @Component AlphaSynapse begin
        DS = α*(1-S)/(1+exp(♉-vpre)*k) - β*S, 0
        synout = gsyn -> g*~S*(v-E)
    end begin
        α = :free
        β = .2
        ♉ = -20
        k = 20
        E = -5
        g = .1
    end
    #test defaults
    c = FitzughCell(:c)
    @test c.couplings[2].targets == []
    @test c.pars[1].val == .5
    @test c.eqns[1].sym == :v
    @test c.eqns[1].val == .0
    #test coupling
    c1 = FitzughCell(:c1, elec = :c2)
    c2 = FitzughCell(:c2, v0 = .1, w0 = .2, a = 3, elec = :c1, to_alpha = :s1)
    s1 = AlphaSynapse(:s1, synout = [:c2])

    @test c1.eqns[1].val == .0
    @test c2.eqns[1].val == .1
    @test c2.pars[2].val == 3

    net = Network([c1,c2,s1])()
end
