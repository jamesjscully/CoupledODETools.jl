abstract type Component end
using MLStyle: @λ

mutable struct Eqn
    sym
    rhs
    val
end
mutable struct Par
    sym
    val
end
mutable struct Coupling
    name
    ex
    flag
    targets
end

eqnhandle = @λ begin
    Expr(:tuple, Expr(:(=), sym, rhs), ic) -> Eqn(Symbol(String(sym)[2:end]),rhs,ic)
    Expr(:(=), sym, Expr(:tuple,rhs,ic)) -> Eqn(Symbol(String(sym)[2:end]),rhs,ic)
    _ -> nothing
end
couplinghandle = @λ begin
    Expr(:(=), name, Expr(:(->), flag, ex)) -> Coupling(name, ex.args[2], flag, Symbol[])
    _ -> nothing
end
parhandle = @λ begin
    Expr(:(=), sym, val) -> Par(sym, eval(val))
    _ -> nothing
end
macro Component(name, eqnex, parex)
    equations = filter(x -> x != nothing, eqnhandle.(eqnex.args))
    parameters = filter(x -> x != nothing, parhandle.(parex.args))
    couplings = filter(x -> x != nothing, couplinghandle.(eqnex.args))
    return quote
        struct $name <: Component
            name
            eqns::Vector{Eqn}
            pars::Vector{Par}
            couplings::Vector{Coupling}
        end
        function $(esc(name))(id; kwargs...)
            eqns = deepcopy($equations)
            pars = deepcopy($parameters)
            coups = deepcopy($couplings)
            # update default initial conditions
            for eqn in eqns
                if Symbol(eqn.sym,0) in keys(kwargs)
                    eqn.val = kwargs[Symbol(eqn.sym,0)]
                end
            end
            for par in pars
                if par.sym in keys(kwargs)
                    par.val = kwargs[par.sym]
                end
            end
            for coup in coups
                if coup.name in keys(kwargs)
                    val = kwargs[coup.name]
                    val isa Array ? coup.targets = val : coup.targets = [val]
                end
            end
            return $(esc(name))(id, eqns, pars, coups)
        end
    end
end
