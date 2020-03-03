struct Network
    components::Vector{Component}
end
Network(varargs...) = Network([varargs...])

struct SharedPar{T}
    name::Symbol
    range::Vector{T}
end

PArray = Union{SharePar, Array}

function (net::Network)(;kwargs...)
    scannedpars = []
    freepars = []
    net = deepcopy(net)
    #deal with couplings first by inserting the expr with correct values
    for c ∈ net.components
        id = Symbol(:_, c.name)
        for coup ∈ c.couplings
            # first step is to replace the ~ denoted values with the source name
            coup.ex = postwalk(x -> @capture(x, ~ var_) ? Symbol(var, id) : x, coup.ex)
            # now update the parameters
            for p in c.pars
                if p.val isa Number
                    #for numerical pars replace sym with val
                    coup.ex = flagreplace(p.sym, coup.ex, p.val)
                else

                    coup.ex = flagreplace(p.sym, coup.ex, Symbol(p.sym,id))
                end
            end
        end
    end
    # insert exprs into couplings only ensures nested couplings are handled
    for sc ∈ net.components
        for scoup ∈ sc.couplings
            for tc ∈ net.components
                if tc.name ∈ scoup.targets
                    for tcoup ∈ tc.couplings
                        tcoup.ex = flaginsert(scoup.flag, tcoup.ex, scoup.ex)
                    end
                end
            end
        end
    end
    # remove flags from couplings to prevent contamination
    flags = [coup.flag for c ∈ net.components for coup ∈ c.couplings]
    for c ∈ net.components
        for coup ∈ c.couplings
            for flag ∈ flags
                if flagin(flag, coup.ex)
                    coup.ex = flagremove(flag, coup.ex)
                end
            end
        end
    end
    #insert couplings into state equations rhs
    for sc ∈ net.components
        for coup ∈ sc.couplings
            for tc ∈ net.components
                if tc.name ∈ coup.targets
                    for eqn ∈ tc.eqns
                        eqn.rhs = flaginsert(coup.flag, eqn.rhs, coup.ex)
                    end
                end
            end
        end
    end
    # remove flags from state eqns rhs
    for c ∈ net.components
        for eqn ∈ c.eqns
            for flag ∈ flags
                if flagin(flag, eqn.rhs)
                    eqn.rhs = flagremove(flag, eqn.rhs)
                end
            end
        end
    end
    #name variables in cells
    for c ∈ net.components
        id = Symbol(:_, c.name)
        for eqn ∈ c.eqns
            for vareqn ∈ c.eqns
                _f(x) = x == :($(vareqn.sym)) ? Symbol(vareqn.sym, id) : x
                eqn.rhs = postwalk( _f, eqn.rhs)
            end
        end
    end
    #handle parameters
    for c in net.components
        for eqn in c.eqns
            for p in c.pars
                newname = Symbol(p.sym, :_, c.name)
                if p.val isa Number
                    #for numerical pars replace sym with val
                    eqn.rhs = flagreplace(p.sym, eqn.rhs, p.val)
                else
                    eqn.rhs = flagreplace(p.sym, eqn.rhs, newname)
                end
                if p.val isa Array
                    push!(scannedpars, (newname,p))
                elseif p.val == :free
                    push!(freepars, newname)
                end
            end
        end
    end
    icarr = [eq.val for c in net.components for eq in c.eqns]
    eqtups = [(Symbol(eq.sym, :_, c.name), eq.rhs) for c in net.components for eq in c.eqns]
    #free parameters: for SLVector version treat scanned as free
    unique!(scannedpars) #get rid of repetitions
    scannednames = [e[1] for e in scannedpars]
    pars = union(freepars, scannednames)
    return (
        icarr = icarr,
        eqtups = eqtups,
        scannedpars = scannedpars,
        scannednames = scannednames,
        pars = pars)
end

function generate(n)
    uType = @SLVector Float64 Tuple([e[1] for e in n.eqtups])
    pType = @SLVector Float64 Tuple(n.pars)
    u0 = uType(n.icarr)
    # prepend pars with p. and vars with u. for out of place
    oeqarr = []
    for e in n.eqtups
        eqex = e[2]
        for v in keys(Dict(n.eqtups))
            eqex = flagreplace(v, eqex, :(u.$v))
        end
        for p in n.pars
            eqex = flagreplace(p, eqex, :(p.$p))
        end
        push!(oeqarr, Expr(:(=), e[1], eqex))
    end
    fcode_inner = Expr(:block, oeqarr..., :($uType($([k for (k,v) in n.eqtups]...))))
    f = :((u,p,t) -> $fcode_inner) |> rmlines
    return (f = f, u0 = u0, pType = pType, uType = uType)
end
function generate_ensemble(n)
    #for scanned parameters
    #make in place equations
    eqs! = []
    idxs = Dict()
    for i in eachindex(n.eqtups)
        eqex = n.eqtups[i][2]
        vs = [var for (var, eq) in n.eqtups]
        for j in eachindex(vs)
            idxs[vs[j]] = j
            eqex = flagreplace(n.eqtups[j][1], eqex, :(u[$j]))
        end
        #deal with shared pars
        shared = [(e.name, e.arr) for e in n.scannedpars if e isa SharedPar]
        j = 1 # insertion index
        # group shared into tuples with unique names
        shared = []
        for e in n.scannedpars
            if !(e isa SharedPar)
                eqex = flagreplace(e[1], eqex, :(p[$j]))
                j+=1
            else
                push!(shared, e)
            end
        end
        done = []
        for e in shared
            if !(e.name in done)
                for e2 in shared
                    if e2.name == e.name
                        eqex = flagreplace(e2[1], eqex, :(p[$j]))
                    end
                end
                push!(done, e.name)
                j+=1
            end
        end
        for e in n.scannedpars
            if e isa SharedPar
                eqex = flagreplace(e[1], eqex, :(p[$j]))
                j+=1
            end
        end

        unshared = [e for e in n.scannedpars if !(e isa SharedPar)]
        ps = n.scannednames
        for j in eachindex(ps)
            eqex = flagreplace(ps[j], eqex, :(p[$j]))
        end
        push!(eqs!, Expr(:(=), :(du[$i]), eqex))
    end
    f = quote
        (du, u, p, t) -> @inbounds $(Expr(:block, eqs!...))
    end |> rmlines
    # create search space
    space = Iterators.product([i[2].val for i in n.scannedpars]...) |> collect
    u0 = Float32[n.icarr...]
    return (f = f, u0 = u0, space = space, idxs = idxs)
end
function generate_cuarray(n)
    space = make_space(n.scannedpars)
    u0 = ArrayPartition(cu.([fill(u0s[i], size(spacecu.x[1])) for i = 1:length(n.eqtups)])...)
    cueqs = map(n.eqtups) do sym, eqn
        :(@__dot__ $(Expr(:(=), Symbol(:d, sym), eqn)))
    end
    f = quote
        (du, u, p, t) -> begin
            $(Expr(:tuple, [e[1] for e in eqtups]...)) = u.x
            $(Expr(:tuple, [Symbol(:d,e[1]) for e in eqtups]...)) = du.x
            $(Expr(:tuple, scannednames...)) = p.x
            @inbounds $(Expr(:block, cueqs...))
        end
    end |> rmlines |> code_to_f32 |> cucode
    return (f = f, u0 = u0, space = space)
end
