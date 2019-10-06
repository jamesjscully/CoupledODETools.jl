struct Network
    components::Vector{Component}
end

function (net::Network)(;kwargs...)
    # add host component name to coupling functions
    scannedpars = []
    freepars = []
    net = deepcopy(net)
    for c ∈ net.components
        id = Symbol(:_, c.name)
        for coup ∈ c.couplings
            coup.ex = postwalk(x -> @capture(x, ~ var_) ? Symbol(var, id) : x, coup.ex)
            for p in c.pars
                if p.val isa Number
                    #for numerical pars replace sym with val
                    coup.ex = flagreplace(p.sym, coup.ex, p.val)
                else
                    coup.ex = flagreplace(p.sym, coup.ex, Symbol(p.sym,id))
                end
                if p.val isa Array
                    push!(scannedpars, (Symbol(p.sym, id),p))
                elseif p.val == :free
                    push!(freepars, Symbol(p.sym, id))
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
    scannednames = [e[1] for e in scannedpars]
    pars = union(freepars, scannednames)
    uType = @SLVector Float32 Tuple(keys(Dict(eqtups)))
    pType = @SLVector Float32 Tuple(pars)
    ics = uType(icarr...)
    # prepend pars with p. and vars with u. for out of place
    oeqarr = []
    for e in eqtups
        eqex = e[2]
        for v in keys(Dict(eqtups))
            eqex = flagreplace(v, eqex, :(u.$v))
        end
        for p in pars
            eqex = flagreplace(p, eqex, :(p.$p))
        end
        push!(oeqarr, Expr(:(=), e[1], eqex))
    end
    fcode_inner = Expr(:block, oeqarr..., :($uType($(keys(Dict(eqtups))...))))
    fcode = Expr(:function, :(f(u,p,t)), fcode_inner)
    #for scanned parameters
    space = Iterators.product([e[2].val for e in scannedpars]...)
    #make in place equations
    eqs! = []
    for i in eachindex(eqtups)
        eqex = eqtups[i][2]
        for j in eachindex([var for (var, eq) in eqtups])
            eqex = flagreplace(eqtups[j][1], eqex, :(u[$j]))
        end
        for j in eachindex(scannednames)
            eqex = flagreplace(scannednames[j], eqex, :(p[$j]))
        end
        push!(eqs!, Expr(:(=), :(u[$i]), eqex))
    end
    fscancode = quote
        function fs(du, u, p, t)
            @inbounds $(Expr(:block, eqs!...))
            nothing
        end
    end
    return (f = fcode, fscan = fscancode, p = pType, u0 = ics)
end
