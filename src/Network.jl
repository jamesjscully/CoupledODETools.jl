struct Network
    components::Vector{Component}
end
Network(varargs...) = Network([varargs...])

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
    scannedpars = [Set(scannedpars)...] #get rid of repetitions
    scannednames = [e[1] for e in scannedpars]
    pars = union(freepars, scannednames)
    uType = @SLVector Float32 Tuple([e[1] for e in eqtups])
    pType = @SLVector Float32 Tuple(pars)
    ics = uType([e[2] for e in eqtups])
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
    #make in place equations
    eqs! = []
    idxs = Dict()
    for i in eachindex(eqtups)
        eqex = eqtups[i][2]
        vs = [var for (var, eq) in eqtups]
        for j in eachindex(vs)
            idxs[vs[j]] = j
            eqex = flagreplace(eqtups[j][1], eqex, :(u[$j]))
        end
        ps = scannednames
        for j in eachindex(ps)
            eqex = flagreplace(ps[j], eqex, :(p[$j]))
        end
        push!(eqs!, Expr(:(=), :(du[$i]), eqex))
    end
    fscode = quote
        function fs(du, u, p, t)
            @inbounds $(Expr(:block, eqs!...))
            nothing
        end
    end
    # create search space
    space = Iterators.product([i[2].val for i in scannedpars]...) |> collect
    u0s = Float32[icarr...]
    return (
        f = fcode,
        fs = fscode,
        pType = pType,
        u0 = ics,
        u0s = u0s,
        space = space,
        idxs = idxs
    )
end
function scan(net::Network, tspan::Tuple{Float32,Float32};
    alg = Tsit5(), saveat = .1f0, output_func = (sol,i) -> (sol, false))
    n = net()
    u0 = n.u0s
    p = n.space[1]
    parr = [[n.space[i]...] for i in 1:length(n.space)]

    prob = ODEProblem(eval(n.fs), n.u0s, tspan, p)
    prob_func = (prob,i,repeat) -> remake(prob,p=parr[i]);
    monteprob = EnsembleProblem(prob, prob_func = prob_func, output_func = output_func)
    sol = solve(monteprob,alg,EnsembleGPUArray(),trajectories=length(parr),saveat=saveat)
    resol = reshape(sol.u, size(n.space)...)
end
