#define metaprogramming functions for inserting couplings
walk(x, inner, outer) = outer(x)
walk(x::Expr, inner, outer) = outer(Expr(x.head, map(inner, x.args)...))
postwalk(f, x) = walk(x, x -> postwalk(f, x), f)
function flaginsert(flag, targetex, insertex)
    _f = @λ begin
        Expr(:call, f, args...) -> begin
            if :($flag) in args
                if f in [:+,:*]
                    println(true)
                    Expr(:call, f, insertex, args...)
                elseif f == :/
                    args = [x == :($flag) ? :($insertex*$x) : x for x in args]
                    Expr(:call, f, args...)
                else
                    args = [x == :($flag) ? :($insertex+$x) : x for x in args]
                    Expr(:call, f, args...)
                end
            else
                Expr(:call, f, args...)
            end
        end
        a -> a
    end
    postwalk(_f,targetex)
end
function flagin(flag, ex)
    result = false
    postwalk(x -> x == flag ? result = true : nothing, ex)
    result
end
function flagremove(flag, targetex)
    _f = @λ begin
        Expr(:call, f, args...) -> begin
            args = filter(x -> x != :($flag), args)
            Expr(:call, f, args...)
        end
        a -> a
    end
    postwalk(_f,targetex)
end
function flagreplace(flag, targetex, insertex)
    _f = @λ begin
        Expr(:call, f, args...) -> begin
            if :($flag) in args
                args = [x == :($flag) ? insertex : x for x in args]
                Expr(:call, f, args...)
            else
                Expr(:call, f, args...)
            end
        end
        a -> a
    end
    postwalk(_f,targetex)
end
function scan(net::Network, tspan::Tuple{Float32,Float32};
    alg = Tsit5(), saveat = .1f0, output_func = (sol,i) -> (sol, false))
    n = net()
    u0 = n.u0s
    p = n.space[1]
    parr = [[n.space[i]...] for i in 1:length(n.space)]

    prob = ODEProblem(eval(n.fs), n.u0s, tspan, p)
    prob_func = (prob,i,repeat) -> remake(prob,p=parr[i]);
    monteprob = EnsembleProblem(prob, prob_func = prob_func)
    sol = solve(monteprob,alg,EnsembleGPUArray(),trajectories=length(parr),saveat=saveat)
    resol = reshape(sol.u, size(n.space)...)
end
