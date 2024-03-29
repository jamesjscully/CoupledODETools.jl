#define metaprogramming functions for inserting couplings
walk(x, inner, outer) = outer(x)
walk(x::Expr, inner, outer) = outer(Expr(x.head, map(inner, x.args)...))
postwalk(f, x) = walk(x, x -> postwalk(f, x), f)
function flaginsert(flag, targetex, insertex)
    _f = @λ begin
        Expr(:call, f, args...) -> begin
            if :($flag) in args
                if f in [:+,:*]
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

function make_space(scannedpars)
    vals = [e[2].val for e in scannedpars]
    size(vals)[1] == 1 ? vals = transpose(vals) : nothing
    arr = Expr[]
    for i = 1:length(vals)
        push!(arr,Expr(
            :comprehension,
                Expr(:generator,
                    Symbol(:v,i), [Expr(:(=), Symbol(:v,j), vals[j]) for j = 1:length(vals)]...
    ))) end
    ArrayPartition(cu.(eval.(arr))...)
end

code_to_f32 = @λ begin
    e::Expr ->  e.head == :ref ? e : Expr(e.head, map(code_to_f32, e.args)...)
    n::Number -> convert(Float32, n)
    a -> a
end

cucode = @λ begin
    e::Expr -> Expr(e.head, map(cucode, e.args)...)
    :exp -> :exp
    a -> a
end
rmlines = @λ begin
    e :: Expr           -> @match e.head begin
        :macrocall => Expr(e.head, map(rmlines, e.args)...)
        a          => Expr(e.head, filter(x -> x !== nothing, map(rmlines, e.args))...)
    end
      :: LineNumberNode -> nothing
    a                   -> a
end
