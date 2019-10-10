using CoupledODETools, DiffEqGPU, OrdinaryDiffEq, CuArrays

CuArrays.allowscalar(false)

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

@Component loren begin
    Dx = σ*(y-x), .1
    Dy = x*(ρ -g*k -z) - y, 0
    Dz = x*y - β*z, 0
    ρout = k -> ~x
end begin
    σ = range(9f0, 11f0, length = 10) |> collect
    β = range(2f0,4f0, length = 10) |> collect
    ρ = 28#range(20f0,35f0, length = 3) |> collect
    g = 0
end
l = loren(:l)
net1 = Network([l])
net1()

sol = scan(net1, (0f0, 100f0))

l1 = loren(:l1, ρout = :l2)
l2 = loren(:l2, ρout = :l1)
net = Network([l1, l2])
net()
