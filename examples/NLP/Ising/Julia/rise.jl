#!/usr/bin/env julia

using ArgParse

using JuMP
using Ipopt

function main(args)
    println("regularizor: $(args["regularizor"])")
    println("loading: $(args["file"])")

    samples = readcsv(args["file"], Int)
    regularizor = args["regularizor"]

    (num_conf, num_row) = size(samples)
    num_spins = num_row - 1

    parameters = Array{Float64}(num_spins, num_spins)

    num_samples = 0
    for k=1:num_conf
      num_samples += samples[k,1]
    end

    println("spins: $(num_spins)")
    println("samples: $(num_samples)")
    println("num conf: $(num_conf)")
    lambda           = regularizor*sqrt(log((num_spins^2)/0.05)/num_samples)
    #RISEobjective(h) = exp(-h)
    #PLobjective(h)   = log(1 + exp(-2h))

    ipopt = IpoptSolver()

    tic()
    for u=1:num_spins
        current_row = u + 1
        nodal_stat  = [ samples[k,current_row] * (j == current_row ? 1 : samples[k,j]) for k=1:num_conf, j=2:num_row]

        m = Model(solver=ipopt)

        @variable(m, x[1:num_spins])
        @variable(m, z[1:num_spins] >= 0)
	@variable(m, obj)
	@NLobjective(m, Min,obj)
        @NLconstraint(m, obj>=sum((samples[k,1]/num_samples)*exp(-sum(x[i]*nodal_stat[k,i] for i=1:num_spins)) for k=1:num_conf) + 
            lambda*sum(z[j] for j=1:num_spins if u!=j)
        )
        # what about taking a log?

        #  Uncomment this line to forbid magnetic field reconstruction
        #  @constraint(m, x[u] == 0)
        
        for j in 1:num_spins
            @constraint(m, z[j] >=  x[j]) #z_plus
            @constraint(m, z[j] >= -x[j]) #z_minus
        end

        # setvalue(x,[0 for j=1:num_spins])
        status = solve(m)

        println(u, " = ", getvalue(x))
        parameters[u,1:num_spins] = deepcopy(getvalue(x))
    end
    time = toc()

    println("DATA, $(args["file"]), $(num_spins), $(num_samples), $(args["regularizor"]), $(time)")

    if args["output"] != nothing
        writecsv(args["output"], parameters)
    else
        println("result:")
        println(parameters)
    end
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--file", "-f"
            help = "the sample input data file (.csv)"
            required = true
        "--regularizor", "-r"
            help = "objective regularizor"
            arg_type = Float64
            default = 0.2
        "--output", "-o"
            help = "the result data ouput file (.csv)"
    end

    return parse_args(s)
end

if isinteractive() == false
    main(parse_commandline())
end
