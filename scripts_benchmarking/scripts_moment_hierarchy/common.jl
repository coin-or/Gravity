using JuMP
using PowerModels

pms_path = joinpath(dirname(pathof(PowerModels)), "..")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file", "-f"
            help = "the power system network data file"
            required = true
        "--time-limit", "-t"
            help = "puts a time limit on the sovler"
            arg_type = Float64
        "--bus-limit", "-b"
            help = "puts a time limit on the sovler"
            arg_type = Int64
            default = 1000
    end

    return parse_args(s)
end



""
function expression_branch_power_yt_from(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    if !haskey(var(pm, nw), :p)
        var(pm, nw)[:p] = Dict{Tuple{Int,Int,Int},Any}()
    end
    if !haskey(var(pm, nw), :q)
        var(pm, nw)[:q] = Dict{Tuple{Int,Int,Int},Any}()
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    expression_branch_power_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function expression_branch_power_yt_to(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    if !haskey(var(pm, nw), :p)
        var(pm, nw)[:p] = Dict{Tuple{Int,Int,Int},Any}()
    end
    if !haskey(var(pm, nw), :q)
        var(pm, nw)[:q] = Dict{Tuple{Int,Int,Int},Any}()
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    expression_branch_power_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end



""
function expression_branch_power_yt_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    var(pm, n, :p)[f_idx] = @NLexpression(pm.model,  (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    var(pm, n, :q)[f_idx] = @NLexpression(pm.model, -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end

""
function expression_branch_power_yt_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    var(pm, n, :p)[t_idx] = @NLexpression(pm.model,  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
    var(pm, n, :q)[t_idx] = @NLexpression(pm.model, -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
end



""
function expression_branch_power_yt_from(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    var(pm, n, :p)[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
    var(pm, n, :q)[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
end


""
function expression_branch_power_yt_to(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    var(pm, n, :p)[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
    var(pm, n, :q)[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
end




"opf with sparce branch power varibales"
function run_opf_sbpv(file, model_type::Type, solver; kwargs...)
    return run_model(file, model_type, solver, build_opf_sbpv; kwargs...)
end


""
function build_opf_sbpv(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)

    for i in ids(pm, :branch)
        expression_branch_power_yt_from(pm, i)
        expression_branch_power_yt_to(pm, i)
    end

    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for (i,branch) in ref(pm, :branch)
        constraint_voltage_angle_difference(pm, i)

        #constraint_thermal_limit_from(pm, i)
        #constraint_thermal_limit_to(pm, i)
        if haskey(branch, "rate_a")
            f_bus_id = branch["f_bus"]
            t_bus_id = branch["t_bus"]
            f_idx = (i, f_bus_id, t_bus_id)
            t_idx = (i, t_bus_id, f_bus_id)

            p_fr = @variable(pm.model, base_name="p_fr", start = 0.0)
            q_fr = @variable(pm.model, base_name="q_fr", start = 0.0)
            p_to = @variable(pm.model, base_name="p_to", start = 0.0)
            q_to = @variable(pm.model, base_name="q_to", start = 0.0)

            if typeof(pm) <: AbstractACPModel
                p_fr_exper = var(pm, :p, f_idx)
                @NLconstraint(pm.model, p_fr_exper == p_fr)
                q_fr_exper = var(pm, :q, f_idx)
                @NLconstraint(pm.model, q_fr_exper == q_fr)

                p_to_exper = var(pm, :p, t_idx)
                @NLconstraint(pm.model, p_to_exper == p_to)
                q_to_exper = var(pm, :q, t_idx)
                @NLconstraint(pm.model, q_to_exper == q_to)
            else
                @constraint(pm.model, var(pm, :p, f_idx) == p_fr)
                @constraint(pm.model, var(pm, :q, f_idx) == q_fr)

                @constraint(pm.model, var(pm, :p, t_idx) == p_to)
                @constraint(pm.model, var(pm, :q, t_idx) == q_to)
            end

            rating = branch["rate_a"]
            JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rating^2)
            JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rating^2)
        end
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

const branch_power_const_tol = 1e-3

function deactivate_rate_a!(network)
    network["active_rates"] = Int[]
    for (i,branch) in network["branch"]
        if haskey(branch, "rate_a")
            branch["rate_a_inactive"] = branch["rate_a"]
            delete!(branch, "rate_a")
        end
    end
end

function activate_rate_a!(network)
    if haskey(network, "active_rates")
        delete!(network, "active_rates")
    end

    for (i,branch) in network["branch"]
        if haskey(branch, "rate_a_inactive")
            branch["rate_a"] = branch["rate_a_inactive"]
            delete!(branch, "rate_a_inactive")
        end
    end
end











function calc_comp_lines_from_points(points)
    line_data = []
    for i in 3:2:length(points)
        x1 = points[i-2]
        y1 = points[i-1]
        x2 = points[i-0]
        y2 = points[i+1]

        m = (y2 - y1)/(x2 - x1)
        b = y1 - m * x1

        push!(line_data, (slope=m, intercept=b))
    end

    for i in 2:length(line_data)
        if line_data[i-1].slope > line_data[i].slope
            Memento.error(_LOGGER, "non-convex pwl function found in points $(points)\nlines: $(line_data)")
        end
    end

    return line_data
end

"computes the cost values that are relivent, based on min/max values"
function get_active_cost_points(component)
    pmin = component["pmin"]
    pmax = component["pmax"]
    ncost = component["ncost"]
    points = component["cost"]

    first_active = 1
    for i in 1:(ncost-1)
        x0 = points[2*i-1]
        x1 = points[2*(i+1)-1]
        if pmin >= x0
            first_active = i
        end
        if pmin <= x1
            break
        end
    end

    last_active = ncost
    for i in 1:(ncost-1)
        x0 = points[end - (2*(i+1)-1)]
        x1 = points[end - (2*i-1)]
        if pmax <= x1
            last_active = ncost - i + 1
        end
        if pmax >= x0
            break
        end
    end

    points = points[2*first_active - 1 : 2*last_active]
    ncost = div(length(points), 2)

    @assert points[1] <= pmin && points[1+2] >= pmin
    @assert points[end-3] <= pmax && points[end-1] >= pmax

    #=
    x1 = points[1]
    y1 = points[2]
    x2 = points[3]
    y2 = points[4]
    m = (y2 - y1)/(x2 - x1)
    if !isnan(m)
        cost_pmin = y2 - m*(x2 - pmin)
    else
        @assert isapprox(y1, y2)
        cost_pmin = y1
    end


    x1 = points[end-3]
    y1 = points[end-2]
    x2 = points[end-1]
    y2 = points[end]
    m = (y2 - y1)/(x2 - x1)
    if !isnan(m)
        cost_pmax = y2 - m*(x2 - pmax)
    else
        @assert isapprox(y1, y2)
        cost_pmax = y1
    end


    return ncost, points, cost_pmin, cost_pmax
    =#
    return ncost, points
end



""
function run_opf_pwl_lambda(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_lambda; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_lambda(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_lambda(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

""
function objective_min_fuel_and_flow_cost_lambda(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_lambda(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

""
function objective_min_fuel_and_flow_cost_pwl_lambda(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_lambda(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_lambda(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            ncost, points = get_active_cost_points(gen)

            pg_cost_lambda = JuMP.@variable(pm.model,
                [i in 1:ncost], base_name="$(n)_pg_cost_lambda",
                lower_bound = 0.0,
                upper_bound = 1.0
            )
            JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

            pg_expr = 0.0
            pg_cost_expr = 0.0
            for i in 1:ncost
                mw = points[2*i-1]
                cost = points[2*i]

                pg_expr += mw*pg_cost_lambda[i]
                pg_cost_expr += cost*pg_cost_lambda[i]
            end
            JuMP.@constraint(pm.model, pg_expr == sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)))
            pg_cost[i] = pg_cost_expr
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end


""
function run_opf_pwl_lambda_cl(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_lambda_cl; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_lambda_cl(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_lambda_cl(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

""
function objective_min_fuel_and_flow_cost_lambda_cl(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_lambda_cl(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

""
function objective_min_fuel_and_flow_cost_pwl_lambda_cl(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_lambda_cl(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_lambda_cl(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            ncost, points = get_active_cost_points(gen)

            if ncost == 2
                if isapprox(points[2], points[4]) # constant function
                    pg_cost[i] = points[2]
                else
                    line = calc_comp_lines_from_points(points)[1]
                    @assert !isnan(line.slope)
                    pg_cost[i] = line.slope*sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) + line.intercept
                end
                continue
            end

            pg_cost_lambda = JuMP.@variable(pm.model,
                [i in 1:ncost], base_name="$(n)_pg_cost_lambda",
                lower_bound = 0.0,
                upper_bound = 1.0
            )
            JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

            pg_expr = 0.0
            pg_cost_expr = 0.0
            for i in 1:ncost
                mw = points[2*i-1]
                cost = points[2*i]

                pg_expr += mw*pg_cost_lambda[i]
                pg_cost_expr += cost*pg_cost_lambda[i]
            end
            JuMP.@constraint(pm.model, pg_expr == sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)))
            pg_cost[i] = pg_cost_expr
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end


""
function run_opf_pwl_psi(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_psi; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_psi(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_psi(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

""
function objective_min_fuel_and_flow_cost_psi(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_psi(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

""
function objective_min_fuel_and_flow_cost_pwl_psi(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_psi(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_psi(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        gen_lines = Dict{Int64,Vector}()
        pg_cost_start = Dict{Int64,Float64}()
        pg_cost_min = Dict{Int64,Float64}()
        pg_cost_max = Dict{Int64,Float64}()

        for (i, gen) in nw_ref[:gen]
            @assert gen["model"] == 1
            ncost, points = get_active_cost_points(gen)

            mws = Float64[]
            costs = Float64[]
            for i in 1:ncost
                push!(mws, points[2*i-1])
                push!(costs, points[2*i])
            end

            gen_lines[i] = calc_comp_lines_from_points(points)
            pg_value = sum(JuMP.start_value(var(pm, n, :pg, i)[c]) for c in conductor_ids(pm, n))
            pg_cost_value = -Inf
            for line in gen_lines[i]
                pg_cost_value = max(pg_cost_value, line.slope*pg_value + line.intercept)
            end
            pg_cost_start[i] = pg_cost_value
            pg_cost_min[i] = costs[1]
            pg_cost_max[i] = costs[end]
        end

        #println(pg_cost_start)
        #println(pg_cost_min)
        #println(pg_cost_max)

        pg_cost = var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in ids(pm, n, :gen)], base_name="$(n)_pg_cost",
            lower_bound = pg_cost_min[i],
            upper_bound = pg_cost_max[i]
        )
        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            for line in gen_lines[i]
                JuMP.@constraint(pm.model, pg_cost[i] >= line.slope*sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) + line.intercept)
            end
        end
    end
end



""
function run_opf_pwl_psi_cl(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_psi_cl; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_psi_cl(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_psi_cl(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

""
function objective_min_fuel_and_flow_cost_psi_cl(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_psi_cl(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

""
function objective_min_fuel_and_flow_cost_pwl_psi_cl(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_psi_cl(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_psi_cl(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i, gen) in nw_ref[:gen]
            @assert gen["model"] == 1

            ncost, points = get_active_cost_points(gen)

            if ncost == 2
                if isapprox(points[2], points[4]) # constant function
                    pg_cost[i] = points[2]
                else
                    line = calc_comp_lines_from_points(points)[1]
                    @assert !isnan(line.slope)
                    pg_cost[i] = line.slope*sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) + line.intercept
                end
                continue
            end

            mws = Float64[]
            costs = Float64[]
            for i in 1:ncost
                push!(mws, points[2*i-1])
                push!(costs, points[2*i])
            end

            pg_cost[i] = JuMP.@variable(pm.model, base_name="$(n)_pg_cost_$(i)",
                lower_bound = costs[1],
                upper_bound = costs[end]
            )

            # gen pwl cost
            gen_lines = calc_comp_lines_from_points(points)
            for line in gen_lines
                JuMP.@constraint(pm.model, pg_cost[i] >= line.slope*sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) + line.intercept)
            end
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end





""
function run_opf_pwl_delta(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_delta; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_delta(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_delta(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end


""
function objective_min_fuel_and_flow_cost_delta(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_delta(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end


""
function objective_min_fuel_and_flow_cost_pwl_delta(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_delta(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost_delta(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            pmin = gen["pmin"]
            mws = Float64[]
            costs = Float64[]

            ncost, points = get_active_cost_points(gen)
            for i in 1:ncost
                push!(mws, points[2*i-1])
                push!(costs, points[2*i])
            end

            #=
            @assert mws[1] <= pmin && mws[2] >= pmin

            #adjust mws[1] to be pmin exactly
            if !isapprox(pmin, mws[1])
                x1 = mws[1]
                y1 = costs[1]
                x2 = mws[2]
                y2 = costs[2]

                m = (y2 - y1)/(x2 - x1)

                if !isnan(m)
                    cost = y2 - m*(x2 - pmin)

                    mws[1] = pmin
                    costs[1] = cost
                else
                    mws[1] = pmin
                end
            end

            # println(gen["pmin"])
            # println(costs)
            # println(mws)


            cost_offset = costs[1]
            mws_offset = mws[1]

            for i in 1:ncost
                mws[i] -= mws_offset
                costs[i] -= cost_offset
            end
            =#

            # println(gen["pmin"])
            # println(costs)
            # println(mws)

            cost_per_mw = Float64[0.0]
            for i in 2:ncost
                x0 = mws[i-1]
                y0 = costs[i-1]
                x1 = mws[i]
                y1 = costs[i]

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    push!(cost_per_mw, m)
                else
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                end
            end

            # println(cost_per_mw)
            # println(mws)
            # println()

            pg_cost_mw = JuMP.@variable(pm.model,
                [i in 2:ncost], base_name="$(n)_pg_cost_mw",
                lower_bound = 0.0,
                upper_bound = mws[i] - mws[i-1]
            )

            JuMP.@constraint(pm.model, mws[1] + sum(pg_cost_mw[i] for i in 2:ncost) == sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)))
            pg_cost[i] = costs[1] + sum(cost_per_mw[i]*pg_cost_mw[i] for i in 2:ncost)
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end



""
function run_opf_pwl_delta_cl(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_delta_cl; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_delta_cl(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_delta_cl(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end


""
function objective_min_fuel_and_flow_cost_delta_cl(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_delta_cl(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end


""
function objective_min_fuel_and_flow_cost_pwl_delta_cl(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_delta_cl(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost_delta_cl(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            pmin = gen["pmin"]
            mws = Float64[]
            costs = Float64[]

            ncost, points = get_active_cost_points(gen)
            for i in 1:ncost
                push!(mws, points[2*i-1])
                push!(costs, points[2*i])
            end

            if ncost == 2
                if isapprox(points[2], points[4]) # constant function
                    pg_cost[i] = points[2]
                else
                    line = calc_comp_lines_from_points(points)[1]
                    @assert !isnan(line.slope)
                    pg_cost[i] = line.slope*sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) + line.intercept
                end
                continue
            end

            #=
            @assert mws[1] <= pmin && mws[2] >= pmin

            #adjust mws[1] to be pmin exactly
            if !isapprox(pmin, mws[1])
                x1 = mws[1]
                y1 = costs[1]
                x2 = mws[2]
                y2 = costs[2]

                m = (y2 - y1)/(x2 - x1)

                if !isnan(m)
                    cost = y2 - m*(x2 - pmin)

                    mws[1] = pmin
                    costs[1] = cost
                else
                    mws[1] = pmin
                end
            end

            # println(gen["pmin"])
            # println(costs)
            # println(mws)


            cost_offset = costs[1]
            mws_offset = mws[1]

            for i in 1:ncost
                mws[i] -= mws_offset
                costs[i] -= cost_offset
            end
            =#

            # println(gen["pmin"])
            # println(costs)
            # println(mws)

            cost_per_mw = Float64[0.0]
            for i in 2:ncost
                x0 = mws[i-1]
                y0 = costs[i-1]
                x1 = mws[i]
                y1 = costs[i]

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    push!(cost_per_mw, m)
                else
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                end
            end

            # println(cost_per_mw)
            # println(mws)
            # println()

            pg_cost_mw = JuMP.@variable(pm.model,
                [i in 2:ncost], base_name="$(n)_pg_cost_mw",
                lower_bound = 0.0,
                upper_bound = mws[i] - mws[i-1]
            )

            JuMP.@constraint(pm.model, mws[1] + sum(pg_cost_mw[i] for i in 2:ncost) == sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)))
            pg_cost[i] = costs[1] + sum(cost_per_mw[i]*pg_cost_mw[i]   for i in 2:ncost)
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end



""
function run_opf_pwl_phi(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_phi; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_phi(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_phi(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

""
function objective_min_fuel_and_flow_cost_phi(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_phi(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end
end

""
function objective_min_fuel_and_flow_cost_pwl_phi(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_phi(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_phi(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            ncost, points = get_active_cost_points(gen)

            mws = Float64[]
            costs = Float64[]

            for j in 1:ncost
                push!(mws, points[2*j-1])
                push!(costs, points[2*j])
            end


            cost_per_mw = Float64[0.0]
            cost_per_mw_b = Float64[0.0]
            for j in 2:ncost
                x0 = mws[j-1]
                y0 = costs[j-1]
                x1 = mws[j]
                y1 = costs[j]

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    b = y1 - m * x1

                    push!(cost_per_mw, m)
                    push!(cost_per_mw_b, b)
                else
                    #println(y0, " ", y1)
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                    push!(cost_per_mw_b, y0)
                end
            end

            # println(cost_per_mw)
            # println(mws)
            # println()

            pmax = gen["pmax"]
            pg_phi = JuMP.@variable(pm.model,
                [j in 3:ncost], base_name="$(n)_pg_phi",
                lower_bound = 0.0,
                upper_bound = pmax - mws[j-1]
            )

            for j in 3:ncost
                JuMP.@constraint(pm.model, pg_phi[j] >= sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) - mws[j-1])
            end

            #pg_cost[i] = cost_per_mw[2] * (sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) - mws[1]) + costs[1]
            pg_cost[i] = cost_per_mw[2] * (sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n))) + cost_per_mw_b[2]
            if ncost > 2
                pg_cost[i] += sum((cost_per_mw[j] - cost_per_mw[j-1])*pg_phi[j] for j in 3:ncost)
            end
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end




""
function run_opf_pwl_phi_cl(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_opf_pwl_phi_cl; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_phi_cl(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_phi_cl(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

""
function objective_min_fuel_and_flow_cost_phi_cl(pm::AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_phi_cl(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end
end

""
function objective_min_fuel_and_flow_cost_pwl_phi_cl(pm::AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_phi_cl(pm; kwargs...)
    objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_phi_cl(pm::AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            ncost, points = get_active_cost_points(gen)

            mws = Float64[]
            costs = Float64[]

            for j in 1:ncost
                push!(mws, points[2*j-1])
                push!(costs, points[2*j])
            end

            if ncost == 2
                if isapprox(points[2], points[4]) # constant function
                    pg_cost[i] = points[2]
                else
                    line = calc_comp_lines_from_points(points)[1]
                    @assert !isnan(line.slope)
                    pg_cost[i] = line.slope*sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) + line.intercept
                end
                continue
            end

            cost_per_mw = Float64[0.0]
            cost_per_mw_b = Float64[0.0]
            for j in 2:ncost
                x0 = mws[j-1]
                y0 = costs[j-1]
                x1 = mws[j]
                y1 = costs[j]

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    b = y1 - m * x1

                    push!(cost_per_mw, m)
                    push!(cost_per_mw_b, b)
                else
                    #println(y0, " ", y1)
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                    push!(cost_per_mw_b, y0)
                end
            end

            # println(cost_per_mw)
            # println(mws)
            # println()

            pmax = gen["pmax"]
            pg_phi = JuMP.@variable(pm.model,
                [j in 3:ncost], base_name="$(n)_pg_phi",
                lower_bound = 0.0,
                upper_bound = pmax - mws[j-1]
            )

            for j in 3:ncost
                JuMP.@constraint(pm.model, pg_phi[j] >= sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) - mws[j-1])
            end

            #pg_cost[i] = cost_per_mw[2] * (sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n)) - mws[1]) + costs[1]
            pg_cost[i] = cost_per_mw[2] * (sum(var(pm, n, :pg, i)[c] for c in conductor_ids(pm, n))) + cost_per_mw_b[2]
            if ncost > 2
                pg_cost[i] += sum((cost_per_mw[j] - cost_per_mw[j-1])*pg_phi[j] for j in 3:ncost)
            end
        end

        report && PowerModels._IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end
