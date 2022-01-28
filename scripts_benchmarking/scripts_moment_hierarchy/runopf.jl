#!/usr/bin/env /Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia
using TSSOS
using ArgParse
using PowerModels

include("common.jl")
include("modelopf.jl") # include the file modelopf.jl
function main(parsed_args) 
case_file = parsed_args["file"]
data = parse_file(case_file)

#data = parse_file("pglib_opf_case" * case * ".m")

# the first order relaxation
model,_ = pop_opf_real(data, normal=true, AngleCons=true, LineLimit="relax")
n = model.n
m = model.m
numeq = model.numeq
supp = model.supp
coe = model.coe
mc = maximum(abs.(coe[1]))
coe[1]=coe[1]./mc

time = @elapsed begin
opt,sol,popd = cs_tssos_first(supp, coe, n, 1, numeq=numeq, CS=false, TS="MF", MomentOne=false)
end
opt *= mc
opt_level1 = opt
time_level1 = time
mb = maximum(maximum.(popd.sb)) # maximal block size
# gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $time, mb = $mb")

# the minimum order relaxation
model,_ = pop_opf_real(data, normal=true, AngleCons=true, LineLimit=true)
n = model.n
m = model.m
numeq = model.numeq
supp = model.supp
coe = model.coe
mc = maximum(abs.(coe[1]))
coe[1] = coe[1]./mc

time = @elapsed begin
opt,sol,popd = cs_tssos_first(supp, coe, n, "min", numeq=numeq, CS="MF", TS="block", MomentOne=true)
end
opt *= mc
opt_level2 = opt
time_level2 = time
maxc = maximum(popd.cliquesize) # maximal clique size
mb = maximum(maximum.(popd.sb)) # maximal block size
# gap = 100*(AC-opt)/AC # optimality gap
println("n = $n, m = $m")
println("mc = $maxc, opt = $opt, time = $time, mb = $mb")
println("$case_file , Solution Level1 = $opt_level1, time1 = $time_level1, Level2 = $opt_level2, time2 = $time_level2")
end

if isinteractive() == false
  main(parse_commandline())
end
