using PowerModels
# using SparseArrays

abstract type AbstractSparsePolyModel end

mutable struct SparsePolyModel <: AbstractSparsePolyModel
    n
    m
    numeq
    nbus
    ng
    nb
    supp
    coe
    dg
end

function normalize(coe)
    mc = maximum(abs.(coe))
    return coe./mc
end

function move_zero!(supp,coe)
    ind=[abs(coe[i])>=1e-8 for i=1:length(coe)]
    return supp[ind],coe[ind]
end

function fl_sum(vector)
    return mapreduce(x->x, +, vector, init = 0.0)
end

function bfind(A, l, a)
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if ndims(A)==2
            temp=A[:, mid]
        else
            temp=A[mid]
        end
        if temp==a
           return mid
        elseif temp<a
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

# function move_zero!(col,row,nz,coe)
#     i=1
#     while i<=length(coe)
#         if abs(coe[i])<=1e-8
#             deleteat!(coe,i)
#             deleteat!(row,col[i]:(col[i+1]-1))
#             deleteat!(nz,col[i]:(col[i+1]-1))
#             lrow=col[i+1]-col[i]
#             deleteat!(col,i)
#             col[i:end].-=lrow
#         else
#             i+=1
#         end
#     end
#     col=convert(Vector{UInt32},col)
#     return col,row,nz,coe
# end

# Voltage only complex formulization
function pop_opf_com(data; normal=true, AngleCons=false, LineLimit=false)
    silence()
    # path = joinpath(dirname(dirname(pathof(PolyOPF))), "pglib", "pglib_opf_" * case * ".m")
    # data = parse_file(path)
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0]
    nbus=length(ref[:bus])
    nb=length(ref[:branch])
    ng=length(ref[:gen])
    n=nbus
    m=4*nbus+2*ng+length(keys(ref[:ref_buses]))
    if AngleCons==true
        m+=2*nb
    end
    if LineLimit==true||LineLimit=="relax"
        m+=2*nb
    end
    numeq=2*nbus-2*ng+length(keys(ref[:ref_buses]))
    dg=2*ones(Int, m)
    supp=Vector{Vector{Vector{Vector{UInt16}}}}(undef, m+1)
    coe=Vector{Vector{complex(Float64)}}(undef, m+1)

    gens=collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    coe[1]=[sum(gen["cost"][3] for (i,gen) in ref[:gen])]
    supp[1]=[[[], []]]
    k=2

    bus=collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    for i=1:nbus
        supp[k]=[[[], []], [[i], [i]]]
        coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1]
        supp[k+1]=[[[], []], [[i], [i]]]
        coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1]
        k+=2
    end

    if AngleCons==true||LineLimit==true||LineLimit=="relax"
        for (i, branch) in ref[:branch]
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            vr = bfind(bus,nbus,branch["f_bus"])
            vt = bfind(bus,nbus,branch["t_bus"])
            srt=sort([vr;vt])
            a1=g+g_fr
            b1=-(b+b_fr)
            c1=-g*tr+b*ti
            d1=b*tr+g*ti
            a2=(g+g_to)*tm^2
            b2=-(b+b_to)*tm^2
            c2=-(g*tr+b*ti)
            d2=-(-b*tr+g*ti)

            # angle differences
            if AngleCons==true
                coe[k]=[tan(branch["angmax"])+im;tan(branch["angmax"])-im]
                supp[k]=[[[vr], [vt]], [[vt], [vr]]]
                coe[k+1]=[-tan(branch["angmin"])-im;-tan(branch["angmin"])+im]
                supp[k+1]=[[[vr], [vt]], [[vt], [vr]]]
                if normal==true
                    coe[k]=normalize(coe[k])
                    coe[k+1]=normalize(coe[k+1])
                end
                k+=2
            end

            # thermal limits
            if LineLimit==true
                coe[k]=[branch["rate_a"]^2*tm^4;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
                supp[k]=[[[], []], [[vr;vr], [vr;vr]], [[vr;vr], srt], [srt, [vr;vr]], [srt, srt]]
                supp[k],coe[k]=move_zero!(supp[k],coe[k])
                coe[k+1]=[branch["rate_a"]^2*tm^4;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
                supp[k+1]=[[[], []], [[vt;vt], [vt;vt]], [srt, [vt;vt]], [[vt;vt], srt], [srt, srt]]
                supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
                dg[k-1:k]=[4;4]
                if normal==true
                    coe[k]=normalize(coe[k])
                    coe[k+1]=normalize(coe[k+1])
                end
                k+=2
            elseif LineLimit=="relax"
                mvr=ref[:bus][bus[vr]]["vmin"]^2
                mvt=ref[:bus][bus[vt]]["vmin"]^2
                coe[k]=[branch["rate_a"]^2*tm^4/mvr;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
                supp[k]=[[[], []], [[vr], [vr]], [[vr], [vt]], [[vt], [vr]], [[vt], [vt]]]
                supp[k],coe[k]=move_zero!(supp[k],coe[k])
                coe[k+1]=[branch["rate_a"]^2*tm^4/mvt;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
                supp[k+1]=[[[], []], [[vt], [vt]], [[vr], [vt]], [[vt], [vr]], [[vr], [vr]]]
                supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
                if normal==true
                    coe[k]=normalize(coe[k])
                    coe[k+1]=normalize(coe[k+1])
                end
                k+=2
            end
        end
    end

    # active/reactive power
    k1=k
    k=k1+4*ng
    for (r, i) in enumerate(bus)
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        ploads=fl_sum(load["pd"] for load in bus_loads)
        qloads=fl_sum(load["qd"] for load in bus_loads)
        sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
        if length(ref[:bus_gens][i])>1
            @error "There are more than one generator at one bus!"
        elseif isempty(ref[:bus_gens][i])
            coe[k]=zeros(complex(Float64), 2*length(ref[:bus_arcs][i])+2)
            supp[k]=Vector{Vector{Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
            supp[k][1:2]=[[UInt16[], UInt16[]], [UInt16[r], UInt16[r]]]
            coe[k][1:2]=[ploads;sgs]
            coe[k+1]=zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
            supp[k+1]=Vector{Vector{Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
            supp[k+1][1:2]=[[UInt16[], UInt16[]], [UInt16[r], UInt16[r]]]
            coe[k+1][1:2]=[qloads;-sbs]
        else
            gen_id=ref[:bus_gens][i][1]
            gen=ref[:gen][gen_id]
            gc2=gen["cost"][2]
            supp[k1:k1+3]=[[[UInt16[], UInt16[]], [UInt16[r], UInt16[r]]] for l=1:4]
            coe[k1]=[gen["pmax"]-ploads, -sgs]
            coe[k1+1]=[-gen["pmin"]+ploads, sgs]
            coe[k1+2]=[gen["qmax"]-qloads, sbs]
            coe[k1+3]=[-gen["qmin"]+qloads, -sbs]
            push!(supp[1], [UInt16[r], UInt16[r]])
            ngen=length(supp[1])
            coe[1][1]+=gc2*ploads
            push!(coe[1], gc2*sgs)
        end
        j=1
        for flow in ref[:bus_arcs][i]
            branch=ref[:branch][flow[1]]
            vr = bfind(bus,nbus,branch["f_bus"])
            vt = bfind(bus,nbus,branch["t_bus"])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            a1=(g+g_fr)/tm^2
            b1=-(b+b_fr)/tm^2
            c1=(-g*tr+b*ti)/tm^2
            d1=(b*tr+g*ti)/tm^2
            a2=g+g_to
            b2=-(b+b_to)
            c2=-(g*tr+b*ti)/tm^2
            d2=-(-b*tr+g*ti)/tm^2

            if isempty(ref[:bus_gens][i])
                supp[k][j+2:j+3]=[[[vr], [vt]], [[vt], [vr]]]
                supp[k+1][j+2:j+3]=[[[vr], [vt]], [[vt], [vr]]]
            else
                push!(supp[1], [[vr], [vt]], [[vt], [vr]])
                push!(supp[k1], [[vr], [vt]], [[vt], [vr]])
                push!(supp[k1+1], [[vr], [vt]], [[vt], [vr]])
                push!(supp[k1+2], [[vr], [vt]], [[vt], [vr]])
                push!(supp[k1+3], [[vr], [vt]], [[vt], [vr]])
            end
            if vr==r
                if isempty(ref[:bus_gens][i])
                    coe[k][2]+=a1
                    coe[k][j+2:j+3]=[(c1+d1*im)/2, (c1-d1*im)/2]
                    coe[k+1][2]+=b1
                    coe[k+1][j+2:j+3]=[(-c1*im+d1)/2, (c1*im+d1)/2]
                else
                    coe[1][ngen]+=gc2*a1
                    push!(coe[1], gc2*(c1+d1*im)/2, gc2*(c1-d1*im)/2)
                    coe[k1][2]-=a1
                    push!(coe[k1], -(c1+d1*im)/2, -(c1-d1*im)/2)
                    coe[k1+1][2]+=a1
                    push!(coe[k1+1], (c1+d1*im)/2, (c1-d1*im)/2)
                    coe[k1+2][2]-=b1
                    push!(coe[k1+2], -(-c1*im+d1)/2, -(c1*im+d1)/2)
                    coe[k1+3][2]+=b1
                    push!(coe[k1+3], (-c1*im+d1)/2, (c1*im+d1)/2)
                end
            else
                if isempty(ref[:bus_gens][i])
                    coe[k][2]+=a2
                    coe[k][j+2:j+3]=[(c2-d2*im)/2, (c2+d2*im)/2]
                    coe[k+1][2]+=b2
                    coe[k+1][j+2:j+3]=[(c2*im+d2)/2, (-c2*im+d2)/2]
                else
                    coe[1][ngen]+=gc2*a2
                    push!(coe[1], gc2*(c2-d2*im)/2, gc2*(c2+d2*im)/2)
                    coe[k1][2]-=a2
                    push!(coe[k1], -(c2-d2*im)/2, -(c2+d2*im)/2)
                    coe[k1+1][2]+=a2
                    push!(coe[k1+1], (c2-d2*im)/2, (c2+d2*im)/2)
                    coe[k1+2][2]-=b2
                    push!(coe[k1+2], -(c2*im+d2)/2, -(-c2*im+d2)/2)
                    coe[k1+3][2]+=b2
                    push!(coe[k1+3], (c2*im+d2)/2, (-c2*im+d2)/2)
                end
            end
            j+=2
        end
        if isempty(ref[:bus_gens][i])
            supp[k],coe[k]=move_zero!(supp[k],coe[k])
            supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
            if normal==true
                coe[k]=normalize(coe[k])
                coe[k+1]=normalize(coe[k+1])
            end
            k+=2
        else
            supp[k1],coe[k1]=move_zero!(supp[k1],coe[k1])
            supp[k1+1],coe[k1+1]=move_zero!(supp[k1+1],coe[k1+1])
            supp[k1+2],coe[k1+2]=move_zero!(supp[k1+2],coe[k1+2])
            supp[k1+3],coe[k1+3]=move_zero!(supp[k1+3],coe[k1+3])
            if normal==true
                coe[k1]=normalize(coe[k1])
                coe[k1+1]=normalize(coe[k1+1])
                coe[k1+2]=normalize(coe[k1+2])
                coe[k1+3]=normalize(coe[k1+3])
            end
            k1+=4
            if gen["cost"][1]>0
                lsupp=length(supp[1])
                gc1=gen["cost"][1]
                coe[1][1]+=gc1*ploads^2
                for l=ngen:lsupp
                    push!(supp[1], [[supp[1][l][1];supp[1][l][1]], [supp[1][l][2];supp[1][l][2]]])
                    push!(coe[1], gc1*coe[1][l]^2/gc2^2)
                    for p=l+1:lsupp
                        push!(supp[1], [sort([supp[1][l][1];supp[1][p][1]]), sort([supp[1][l][2];supp[1][p][2]])])
                        push!(coe[1], 2*gc1*coe[1][l]*coe[1][p]/gc2^2)
                    end
                    coe[1][l]+=2*gc1*ploads*coe[1][l]/gc2
                end
            end
        end
    end
    supp[1],coe[1]=move_zero!(supp[1],coe[1])
    supp[1],coe[1]=resort(supp[1],coe[1],field="complex")

    # reference voltage
    for key in keys(ref[:ref_buses])
        i=bfind(bus,nbus,key)
        supp[k]=[[[i;i], []], [[i], [i]], [[], [i;i]]]
        coe[k]=[1;-2;1]
        k+=1
    end
    return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg)
end

# Voltage only real formulization
# function pop_opf_real(data; normal=true, AngleCons=false, LineLimit=false)
#     silence()
#     PowerModels.standardize_cost_terms!(data, order=2)
#     ref = PowerModels.build_ref(data)[:nw][0]
#     nbus=length(ref[:bus])
#     nb=length(ref[:branch])
#     ng=length(ref[:gen])
#     n=2*nbus
#     m=4*nbus+2*ng+length(keys(ref[:ref_buses]))
#     if AngleCons==true
#         m+=2*nb
#     end
#     if LineLimit==true||LineLimit=="relax"
#         m+=2*nb
#     end
#     numeq=2*nbus-2*ng+length(keys(ref[:ref_buses]))
#     dg=2*ones(Int, m)
#     startpoint=zeros(n)
#     supp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
#     coe=Vector{Vector{Float64}}(undef, m+1)
#
#     gens=collect(keys(ref[:gen]))
#     sort!(gens)
#     # objective function
#     coe[1]=[sum(gen["cost"][3] for (i,gen) in ref[:gen])]
#     supp[1]=[[]]
#     k=2
#
#     bus=collect(keys(ref[:bus]))
#     sort!(bus)
#     # voltage magnitude constraints
#     for i=1:nbus
#         supp[k]=[[], [i;i], [i+nbus;i+nbus]]
#         coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
#         supp[k+1]=[[], [i;i], [i+nbus;i+nbus]]
#         coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
#         startpoint[i]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
#         startpoint[i+nbus]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
#         k+=2
#     end
#
#     if AngleCons==true||LineLimit==true||LineLimit=="relax"
#         for (i, branch) in ref[:branch]
#             g, b = PowerModels.calc_branch_y(branch)
#             tr, ti = PowerModels.calc_branch_t(branch)
#             g_fr = branch["g_fr"]
#             b_fr = branch["b_fr"]
#             g_to = branch["g_to"]
#             b_to = branch["b_to"]
#             tm = branch["tap"]
#             vr = bfind(bus,nbus,branch["f_bus"])
#             vt = bfind(bus,nbus,branch["t_bus"])
#             srt=sort([vr;vt])
#             dsrt=sort([vr;vr;vt;vt])
#             ab1=(g+g_fr)^2+(b+b_fr)^2
#             cd1=(-g*tr+b*ti)^2+(b*tr+g*ti)^2
#             acbd1=(g+g_fr)*(-g*tr+b*ti)-(b+b_fr)*(b*tr+g*ti)
#             bcad1=-(b+b_fr)*(-g*tr+b*ti)-(g+g_fr)*(b*tr+g*ti)
#             ab2=(g+g_to)^2*tm^4+(b+b_to)^2*tm^4
#             cd2=(g*tr+b*ti)^2+(-b*tr+g*ti)^2
#             acbd2=-(g+g_to)*tm^2*(g*tr+b*ti)+(b+b_to)*tm^2*(-b*tr+g*ti)
#             bcad2=(b+b_to)*tm^2*(g*tr+b*ti)+(g+g_to)*tm^2*(-b*tr+g*ti)
#
#             # angle differences
#             if AngleCons==true
#                 coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
#                 supp[k]=[srt, srt.+nbus, [vt;vr+nbus], [vr;vt+nbus]]
#                 coe[k+1]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
#                 supp[k+1]=[[vt;vr+nbus], [vr;vt+nbus], srt, srt.+nbus]
#                 if normal==true
#                     coe[k]=normalize(coe[k])
#                     coe[k+1]=normalize(coe[k+1])
#                 end
#                 k+=2
#             end
#
#             # thermal limits
#             if LineLimit==true
#                 coe[k]=[branch["rate_a"]^2*tm^4;-ab1;-2*ab1;-ab1;-cd1;-cd1;-cd1;-cd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*bcad1;2*bcad1;-2*bcad1;2*bcad1]
#                 supp[k]=[[], [vr;vr;vr;vr], [vr;vr;vr+nbus;vr+nbus], [vr+nbus;vr+nbus;vr+nbus;vr+nbus], [vt;vt;vr+nbus;vr+nbus], [vr;vr;vt+nbus;vt+nbus],
#                 dsrt, dsrt.+nbus, sort([vr;vr;vr;vt]), [vr;vr;sort([vr+nbus;vt+nbus])], [srt;vr+nbus;vr+nbus], sort([vr;vr;vr;vt]).+nbus, [sort([vr;vr;vt]);vr+nbus],
#                 [vr;vr;vr;vt+nbus], [vt;vr+nbus;vr+nbus;vr+nbus], [vr;sort([vr+nbus;vr+nbus;vt+nbus])]]
#                 supp[k],coe[k]=move_zero!(supp[k],coe[k])
#                 coe[k+1]=[branch["rate_a"]^2*tm^4;-ab2;-2*ab2;-ab2;-cd2;-cd2;-cd2;-cd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*bcad2;2*bcad2;-2*bcad2;2*bcad2]
#                 supp[k+1]=[[], [vt;vt;vt;vt], [vt;vt;vt+nbus;vt+nbus], [vt+nbus;vt+nbus;vt+nbus;vt+nbus], [vr;vr;vt+nbus;vt+nbus], [vt;vt;vr+nbus;vr+nbus],
#                 dsrt, dsrt.+nbus, sort([vt;vt;vt;vr]), [vt;vt;sort([vt+nbus;vr+nbus])], [srt;vt+nbus;vt+nbus], sort([vt;vt;vt;vr]).+nbus, [sort([vt;vt;vr]);vt+nbus],
#                 [vt;vt;vt;vr+nbus], [vr;vt+nbus;vt+nbus;vt+nbus], [vt;sort([vt+nbus;vt+nbus;vr+nbus])]]
#                 supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
#                 dg[k-1:k]=[4;4]
#                 if normal==true
#                     coe[k]=normalize(coe[k])
#                     coe[k+1]=normalize(coe[k+1])
#                 end
#                 k+=2
#             elseif LineLimit=="relax"
#                 mvr=ref[:bus][bus[vr]]["vmin"]^2
#                 mvt=ref[:bus][bus[vt]]["vmin"]^2
#                 coe[k]=[branch["rate_a"]^2*tm^4/mvr;-ab1;-ab1;-cd1;-cd1;-2*acbd1;2*bcad1;-2*bcad1;-2*acbd1]
#                 supp[k]=[[], [vr;vr], [vr+nbus;vr+nbus], [vt;vt], [vt+nbus;vt+nbus], srt, [vr;vt+nbus], [vt;vr+nbus], srt.+nbus]
#                 supp[k],coe[k]=move_zero!(supp[k],coe[k])
#                 coe[k+1]=[branch["rate_a"]^2*tm^4/mvt;-ab2;-ab2;-cd2;-cd2;-2*acbd2;2*bcad2;-2*bcad2;-2*acbd2]
#                 supp[k+1]=[[], [vt;vt], [vt+nbus;vt+nbus], [vr;vr], [vr+nbus;vr+nbus], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
#                 supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
#                 if normal==true
#                     coe[k]=normalize(coe[k])
#                     coe[k+1]=normalize(coe[k+1])
#                 end
#                 k+=2
#             end
#         end
#     end
#
#     # active/reactive power
#     k1=k
#     k=k1+4*ng
#     for (r, i) in enumerate(bus)
#         bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
#         bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
#         ploads=fl_sum(load["pd"] for load in bus_loads)
#         qloads=fl_sum(load["qd"] for load in bus_loads)
#         sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
#         sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
#         if length(ref[:bus_gens][i])>1
#             @error "There are more than one generator at one bus!"
#         elseif isempty(ref[:bus_gens][i])
#             coe[k]=zeros(Float64, 4*length(ref[:bus_arcs][i])+3)
#             supp[k]=Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
#             supp[k][1:3]=[[], [r;r], [r+nbus;r+nbus]]
#             coe[k][1:3]=[ploads, sgs, sgs]
#             coe[k+1]=zeros(Float64, 4*length(ref[:bus_arcs][i])+3)
#             supp[k+1]=Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
#             supp[k+1][1:3]=[[], [r;r], [r+nbus;r+nbus]]
#             coe[k+1][1:3]=[qloads, -sbs, -sbs]
#         else
#             gen_id=ref[:bus_gens][i][1]
#             gen=ref[:gen][gen_id]
#             gc2=gen["cost"][2]
#             supp[k1:k1+3]=[[[], [r;r], [r+nbus;r+nbus]] for l=1:4]
#             coe[k1]=[gen["pmax"]-ploads, -sgs, -sgs]
#             coe[k1+1]=[-gen["pmin"]+ploads, sgs, sgs]
#             coe[k1+2]=[gen["qmax"]-qloads, sbs, sbs]
#             coe[k1+3]=[-gen["qmin"]+qloads, -sbs, -sbs]
#             push!(supp[1], [r;r], [r+nbus;r+nbus])
#             ngen=length(supp[1])
#             coe[1][1]+=gc2*ploads
#             push!(coe[1], gc2*sgs, gc2*sgs)
#         end
#         j=1
#         for flow in ref[:bus_arcs][i]
#             branch=ref[:branch][flow[1]]
#             vr = bfind(bus,nbus,branch["f_bus"])
#             vt = bfind(bus,nbus,branch["t_bus"])
#             srt=sort([vr;vt])
#             g, b = PowerModels.calc_branch_y(branch)
#             tr, ti = PowerModels.calc_branch_t(branch)
#             g_fr = branch["g_fr"]
#             b_fr = branch["b_fr"]
#             g_to = branch["g_to"]
#             b_to = branch["b_to"]
#             tm = branch["tap"]
#             a1=(g+g_fr)/tm^2
#             b1=-(b+b_fr)/tm^2
#             c1=(-g*tr+b*ti)/tm^2
#             d1=(b*tr+g*ti)/tm^2
#             a2=g+g_to
#             b2=-(b+b_to)
#             c2=-(g*tr+b*ti)/tm^2
#             d2=-(-b*tr+g*ti)/tm^2
#
#             if isempty(ref[:bus_gens][i])
#                 supp[k][j+3:j+6]=[srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
#                 supp[k+1][j+3:j+6]=[srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
#             else
#                 push!(supp[1], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus)
#                 push!(supp[k1], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus)
#                 push!(supp[k1+1], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus)
#                 push!(supp[k1+2], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus)
#                 push!(supp[k1+3], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus)
#             end
#             if vr==r
#                 if isempty(ref[:bus_gens][i])
#                     coe[k][2:3].+=a1
#                     coe[k][j+3:j+6]=[c1;-d1;d1;c1]
#                     coe[k+1][2:3].+=b1
#                     coe[k+1][j+3:j+6]=[d1;c1;-c1;d1]
#                 else
#                     coe[1][ngen-1:ngen].+=gc2*a1
#                     push!(coe[1], gc2*c1, -gc2*d1, gc2*d1, gc2*c1)
#                     coe[k1][2:3].-=a1
#                     push!(coe[k1], -c1, d1, -d1, -c1)
#                     coe[k1+1][2:3].+=a1
#                     push!(coe[k1+1], c1, -d1, d1, c1)
#                     coe[k1+2][2:3].-=b1
#                     push!(coe[k1+2], -d1, -c1, c1, -d1)
#                     coe[k1+3][2:3].+=b1
#                     push!(coe[k1+3], d1, c1, -c1, d1)
#                 end
#             else
#                 if isempty(ref[:bus_gens][i])
#                     coe[k][2:3].+=a2
#                     coe[k][j+3:j+6]=[c2;d2;-d2;c2]
#                     coe[k+1][2:3].+=b2
#                     coe[k+1][j+3:j+6]=[d2;-c2;c2;d2]
#                 else
#                     coe[1][ngen-1:ngen].+=gc2*a2
#                     push!(coe[1], gc2*c2, gc2*d2, -gc2*d2, gc2*c2)
#                     coe[k1][2:3].-=a2
#                     push!(coe[k1], -c2, -d2, d2, -c2)
#                     coe[k1+1][2:3].+=a2
#                     push!(coe[k1+1], c2, d2, -d2, c2)
#                     coe[k1+2][2:3].-=b2
#                     push!(coe[k1+2], -d2, c2, -c2, -d2)
#                     coe[k1+3][2:3].+=b2
#                     push!(coe[k1+3], d2, -c2, c2, d2)
#                 end
#             end
#             j+=4
#         end
#         if isempty(ref[:bus_gens][i])
#             supp[k],coe[k]=move_zero!(supp[k],coe[k])
#             supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
#             if normal==true
#                 coe[k]=normalize(coe[k])
#                 coe[k+1]=normalize(coe[k+1])
#             end
#             k+=2
#         else
#             supp[k1],coe[k1]=move_zero!(supp[k1],coe[k1])
#             supp[k1+1],coe[k1+1]=move_zero!(supp[k1+1],coe[k1+1])
#             supp[k1+2],coe[k1+2]=move_zero!(supp[k1+2],coe[k1+2])
#             supp[k1+3],coe[k1+3]=move_zero!(supp[k1+3],coe[k1+3])
#             if normal==true
#                 coe[k1]=normalize(coe[k1])
#                 coe[k1+1]=normalize(coe[k1+1])
#                 coe[k1+2]=normalize(coe[k1+2])
#                 coe[k1+3]=normalize(coe[k1+3])
#             end
#             k1+=4
#             if gen["cost"][1]>0
#                 lsupp=length(supp[1])
#                 gc1=gen["cost"][1]
#                 coe[1][1]+=gc1*ploads^2
#                 for l=ngen-1:lsupp
#                     push!(supp[1], sadd(supp[1][l], supp[1][l]))
#                     push!(coe[1], gc1*coe[1][l]^2/gc2^2)
#                     for p=l+1:lsupp
#                         push!(supp[1], sadd(supp[1][l], supp[1][p]))
#                         push!(coe[1], 2*gc1*coe[1][l]*coe[1][p]/gc2^2)
#                     end
#                     coe[1][l]+=2*gc1*ploads*coe[1][l]/gc2
#                 end
#             end
#         end
#     end
#     supp[1],coe[1]=move_zero!(supp[1],coe[1])
#     supp[1],coe[1]=resort(supp[1],coe[1])
#
#     # reference voltage
#     for key in keys(ref[:ref_buses])
#         i=bfind(bus,nbus,key)
#         supp[k]=[[i+nbus;i+nbus]]
#         coe[k]=[1]
#         k+=1
#     end
#     return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg),startpoint
# end

function sadd(a, b)
    c=[a;b]
    return sort!(c)
end

function resort(supp, coe; field="real")
    nsupp=copy(supp)
    sort!(nsupp)
    unique!(nsupp)
    l=length(nsupp)
    if field=="real"
        ncoe=zeros(l)
    else
        ncoe=zeros(ComplexF64, l)
    end
    for i=1:length(supp)
        locb=bfind(nsupp, l, supp[i])
        ncoe[locb]+=coe[i]
    end
    return nsupp,ncoe
end

# function pop_opf_com(data::Dict{String, Any}; normal=true, AngleCons=false, LineLimit=false)
#     silence()
#     PowerModels.standardize_cost_terms!(data, order=2)
#     ref = PowerModels.build_ref(data)[:nw][0]
#     nbus=length(ref[:bus])
#     nb=length(ref[:branch])
#     ng=length(ref[:gen])
#     n=nbus+ng
#     m=4*nbus+2*ng+length(keys(ref[:ref_buses]))
#     if AngleCons==true
#         m+=2*nb
#     end
#     if LineLimit==true||LineLimit=="relax"
#         m+=2*nb
#     end
#     numeq=2*nbus+length(keys(ref[:ref_buses]))
#     dg=2*ones(Int, m)
#     supp=Vector{Vector{Vector{Vector{UInt16}}}}(undef, m+1)
#     coe=Vector{Vector{complex(Float64)}}(undef, m+1)
#
#     gens=collect(keys(ref[:gen]))
#     sort!(gens)
#     # objective function
#     nc=5*ng+1
#     coe[1]=Vector{complex(Float64)}(undef, nc)
#     supp[1]=Vector{Vector{Vector{UInt16}}}(undef, nc)
#     coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
#     supp[1][1]=[[], []]
#     for i=1:ng
#         gen=ref[:gen][gens[i]]
#         coe[1][5*(i-1)+2:5*(i-1)+3]=[0.5*gen["cost"][2], 0.5*gen["cost"][2]]
#         coe[1][5*(i-1)+4:5*(i-1)+6]=[0.25*gen["cost"][1], 0.5*gen["cost"][1], 0.25*gen["cost"][1]]
#         supp[1][5*(i-1)+2:5*(i-1)+6]=[[[nbus+i], []], [[], [nbus+i]], [[nbus+i;nbus+i], []], [[nbus+i], [nbus+i]], [[], [nbus+i;nbus+i]]]
#     end
#     supp[1],coe[1]=move_zero!(supp[1],coe[1])
#     k=2
#
#     bus=collect(keys(ref[:bus]))
#     sort!(bus)
#     # voltage magnitude constraints
#     for i=1:nbus
#         supp[k]=[[[], []], [[i], [i]]]
#         coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1]
#         supp[k+1]=[[[], []], [[i], [i]]]
#         coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1]
#         k+=2
#     end
#
#     if AngleCons==true||LineLimit==true||LineLimit=="relax"
#         for (i, branch) in ref[:branch]
#             g, b = PowerModels.calc_branch_y(branch)
#             tr, ti = PowerModels.calc_branch_t(branch)
#             g_fr = branch["g_fr"]
#             b_fr = branch["b_fr"]
#             g_to = branch["g_to"]
#             b_to = branch["b_to"]
#             tm = branch["tap"]
#             vr = bfind(bus,nbus,branch["f_bus"])
#             vt = bfind(bus,nbus,branch["t_bus"])
#             srt=sort([vr, vt])
#             a1=g+g_fr
#             b1=-(b+b_fr)
#             c1=-g*tr+b*ti
#             d1=b*tr+g*ti
#             a2=(g+g_to)*tm^2
#             b2=-(b+b_to)*tm^2
#             c2=-(g*tr+b*ti)
#             d2=-(-b*tr+g*ti)
#
#             # angle differences
#             if AngleCons==true
#                 coe[k]=[tan(branch["angmax"])+im;tan(branch["angmax"])-im]
#                 supp[k]=[[[vr], [vt]], [[vt], [vr]]]
#                 coe[k+1]=[-tan(branch["angmin"])-im;-tan(branch["angmin"])+im]
#                 supp[k+1]=[[[vr], [vt]], [[vt], [vr]]]
#                 if normal==true
#                     coe[k]=normalize(coe[k])
#                     coe[k+1]=normalize(coe[k+1])
#                 end
#                 k+=2
#             end
#
#             # thermal limits
#             if LineLimit==true
#                 coe[k]=[branch["rate_a"]^2*tm^4;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
#                 supp[k]=[[[], []], [[vr;vr], [vr;vr]], [[vr;vr], srt], [srt, [vr;vr]], [srt, srt]]
#                 supp[k],coe[k]=move_zero!(supp[k],coe[k])
#                 coe[k+1]=[branch["rate_a"]^2*tm^4;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
#                 supp[k+1]=[[[], []], [[vt;vt], [vt;vt]], [srt, [vt;vt]], [[vt;vt], srt], [srt, srt]]
#                 supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
#                 dg[k-1:k]=[4;4]
#                 if normal==true
#                     coe[k]=normalize(coe[k])
#                     coe[k+1]=normalize(coe[k+1])
#                 end
#                 k+=2
#             elseif LineLimit=="relax"
#                 mvr=ref[:bus][bus[vr]]["vmin"]^2
#                 mvt=ref[:bus][bus[vt]]["vmin"]^2
#                 coe[k]=[branch["rate_a"]^2*tm^4/mvr;-(a1^2+b1^2);-(a1*c1+b1*d1)+(b1*c1-a1*d1)*im;-(a1*c1+b1*d1)+(a1*d1-b1*c1)*im;-(c1^2+d1^2)]
#                 supp[k]=[[[], []], [[vr], [vr]], [[vr], [vt]], [[vt], [vr]], [[vt], [vt]]]
#                 supp[k],coe[k]=move_zero!(supp[k],coe[k])
#                 coe[k+1]=[branch["rate_a"]^2*tm^4/mvt;-(a2^2+b2^2);-(a2*c2+b2*d2)+(a2*d2-b2*c2)*im;-(a2*c2+b2*d2)+(b2*c2-a2*d2)*im;-(c2^2+d2^2)]
#                 supp[k+1]=[[[], []], [[vt], [vt]], [[vr], [vt]], [[vt], [vr]], [[vr], [vr]]]
#                 supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
#                 if normal==true
#                     coe[k]=normalize(coe[k])
#                     coe[k+1]=normalize(coe[k+1])
#                 end
#                 k+=2
#             end
#         end
#     end
#
#     # power generation bound
#     zero_pgen=UInt16[]
#     for i=1:ng
#         gen=ref[:gen][gens[i]]
#         if gen["pmax"]>=1e-6
#             coe[k]=[-4*gen["pmin"]*gen["pmax"];2*gen["pmin"]+2*gen["pmax"];2*gen["pmin"]+2*gen["pmax"];-1;-2;-1]
#             if normal==true
#                 coe[k]=normalize(coe[k])
#             end
#             supp[k]=[[[], []], [[i+nbus], []], [[], [i+nbus]], [[i+nbus;i+nbus], []], [[i+nbus], [i+nbus]], [[], [i+nbus;i+nbus]]]
#             supp[k],coe[k]=move_zero!(supp[k],coe[k])
#             k+=1
#         else
#             push!(zero_pgen, i)
#             numeq+=1
#         end
#         coe[k]=[-4*gen["qmin"]*gen["qmax"];-2*gen["qmin"]*im-2*gen["qmax"]*im;2*gen["qmin"]*im+2*gen["qmax"]*im;1;-2;1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         supp[k]=[[[], []], [[i+nbus], []], [[], [i+nbus]], [[i+nbus;i+nbus], []], [[i+nbus], [i+nbus]], [[], [i+nbus;i+nbus]]]
#         supp[k],coe[k]=move_zero!(supp[k],coe[k])
#         k+=1
#     end
#
#     # active/reactive power
#     for (r, i) in enumerate(bus)
#         bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
#         bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
#         coe[k]=zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
#         supp[k]=Vector{Vector{Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
#         supp[k][1:2]=[[[], []], [[r], [r]]]
#         coe[k][1]=fl_sum(load["pd"] for load in bus_loads)
#         coe[k][2]=fl_sum(shunt["gs"] for shunt in bus_shunts)
#         coe[k+1]=zeros(ComplexF64, 2*length(ref[:bus_arcs][i])+2)
#         supp[k+1]=Vector{Vector{Vector{UInt16}}}(undef, 2*length(ref[:bus_arcs][i])+2)
#         supp[k+1][1:2]=[[[], []], [[r], [r]]]
#         coe[k+1][1]=fl_sum(load["qd"] for load in bus_loads)
#         coe[k+1][2]=-fl_sum(shunt["bs"] for shunt in bus_shunts)
#         j=1
#         for flow in ref[:bus_arcs][i]
#             branch=ref[:branch][flow[1]]
#             vr = bfind(bus,nbus,branch["f_bus"])
#             vt = bfind(bus,nbus,branch["t_bus"])
#             g, b = PowerModels.calc_branch_y(branch)
#             tr, ti = PowerModels.calc_branch_t(branch)
#             g_fr = branch["g_fr"]
#             b_fr = branch["b_fr"]
#             g_to = branch["g_to"]
#             b_to = branch["b_to"]
#             tm = branch["tap"]
#             a1=(g+g_fr)/tm^2
#             b1=-(b+b_fr)/tm^2
#             c1=(-g*tr+b*ti)/tm^2
#             d1=(b*tr+g*ti)/tm^2
#             a2=g+g_to
#             b2=-(b+b_to)
#             c2=-(g*tr+b*ti)/tm^2
#             d2=-(-b*tr+g*ti)/tm^2
#
#             supp[k][j+2:j+3]=[[[vr], [vt]], [[vt], [vr]]]
#             supp[k+1][j+2:j+3]=[[[vr], [vt]], [[vt], [vr]]]
#             if vr==r
#                 coe[k][2]+=a1
#                 coe[k][j+2:j+3]=[(c1+d1*im)/2, (c1-d1*im)/2]
#                 coe[k+1][2]+=b1
#                 coe[k+1][j+2:j+3]=[(-c1*im+d1)/2, (c1*im+d1)/2]
#             else
#                 coe[k][2]+=a2
#                 coe[k][j+2:j+3]=[(c2-d2*im)/2, (c2+d2*im)/2]
#                 coe[k+1][2]+=b2
#                 coe[k+1][j+2:j+3]=[(c2*im+d2)/2, (-c2*im+d2)/2]
#             end
#             j+=2
#         end
#         if !isempty(ref[:bus_gens][i])
#             for gen_id in ref[:bus_gens][i]
#                 gen=bfind(gens, ng, gen_id)
#                 push!(supp[k], [[gen+nbus], []], [[], [gen+nbus]])
#                 push!(coe[k], -0.5, -0.5)
#                 push!(supp[k+1], [[gen+nbus], []], [[], [gen+nbus]])
#                 push!(coe[k+1], 0.5*im, -0.5*im)
#             end
#         end
#         if normal==true
#             coe[k]=normalize(coe[k])
#             coe[k+1]=normalize(coe[k+1])
#         end
#         supp[k],coe[k]=move_zero!(supp[k],coe[k])
#         supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
#         k+=2
#     end
#
#     # reference voltage
#     for key in keys(ref[:ref_buses])
#         i=bfind(bus,nbus,key)
#         supp[k]=[[[i;i], []], [[i], [i]], [[], [i;i]]]
#         coe[k]=[1;-2;1]
#         k+=1
#     end
#
#     # zero power generation
#     for i in zero_pgen
#         supp[k]=[[[i+nbus], []], [[], [i+nbus]]]
#         coe[k]=[1;1]
#         dg[k-1]=1
#         k+=1
#     end
#
#     return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg)
# end

function pop_opf_real(data::Dict{String, Any}; normal=true, AngleCons=false, LineLimit=false)
    silence()
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0]
    nbus=length(ref[:bus])
    nb=length(ref[:branch])
    ng=length(ref[:gen])
    n=2*nbus+2*ng
    m=4*nbus+2*ng+length(keys(ref[:ref_buses]))
    if AngleCons==true
        m+=2*nb
    end
    if LineLimit==true||LineLimit=="relax"
        m+=2*nb
    end
    numeq=2*nbus+length(keys(ref[:ref_buses]))
    dg=2*ones(Int, m)
    startpoint=zeros(n)
    supp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
    coe=Vector{Vector{Float64}}(undef, m+1)

    gens=collect(keys(ref[:gen]))
    sort!(gens)
    # objective function
    nc=2*ng+1
    coe[1]=Vector{Float64}(undef, nc)
    supp[1]=Vector{Vector{UInt16}}(undef, nc)
    coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
    supp[1][1]=[]
    for i=1:ng
        gen=ref[:gen][gens[i]]
        coe[1][2*(i-1)+2:2*(i-1)+3]=[gen["cost"][2];gen["cost"][1]]
        supp[1][2*(i-1)+2:2*(i-1)+3]=[[2*nbus+i], [2*nbus+i;2*nbus+i]]
    end
    supp[1],coe[1]=move_zero!(supp[1],coe[1])
    k=2

    bus=collect(keys(ref[:bus]))
    sort!(bus)
    # voltage magnitude constraints
    for i=1:nbus
        supp[k]=[[], [i;i], [i+nbus;i+nbus]]
        coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
        supp[k+1]=[[], [i;i], [i+nbus;i+nbus]]
        coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
        startpoint[i]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        startpoint[i+nbus]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        k+=2
    end

    if AngleCons==true||LineLimit==true||LineLimit=="relax"
        for (i, branch) in ref[:branch]
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            vr = bfind(bus,nbus,branch["f_bus"])
            vt = bfind(bus,nbus,branch["t_bus"])
            srt=sort([vr;vt])
            dsrt=sort([vr;vr;vt;vt])
            ab1=(g+g_fr)^2+(b+b_fr)^2
            cd1=(-g*tr+b*ti)^2+(b*tr+g*ti)^2
            acbd1=(g+g_fr)*(-g*tr+b*ti)-(b+b_fr)*(b*tr+g*ti)
            bcad1=-(b+b_fr)*(-g*tr+b*ti)-(g+g_fr)*(b*tr+g*ti)
            ab2=(g+g_to)^2*tm^4+(b+b_to)^2*tm^4
            cd2=(g*tr+b*ti)^2+(-b*tr+g*ti)^2
            acbd2=-(g+g_to)*tm^2*(g*tr+b*ti)+(b+b_to)*tm^2*(-b*tr+g*ti)
            bcad2=(b+b_to)*tm^2*(g*tr+b*ti)+(g+g_to)*tm^2*(-b*tr+g*ti)

            # angle differences
            if AngleCons==true
                coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
                supp[k]=[srt, srt.+nbus, [vt;vr+nbus], [vr;vt+nbus]]
                coe[k+1]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
                supp[k+1]=[[vt;vr+nbus], [vr;vt+nbus], srt, srt.+nbus]
                if normal==true
                    coe[k]=normalize(coe[k])
                    coe[k+1]=normalize(coe[k+1])
                end
                k+=2
            end

            # thermal limits
            if LineLimit==true
                coe[k]=[branch["rate_a"]^2*tm^4;-ab1;-2*ab1;-ab1;-cd1;-cd1;-cd1;-cd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*bcad1;2*bcad1;-2*bcad1;2*bcad1]
                supp[k]=[[], [vr;vr;vr;vr], [vr;vr;vr+nbus;vr+nbus], [vr+nbus;vr+nbus;vr+nbus;vr+nbus], [vt;vt;vr+nbus;vr+nbus], [vr;vr;vt+nbus;vt+nbus],
                dsrt, dsrt.+nbus, sort([vr;vr;vr;vt]), [vr;vr;sort([vr+nbus;vt+nbus])], [srt;vr+nbus;vr+nbus], sort([vr;vr;vr;vt]).+nbus, [sort([vr;vr;vt]);vr+nbus],
                [vr;vr;vr;vt+nbus], [vt;vr+nbus;vr+nbus;vr+nbus], [vr;sort([vr+nbus;vr+nbus;vt+nbus])]]
                supp[k],coe[k]=move_zero!(supp[k],coe[k])
                coe[k+1]=[branch["rate_a"]^2*tm^4;-ab2;-2*ab2;-ab2;-cd2;-cd2;-cd2;-cd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*bcad2;2*bcad2;-2*bcad2;2*bcad2]
                supp[k+1]=[[], [vt;vt;vt;vt], [vt;vt;vt+nbus;vt+nbus], [vt+nbus;vt+nbus;vt+nbus;vt+nbus], [vr;vr;vt+nbus;vt+nbus], [vt;vt;vr+nbus;vr+nbus],
                dsrt, dsrt.+nbus, sort([vt;vt;vt;vr]), [vt;vt;sort([vt+nbus;vr+nbus])], [srt;vt+nbus;vt+nbus], sort([vt;vt;vt;vr]).+nbus, [sort([vt;vt;vr]);vt+nbus],
                [vt;vt;vt;vr+nbus], [vr;vt+nbus;vt+nbus;vt+nbus], [vt;sort([vt+nbus;vt+nbus;vr+nbus])]]
                supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
                dg[k-1:k]=[4;4]
                if normal==true
                    coe[k]=normalize(coe[k])
                    coe[k+1]=normalize(coe[k+1])
                end
                k+=2
            elseif LineLimit=="relax"
                mvr=ref[:bus][bus[vr]]["vmin"]^2
                mvt=ref[:bus][bus[vt]]["vmin"]^2
                coe[k]=[branch["rate_a"]^2*tm^4/mvr;-ab1;-ab1;-cd1;-cd1;-2*acbd1;2*bcad1;-2*bcad1;-2*acbd1]
                supp[k]=[[], [vr;vr], [vr+nbus;vr+nbus], [vt;vt], [vt+nbus;vt+nbus], srt, [vr;vt+nbus], [vt;vr+nbus], srt.+nbus]
                supp[k],coe[k]=move_zero!(supp[k],coe[k])
                coe[k+1]=[branch["rate_a"]^2*tm^4/mvt;-ab2;-ab2;-cd2;-cd2;-2*acbd2;2*bcad2;-2*bcad2;-2*acbd2]
                supp[k+1]=[[], [vt;vt], [vt+nbus;vt+nbus], [vr;vr], [vr+nbus;vr+nbus], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
                supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
                if normal==true
                    coe[k]=normalize(coe[k])
                    coe[k+1]=normalize(coe[k+1])
                end
                k+=2
            end
        end
    end

    # power generation bound
    zero_pgen=UInt16[]
    for i=1:ng
        gen=ref[:gen][gens[i]]
        if gen["pmax"]>=1e-6
            coe[k]=[-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
            if normal==true
                coe[k]=normalize(coe[k])
            end
            supp[k]=[[], [i+2*nbus], [i+2*nbus;i+2*nbus]]
            supp[k],coe[k]=move_zero!(supp[k],coe[k])
            startpoint[i+2*nbus]=(gen["pmin"]+gen["pmax"])/2
            k+=1
        else
            push!(zero_pgen, i)
            startpoint[i+2*nbus]=0
            numeq+=1
        end
        coe[k]=[-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
        if normal==true
            coe[k]=normalize(coe[k])
        end
        supp[k]=[[], [i+2*nbus+ng], [i+2*nbus+ng;i+2*nbus+ng]]
        supp[k],coe[k]=move_zero!(supp[k],coe[k])
        startpoint[i+2*nbus+ng]=0
        k+=1
    end

    # active/reactive power
    for (r, i) in enumerate(bus)
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        coe[k]=zeros(4*length(ref[:bus_arcs][i])+3)
        supp[k]=Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
        supp[k][1:3]=[[], [r;r], [r+nbus;r+nbus]]
        coe[k][1]=fl_sum(load["pd"] for load in bus_loads)
        coe[k+1]=zeros(4*length(ref[:bus_arcs][i])+3)
        supp[k+1]=Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
        supp[k+1][1:3]=[[], [r;r], [r+nbus;r+nbus]]
        coe[k+1][1]=fl_sum(load["qd"] for load in bus_loads)
        sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
        coe[k][2:3]=[sgs;sgs]
        coe[k+1][2:3]=[-sbs;-sbs]
        j=1
        for flow in ref[:bus_arcs][i]
            branch=ref[:branch][flow[1]]
            vr = bfind(bus,nbus,branch["f_bus"])
            vt = bfind(bus,nbus,branch["t_bus"])
            srt=sort([vr;vt])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            a1=(g+g_fr)/tm^2
            b1=-(b+b_fr)/tm^2
            c1=(-g*tr+b*ti)/tm^2
            d1=(b*tr+g*ti)/tm^2
            a2=g+g_to
            b2=-(b+b_to)
            c2=-(g*tr+b*ti)/tm^2
            d2=-(-b*tr+g*ti)/tm^2

            supp[k][j+3:j+6]=[srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
            supp[k+1][j+3:j+6]=[srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
            if vr==r
                coe[k][2:3].+=a1
                coe[k][j+3:j+6]=[c1;-d1;d1;c1]
                coe[k+1][2:3].+=b1
                coe[k+1][j+3:j+6]=[d1;c1;-c1;d1]
            else
                coe[k][2:3].+=a2
                coe[k][j+3:j+6]=[c2;d2;-d2;c2]
                coe[k+1][2:3].+=b2
                coe[k+1][j+3:j+6]=[d2;-c2;c2;d2]
            end
            j+=4
        end
        if !isempty(ref[:bus_gens][i])
            for gen_id in ref[:bus_gens][i]
                gen=bfind(gens, ng, gen_id)
                push!(supp[k], [gen+2*nbus])
                push!(coe[k], -1)
                push!(supp[k+1], [gen+2*nbus+ng])
                push!(coe[k+1], -1)
            end
        end
        if normal==true
            coe[k]=normalize(coe[k])
            coe[k+1]=normalize(coe[k+1])
        end
        supp[k],coe[k]=move_zero!(supp[k],coe[k])
        supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
        k+=2
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        i=bfind(bus,nbus,key)
        supp[k]=[[i+nbus;i+nbus]]
        coe[k]=[1]
        k+=1
    end

    # zero power generation
    for i in zero_pgen
        supp[k]=[[i+2*nbus]]
        coe[k]=[1]
        dg[k-1]=1
        k+=1
    end

    return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg),startpoint
end

# function clique_opf_two(n,m,nbus,supp;alg="MF",minimize=true)
#     G=SimpleGraph(n)
#     for i=1:m+1, j = 1:supp[i].n
#         add_clique!(G,supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)])
#     end
#     for i=1:nbus
#         add_edge!(G,i,i+nbus)
#     end
#     if alg=="NC"
#         cliques,cql,cliquesize=max_cliques(G)
#     else
#         cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
#     end
#     uc=unique(cliquesize)
#     sizes=[sum(cliquesize.== i) for i in uc]
#     println("The clique sizes of varibles:\n$uc\n$sizes")
#     return cliques,cql,cliquesize
# end

# function clique_opf_four(n,m,nbus,nb,supp;alg="MF",minimize=true)
#     G=SimpleGraph(n)
#     for i=2:2*nbus+4*nb+1, j = 1:supp[i].n
#         add_clique!(G,supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)])
#     end
#     if alg=="NC"
#         cliques,cql,cliquesize=max_cliques(G)
#     else
#         cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
#     end
#     for i=m-numeq+1:m-numeq+2*nbus
#         cql+=1
#         rind=sort(unique(supp[i+1].rowval))
#         push!(cliques, rind)
#         push!(cliquesize, length(rind))
#     end
#     uc=unique(cliquesize)
#     sizes=[sum(cliquesize.== i) for i in uc]
#     println("The clique sizes of varibles:\n$uc\n$sizes")
#     return cliques,cql,cliquesize
# end

# function clique_opf_four(data,nbus,ng;alg="MF",minimize=true)
#     ref = PowerModels.build_ref(data)[:nw][0]
#     bus=collect(keys(ref[:bus]))
#     sort!(bus)
#     gens=collect(keys(ref[:gen]))
#     sort!(gens)
#     G=SimpleGraph(nbus)
#     for (i, branch) in ref[:branch]
#         vr = bfind(bus,nbus,branch["f_bus"])
#         vt = bfind(bus,nbus,branch["t_bus"])
#         add_edge!(G, vr, vt)
#     end
#     for i in bus
#         temp=UInt16[]
#         for flow in ref[:bus_arcs][i]
#             branch=ref[:branch][flow[1]]
#             vr = bfind(bus,nbus,branch["f_bus"])
#             vt = bfind(bus,nbus,branch["t_bus"])
#             push!(temp, vr, vt)
#         end
#         unique!(temp)
#         add_clique!(G, temp)
#     end
#     # for i=2*nbus+2:4:2*nbus+4*nb-2
#     #     add_clique!(G, [supp[i][1][1];supp[i][1][2]])
#     # end
#     if alg=="NC"
#         cliques,cql,cliquesize=max_cliques(G)
#     else
#         cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
#     end
#     for clique in cliques
#         temp=UInt16[]
#         for r in clique
#             if !isempty(ref[:bus_gens][bus[r]])
#                 for gen_id in ref[:bus_gens][bus[r]]
#                     i=bfind(gens, ng, gen_id)
#                     push!(temp, nbus+i)
#                 end
#             end
#         end
#         if !isempty(temp)
#             sort!(temp)
#             append!(clique, temp)
#         end
#     end
#     # for i in genlabel
#     #     cql+=1
#     #     rind=supp[i][1][1]
#     #     for j=2:length(supp[i])
#     #         append!(rind, supp[i][j][1])
#     #     end
#     #     rind=sort(unique(rind))
#     #     if all(clique -> !(rind ⊆ clique), cliques)
#     #         push!(cliques, rind)
#     #         push!(cliquesize, length(rind))
#     #     end
#     # end
#     cliquesize=length.(cliques)
#     uc=unique(cliquesize)
#     sizes=[sum(cliquesize.== i) for i in uc]
#     println("The clique sizes of varibles:\n$uc\n$sizes")
#     return cliques,cql,cliquesize
# end

# function approx_sol_opf(moment,n,cliques,cql,cliquesize)
#     qsol=Float64[]
#     lcq=sum(cliquesize)
#     A=zeros(lcq,n)
#     q=1
#     for k=1:cql
#         cqs=cliquesize[k]
#         if cqs==1
#             push!(qsol, moment[k])
#         else
#             F=eigen(moment[k], cqs:cqs)
#             temp=sqrt(F.values[1])*F.vectors[:,1]
#             for l=1:k-1
#                 inter=intersect(cliques[l], cliques[k])
#                 if inter!=[]
#                     flag=0
#                     for r in inter
#                         if l==1
#                             ind1=bfind(cliques[l],cliquesize[l],r)
#                         else
#                             ind1=sum(cliquesize[s] for s=1:l-1)+bfind(cliques[l],cliquesize[l],r)
#                         end
#                         ind2=bfind(cliques[k],cqs,r)
#                         if (qsol[ind1]>=1e-3&&temp[ind2]<=-1e-3)||(qsol[ind1]<=-1e-3&&temp[ind2]>=1e-3)
#                             temp=-temp
#                             flag=1
#                             break
#                         elseif (qsol[ind1]>=1e-3&&temp[ind2]>=1e-3)||(qsol[ind1]<=-1e-3&&temp[ind2]<=-1e-3)
#                             flag=1
#                             break
#                         end
#                     end
#                     if flag==1
#                         break
#                     end
#                 end
#             end
#             append!(qsol, temp)
#         end
#         for j=1:cqs
#             A[q,cliques[k][j]]=1
#             q+=1
#         end
#     end
#     return (A'*A)\(A'*qsol)
# end

# Voltage-power real formulization
# function pop_opf(case::String; normal=true)
#     silence()
#     path = joinpath(dirname(dirname(pathof(PolyOPF))), "pglib", "pglib_opf_" * case * ".m")
#     data = parse_file(path)
#     PowerModels.standardize_cost_terms!(data, order=2)
#     ref = PowerModels.build_ref(data)[:nw][0]
#     nbus=length(ref[:bus])
#     nb=length(ref[:branch])
#     ng=length(ref[:gen])
#     n=2*nbus+2*ng
#     m=4*nbus+4*nb+2*ng+length(keys(ref[:ref_buses]))
#     numeq=2*nbus+length(keys(ref[:ref_buses]))
#     dg=2*ones(Int, m)
#     supp=Vector{SparseMatrixCSC{UInt8,UInt32}}(undef, m+1)
#     coe=Vector{Vector{Float64}}(undef, m+1)
#
#     gens=collect(keys(ref[:gen]))
#     sort!(gens)
#     # objective function
#     nc=2*ng+1
#     col=UInt32[i for i=1:nc]
#     col=[1;col]
#     row=UInt32[i+2*nbus for i=1:ng]
#     append!(row,row)
#     nz=ones(UInt8,ng)
#     append!(nz,2*nz)
#     coe[1]=Vector{Float64}(undef, nc)
#     coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
#     for i=1:ng
#         gen=ref[:gen][gens[i]]
#         coe[1][i+1]=gen["cost"][2]
#         coe[1][i+ng+1]=gen["cost"][1]
#     end
#     col,row,nz,coe[1]=move_zero!(col,row,nz,coe[1])
#     supp[1]=SparseMatrixCSC(n,length(coe[1]),col,row,nz)
#     k=2
#
#     bus=collect(keys(ref[:bus]))
#     sort!(bus)
#     # voltage magnitude constraints
#     col=UInt32[1;1;2;3]
#     nz=UInt8[2;2]
#     for i=1:nbus
#         row=UInt32[i;i+nbus]
#         supp[k]=SparseMatrixCSC(n,3,col,row,nz)
#         coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
#         supp[k+1]=SparseMatrixCSC(n,3,col,row,nz)
#         coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
#         k+=2
#     end
#
#     for (i, branch) in ref[:branch]
#         g, b = PowerModels.calc_branch_y(branch)
#         tr, ti = PowerModels.calc_branch_t(branch)
#         g_fr = branch["g_fr"]
#         b_fr = branch["b_fr"]
#         g_to = branch["g_to"]
#         b_to = branch["b_to"]
#         tm = branch["tap"]
#         ab1=-(g+g_fr)^2-(b+b_fr)^2
#         cd1=-(-g*tr+b*ti)^2-(-b*tr-g*ti)^2
#         acbd1=-2*((g+g_fr)*(-g*tr+b*ti)+(b+b_fr)*(-b*tr-g*ti))
#         bcad1=-2*((b+b_fr)*(-g*tr+b*ti)-(g+g_fr)*(-b*tr-g*ti))
#         ab2=-(g+g_to)^2-(b+b_to)^2
#         cd2=-(-g*tr-b*ti)^2/tm^4-(-b*tr+g*ti)^2/tm^4
#         acbd2=-2*((g+g_to)*(-g*tr-b*ti)+(b+b_to)*(-b*tr+g*ti))/tm^2
#         bcad2=-2*((b+b_to)*(-g*tr-b*ti)-(g+g_to)*(-b*tr+g*ti))/tm^2
#         vr_fr = bfind(bus,nbus,branch["f_bus"])
#         vr_to = bfind(bus,nbus,branch["t_bus"])
#         vi_fr = vr_fr+nbus
#         vi_to = vr_to+nbus
#         svr=sort([vr_fr;vr_to])
#         svi=sort([vi_fr;vi_to])
#
#         # angle differences
#         col=UInt32[1;3;5;7;9]
#         nz=ones(UInt8,8)
#         row=UInt32[svr;svi;vr_to;vi_fr;vr_fr;vi_to]
#         coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#         col=UInt32[1;3;5;7;9]
#         nz=ones(UInt8,8)
#         row=UInt32[vr_to;vi_fr;vr_fr;vi_to;svr;svi]
#         coe[k]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#
#         # thermal limits
#         col=UInt32[1;1;2;4;6;8;11;13;16;18;21;24;26;28;29;31;33]
#         row=UInt32[vr_fr;svr;vr_fr;vi_to;svr;svr;vi_fr;vr_fr;vi_fr;vr_fr;svi;vr_fr;vi_to;svr;vi_fr;vr_fr;svi;vr_to;vi_fr;vr_to;vi_fr;vi_fr;svi;svi]
#         if vr_fr<vr_to
#             nz=UInt8[4;3;1;3;1;2;2;2;1;1;2;2;2;1;1;2;2;1;1;2;1;2;1;2;2;1;3;4;3;1;2;2]
#         else
#             nz=UInt8[4;1;3;3;1;2;2;1;2;1;2;2;2;1;1;2;2;1;1;2;1;1;2;2;2;1;3;4;1;3;2;2]
#         end
#         coe[k]=[branch["rate_a"]^2*tm^4;ab1;acbd1;bcad1;cd1;-bcad1;2*ab1;acbd1;cd1;acbd1;bcad1;cd1;-bcad1;ab1;acbd1;cd1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         dg[k-1]=4
#         k+=1
#         col=UInt32[1;1;3;5;7;10;13;15;16;18;20;23;25;28;30;32;33]
#         row=UInt32[svr;vr_fr;vi_to;svr;svr;vi_to;svr;vi_to;vr_fr;vi_to;vr_to;vr_to;vi_fr;vr_to;vi_fr;vr_to;svi;vr_to;vi_to;vr_to;svi;svi;svi;vi_to]
#         if vr_fr<vr_to
#             nz=UInt8[2;2;2;2;1;3;1;2;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;1;2;2;2;1;3;4]
#         else
#             nz=UInt8[2;2;2;2;3;1;2;1;1;1;1;2;1;3;4;3;1;2;2;2;1;1;2;2;1;2;1;2;2;3;1;4]
#         end
#         coe[k]=[branch["rate_a"]^2;cd2;cd2;acbd2;-bcad2;acbd2;-bcad2;ab2;bcad2;cd2;acbd2;2*ab2;bcad2;cd2;acbd2;ab2]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         dg[k-1]=4
#         k+=1
#     end
#
#     # power generation bound
#     zero_pgen=UInt16[]
#     for i=1:ng
#         gen=ref[:gen][gens[i]]
#         if gen["pmax"]>=1e-6
#             col=UInt32[1;1;2;3]
#             nz=UInt8[1;2]
#             row=UInt32[i+2*nbus;i+2*nbus]
#             coe[k]=[-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
#             if normal==true
#                 coe[k]=normalize(coe[k])
#             end
#             col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#             supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#             k+=1
#         else
#             push!(zero_pgen, i)
#             numeq+=1
#         end
#         col=UInt32[1;1;2;3]
#         nz=UInt8[1;2]
#         row=UInt32[i+2*nbus+ng;i+2*nbus+ng]
#         coe[k]=[-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#     end
#
#     # active/reactive power
#     for (r, i) in enumerate(bus)
#         bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
#         bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
#         add_col=UInt32[2;4;6;8]
#         add_nz=UInt8[1;1;1;1;1;1;1;1]
#         col=UInt32[1;1;2;3]
#         row=UInt32[r;r+nbus]
#         nz=UInt8[2;2]
#         coe[k]=zeros(Float64,4*length(ref[:bus_arcs][i])+3)
#         coe[k+1]=zeros(Float64,4*length(ref[:bus_arcs][i])+3)
#         coe[k][1]=fl_sum(load["pd"] for load in bus_loads)
#         coe[k+1][1]=fl_sum(load["qd"] for load in bus_loads)
#         sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
#         sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
#         coe[k][2:3]=[sgs;sgs]
#         coe[k+1][2:3]=[-sbs;-sbs]
#         j=1
#         for flow in ref[:bus_arcs][i]
#             append!(col,col[end].+add_col)
#             append!(nz,add_nz)
#             branch=ref[:branch][flow[1]]
#             vr_fr = bfind(bus,nbus,branch["f_bus"])
#             vr_to = bfind(bus,nbus,branch["t_bus"])
#             vi_fr = vr_fr+nbus
#             vi_to = vr_to+nbus
#             g, b = PowerModels.calc_branch_y(branch)
#             tr, ti = PowerModels.calc_branch_t(branch)
#             g_fr = branch["g_fr"]
#             b_fr = branch["b_fr"]
#             g_to = branch["g_to"]
#             b_to = branch["b_to"]
#             tm = branch["tap"]
#             temp1=(g+g_fr)/tm^2
#             temp2=g+g_to
#             temp3=(-g*tr+b*ti)/tm^2
#             temp4=(-b*tr-g*ti)/tm^2
#             temp5=(-g*tr-b*ti)/tm^2
#             temp6=(-b*tr+g*ti)/tm^2
#             temp7=-(b+b_fr)/tm^2
#             temp8=-(b+b_to)
#             if vr_fr==r
#                 coe[k][2:3].+=temp1
#                 coe[k][4+4*(j-1):7+4*(j-1)]=[temp3;temp3;temp4;-temp4]
#                 coe[k+1][2:3].+=temp7
#                 coe[k+1][4+4*(j-1):7+4*(j-1)]=[-temp4;-temp4;temp3;-temp3]
#             else
#                 coe[k][2:3].+=temp2
#                 coe[k][4+4*(j-1):7+4*(j-1)]=[temp5;temp5;-temp6;temp6]
#                 coe[k+1][2:3].+=temp8
#                 coe[k+1][4+4*(j-1):7+4*(j-1)]=[-temp6;-temp6;-temp5;temp5]
#             end
#             if vr_fr<vr_to
#                 append!(row,[vr_fr;vr_to;vi_fr;vi_to;vr_to;vi_fr;vr_fr;vi_to])
#             else
#                 append!(row,[vr_to;vr_fr;vi_to;vi_fr;vr_to;vi_fr;vr_fr;vi_to])
#             end
#             j+=1
#         end
#         qrow=copy(row)
#         if !isempty(ref[:bus_gens][i])
#             bus_gen=UInt32[]
#             for gen_id in ref[:bus_gens][i]
#                 push!(bus_gen, bfind(gens,ng,gen_id))
#             end
#             lgen=length(bus_gen)
#             append!(coe[k],-ones(Float64,lgen))
#             append!(col,UInt32[l for l=1:lgen].+col[end])
#             append!(row,bus_gen.+2*nbus)
#             append!(nz,ones(UInt8,lgen))
#             append!(coe[k+1],-ones(Float64,lgen))
#             append!(qrow,bus_gen.+(2*nbus+ng))
#         end
#         qcol=copy(col)
#         qnz=copy(nz)
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         if normal==true
#             coe[k+1]=normalize(coe[k+1])
#         end
#         qcol,qrow,qnz,coe[k+1]=move_zero!(qcol,qrow,qnz,coe[k+1])
#         supp[k+1]=SparseMatrixCSC(n,length(coe[k+1]),qcol,qrow,qnz)
#         k+=2
#     end
#
#     # reference voltage
#     for key in keys(ref[:ref_buses])
#         i=bfind(bus,nbus,key)
#         supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+nbus],UInt8[2])
#         coe[k]=[1]
#         k+=1
#     end
#
#     # zero power generation
#     for i in zero_pgen
#         supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+2*nbus],UInt8[1])
#         coe[k]=[1]
#         dg[k-1]=1
#         k+=1
#     end
#     return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg)
# end

# Voltage-power real QCQP formulization
# function pop_opf_two(case::String; normal=true)
#     silence()
#     path = joinpath(dirname(dirname(pathof(PolyOPF))), "pglib", "pglib_opf_" * case * ".m")
#     data = parse_file(path)
#     PowerModels.standardize_cost_terms!(data, order=2)
#     ref = PowerModels.build_ref(data)[:nw][0]
#     nbus=length(ref[:bus])
#     ng=length(ref[:gen])
#     nb=length(ref[:branch])
#     n=2*nbus+2*ng+4*nb
#     m=4*nbus+8*nb+2*ng+length(ref[:ref_buses])
#     numeq=2*nbus+4*nb+length(ref[:ref_buses])
#     dg=2*ones(Int, m)
#     supp=Vector{SparseMatrixCSC{UInt8,UInt32}}(undef, m+1)
#     coe=Vector{Vector{Float64}}(undef, m+1)
#
#     gens=collect(keys(ref[:gen]))
#     sort!(gens)
#     # objective function
#     nc=2*ng+1
#     col=UInt32[i for i=1:nc]
#     col=[1;col]
#     row=UInt32[i+2*nbus for i=1:ng]
#     append!(row,row)
#     nz=ones(UInt8,ng)
#     append!(nz,2*nz)
#     coe[1]=Vector{Float64}(undef, nc)
#     coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
#     for i=1:ng
#         gen=ref[:gen][gens[i]]
#         coe[1][i+1]=gen["cost"][2]
#         coe[1][i+ng+1]=gen["cost"][1]
#     end
#     col,row,nz,coe[1]=move_zero!(col,row,nz,coe[1])
#     supp[1]=SparseMatrixCSC(n,length(coe[1]),col,row,nz)
#
#     bus=collect(keys(ref[:bus]))
#     sort!(bus)
#     # voltage magnitude constraints
#     k=2
#     col=UInt32[1;1;2;3]
#     nz=UInt8[2;2]
#     for i=1:nbus
#         row=UInt32[i;i+nbus]
#         supp[k]=SparseMatrixCSC(n,3,col,row,nz)
#         coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
#         supp[k+1]=SparseMatrixCSC(n,3,col,row,nz)
#         coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
#         k+=2
#     end
#
#     branchs=collect(keys(ref[:branch]))
#     sort!(branchs)
#     for i=1:nb
#         branch=ref[:branch][branchs[i]]
#         vr_fr = bfind(bus,nbus,branch["f_bus"])
#         vr_to = bfind(bus,nbus,branch["t_bus"])
#         vi_fr = vr_fr+nbus
#         vi_to = vr_to+nbus
#         svr=sort([vr_fr;vr_to])
#         svi=sort([vi_fr;vi_to])
#
#         # angle differences
#         col=UInt32[1;3;5;7;9]
#         nz=ones(UInt8,8)
#         row=UInt32[svr;svi;vr_to;vi_fr;vr_fr;vi_to]
#         coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#         col=UInt32[1;3;5;7;9]
#         nz=ones(UInt8,8)
#         row=UInt32[vr_to;vi_fr;vr_fr;vi_to;svr;svi]
#         coe[k]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#
#         # thermal limits
#         col=UInt32[1;1;2;3]
#         row=UInt32[2*nbus+2*ng+i;2*nbus+2*ng+nb+i]
#         nz=[2;2]
#         coe[k]=[branch["rate_a"]^2;-1;-1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         supp[k]=SparseMatrixCSC(n,3,col,row,nz)
#         coe[k+1]=coe[k]
#         row=UInt32[2*nbus+2*ng+2*nb+i;2*nbus+2*ng+3*nb+i]
#         supp[k+1]=SparseMatrixCSC(n,3,col,row,nz)
#         k+=2
#    end
#
#     # power generation bound
#     zero_pgen=UInt16[]
#     for i=1:ng
#         gen=ref[:gen][gens[i]]
#         if gen["pmax"]>=1e-6
#             col=UInt32[1;1;2;3]
#             nz=UInt8[1;2]
#             row=UInt32[i+2*nbus;i+2*nbus]
#             coe[k]=[-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
#             if normal==true
#                 coe[k]=normalize(coe[k])
#             end
#             col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#             supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#             k+=1
#         else
#             push!(zero_pgen, i)
#             numeq+=1
#         end
#         col=UInt32[1;1;2;3]
#         nz=UInt8[1;2]
#         row=UInt32[i+2*nbus+ng;i+2*nbus+ng]
#         coe[k]=[-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#     end
#
#     for i=1:nb
#         branch=ref[:branch][branchs[i]]
#         g, b = PowerModels.calc_branch_y(branch)
#         tr, ti = PowerModels.calc_branch_t(branch)
#         g_fr = branch["g_fr"]
#         b_fr = branch["b_fr"]
#         g_to = branch["g_to"]
#         b_to = branch["b_to"]
#         tm = branch["tap"]
#         vr_fr = bfind(bus,nbus,branch["f_bus"])
#         vr_to = bfind(bus,nbus,branch["t_bus"])
#         vi_fr = vr_fr+nbus
#         vi_to = vr_to+nbus
#         svr=sort([vr_fr;vr_to])
#         svi=sort([vi_fr;vi_to])
#
#         # Line Flow
#         col=UInt32[1;2;3;5;7;9;11;12]
#         nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
#         row=UInt32[vr_fr;vi_fr;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+i]
#         coe[k]=[g+g_fr;g+g_fr;-g*tr+b*ti;-g*tr+b*ti;-b*tr-g*ti;b*tr+g*ti;-tm^2]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#         col=UInt32[1;2;3;5;7;9;11;12]
#         nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
#         row=UInt32[vr_fr;vi_fr;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+nb+i]
#         coe[k]=[-(b+b_fr);-(b+b_fr);b*tr+g*ti;b*tr+g*ti;-g*tr+b*ti;g*tr-b*ti;-tm^2]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#         col=UInt32[1;2;3;5;7;9;11;12]
#         nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
#         row=UInt32[vr_to;vi_to;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+2*nb+i]
#         coe[k]=[(g+g_to)*tm^2;(g+g_to)*tm^2;-g*tr-b*ti;-g*tr-b*ti;b*tr-g*ti;-b*tr+g*ti;-tm^2]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#         col=UInt32[1;2;3;5;7;9;11;12]
#         nz=UInt8[2;2;1;1;1;1;1;1;1;1;1]
#         row=UInt32[vr_to;vi_to;svr;svi;vr_to;vi_fr;vr_fr;vi_to;2*nbus+2*ng+3*nb+i]
#         coe[k]=[-(b+b_to)*tm^2;-(b+b_to)*tm^2;b*tr-g*ti;b*tr-g*ti;g*tr+b*ti;-g*tr-b*ti;-tm^2]
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         k+=1
#    end
#
#     # active/reactive power
#     for r=1:nbus
#         i=bus[r]
#         bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
#         bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
#         lg=length(ref[:bus_arcs][i])+length(ref[:bus_gens][i])+2
#         col=UInt32[i for i=1:lg+1]
#         col=UInt32[1;col]
#         row=zeros(UInt32,lg)
#         row[1]=r
#         row[2]=r+nbus
#         qrow=copy(row)
#         nz=ones(UInt8,lg)
#         nz[1]=2
#         nz[2]=2
#         coe[k]=ones(Float64,lg+1)
#         coe[k+1]=ones(Float64,lg+1)
#         coe[k][1]=fl_sum(load["pd"] for load in bus_loads)
#         coe[k+1][1]=fl_sum(load["qd"] for load in bus_loads)
#         sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
#         sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)
#         coe[k][2]=sgs
#         coe[k][3]=sgs
#         coe[k+1][2]=-sbs
#         coe[k+1][3]=-sbs
#         j=3
#         for flow in ref[:bus_arcs][i]
#             branch=ref[:branch][flow[1]]
#             s=bfind(branchs,nb,flow[1])
#             vr_fr = bfind(bus,nbus,branch["f_bus"])
#             if vr_fr==r
#                 row[j]=2*nbus+2*ng+s
#                 qrow[j]=2*nbus+2*ng+nb+s
#             else
#                 row[j]=2*nbus+2*ng+2*nb+s
#                 qrow[j]=2*nbus+2*ng+3*nb+s
#             end
#             j+=1
#         end
#         for gen_id in ref[:bus_gens][i]
#             gen=bfind(gens,ng,gen_id)
#             coe[k][j+1]=-1
#             coe[k+1][j+1]=-1
#             row[j]=gen+2*nbus
#             qrow[j]=gen+2*nbus+ng
#             j+=1
#         end
#         qcol=copy(col)
#         qnz=copy(nz)
#         if normal==true
#             coe[k]=normalize(coe[k])
#         end
#         col,row,nz,coe[k]=move_zero!(col,row,nz,coe[k])
#         supp[k]=SparseMatrixCSC(n,length(coe[k]),col,row,nz)
#         if normal==true
#             coe[k+1]=normalize(coe[k+1])
#         end
#         qcol,qrow,qnz,coe[k+1]=move_zero!(qcol,qrow,qnz,coe[k+1])
#         supp[k+1]=SparseMatrixCSC(n,length(coe[k+1]),qcol,qrow,qnz)
#         k+=2
#     end
#
#     # reference voltage
#     for key in keys(ref[:ref_buses])
#         i=bfind(bus,nbus,key)
#         supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+nbus],UInt8[2])
#         coe[k]=[1]
#         k+=1
#     end
#
#     # zero power generation
#     for i in zero_pgen
#         supp[k]=SparseMatrixCSC(n,1,UInt32[1;2],UInt32[i+2*nbus],UInt8[1])
#         coe[k]=[1]
#         dg[k-1]=1
#         k+=1
#     end
#     return SparsePolyModel(n,m,numeq,nbus,ng,nb,supp,coe,dg)
# end