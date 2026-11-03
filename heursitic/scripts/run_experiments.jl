# So sánh TSO và TSR theo giao thức của bài báo (Mục 5.2, Table 5).
#
#   julia scripts/run_experiments.jl lilim  [iters] [time_limit_s] [số_instance]
#   julia scripts/run_experiments.jl random [iters] [time_limit_s] [số_instance] [n1,n2,...]
#
# Ví dụ:
#   julia scripts/run_experiments.jl lilim 5000 600 10
#   julia scripts/run_experiments.jl random 5000 3600 10 100,200,300
#
# Kết quả ghi vào results/: results_<set>.csv (từng instance),
# summary_<set>.md (bảng kiểu Table 5), history_<set>.csv (đường hội tụ, kiểu Hình 12).
using Random, Printf, Statistics
include(joinpath(@__DIR__, "..", "src", "PDPCR.jl"))
using .PDPCR

const ROOT = joinpath(@__DIR__, "..")
set     = length(ARGS) >= 1 ? ARGS[1] : "lilim"
iters   = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5000
tlimit  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : Inf
nins    = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 10
sizes   = length(ARGS) >= 5 ? parse.(Int, split(ARGS[5], ",")) : [100, 200]
const M = 10

# ---- danh sách instance (cùng seed để cả nhóm dùng chung) -------------------
instances = Instance[]
if set == "lilim"
    files = sort(readdir(joinpath(ROOT, "data", "lilim100"); join = true))
    for f in files[1:min(nins, length(files))]
        push!(instances, read_lilim(f; m = M, seed = 1))
    end
else
    for n in sizes, k in 1:nins
        push!(instances, random_instance(n, M; seed = k))
    end
end
mkpath(joinpath(ROOT, "results"))
mkpath(joinpath(ROOT, "data", "generated"))
for inst in instances     # lưu lại để người làm MILP dùng đúng cùng dữ liệu
    write_instance(joinpath(ROOT, "data", "generated", inst.name * ".txt"), inst)
end

# ---- chạy -------------------------------------------------------------------
rows = []
hist = Dict(:TSO => Vector{Vector{Float64}}(), :TSR => Vector{Vector{Float64}}())
for inst in instances
    init = greedy_solution(inst)               # cùng lời giải ban đầu cho cả hai
    out = Dict{Symbol,Any}()
    for v in (:TSO, :TSR)
        P = TabuParams(v, inst.n; max_iter = iters, time_limit = tlimit)
        r = tabu_search(inst, P; init = init, rng = Random.Xoshiro(2024))
        @assert check_feasible(inst, r.best)
        out[v] = r
        push!(hist[v], r.history ./ r.initial_value .* 100)
    end
    a, b = out[:TSO], out[:TSR]
    diff = 100 * (a.best_value - b.best_value) / a.best_value    # % TSR tốt hơn TSO
    push!(rows, (name = inst.name, n = inst.n, cls = 50 * round(Int, inst.n / 50), m = inst.m, init = a.initial_value,
                 tso = a.best_value, tso_t = a.time, tso_it = a.iterations,
                 tsr = b.best_value, tsr_t = b.time, tsr_it = b.iterations, diff = diff))
    @printf("%-22s n=%3d  init=%9.2f  TSO=%9.2f (%6.1fs)  TSR=%9.2f (%6.1fs)  diff=%+6.2f%%\n",
            inst.name, inst.n, a.initial_value, a.best_value, a.time, b.best_value, b.time, diff)
end

# ---- ghi kết quả --------------------------------------------------------------
open(joinpath(ROOT, "results", "results_$(set).csv"), "w") do io
    println(io, "instance,n,m,init,tso,tso_time,tso_iters,tsr,tsr_time,tsr_iters,diff_pct")
    for r in rows
        println(io, join((r.name, r.n, r.m, r.init, r.tso, r.tso_t, r.tso_it, r.tsr, r.tsr_t, r.tsr_it, r.diff), ","))
    end
end

open(joinpath(ROOT, "results", "summary_$(set).md"), "w") do io
    println(io, "| #ins | n | m | TSO time(s) | TSO impr. vs init (%) | TSR impr. vs init (%) | diff avg | diff min | diff max | TSR time(s) | #better |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
    for n in sort(unique(r.cls for r in rows))          # nhóm theo cỡ danh nghĩa (Li & Lim: n ≈ 50)
        g = [r for r in rows if r.cls == n]
        d = [r.diff for r in g]
        @printf(io, "| %d | %d | %d | %.1f | %.2f | %.2f | %.2f | %.2f | %.2f | %.1f | %d |\n",
                length(g), n, M, mean(r.tso_t for r in g),
                mean(100 * (1 - r.tso / r.init) for r in g), mean(100 * (1 - r.tsr / r.init) for r in g),
                mean(d), minimum(d), maximum(d), mean(r.tsr_t for r in g), count(>(1e-9), d))
    end
end

open(joinpath(ROOT, "results", "history_$(set).csv"), "w") do io
    L = minimum(length, vcat(hist[:TSO], hist[:TSR]))
    println(io, "iteration,TSO,TSR")
    for k in 1:L
        @printf(io, "%d,%.4f,%.4f\n", k, mean(h[k] for h in hist[:TSO]), mean(h[k] for h in hist[:TSR]))
    end
end
println("\nĐã ghi kết quả vào ", joinpath(ROOT, "results"))
print(read(joinpath(ROOT, "results", "summary_$(set).md"), String))
