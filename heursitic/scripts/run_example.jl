# Chạy minh họa: Ví dụ 1 của bài báo + một instance Li & Lim.
#   julia scripts/run_example.jl
using Random
include(joinpath(@__DIR__, "..", "src", "PDPCR.jl"))
using .PDPCR

println("="^70, "\nVí dụ 1 (n = 3, m = 2)\n", "="^70)
inst = example1()
g = greedy_solution(inst)
println("Lời giải tham lam:")
print_schedule(inst, g)
res = tabu_search(inst, TabuParams(:TSR, inst.n; max_iter = 200))
println("\nSau tabu search (TSR):")
print_schedule(inst, res.best)

println("\n", "="^70, "\nInstance Li & Lim lc101 (53 task, m = 10, r_i ~ U{1..5})\n", "="^70)
inst = read_lilim(joinpath(@__DIR__, "..", "data", "lilim100", "lc101.txt"); m = 10, seed = 1)
init = greedy_solution(inst)
for v in (:TSO, :TSR)
    r = tabu_search(inst, TabuParams(v, inst.n; max_iter = 2000); init = init,
                    rng = Random.Xoshiro(42), verbose = true)
    println(rpad(string(v), 4), ": ban đầu = ", round(r.initial_value; digits = 2),
            "  best = ", round(r.best_value; digits = 2),
            "  cải thiện = ", round(100 * (1 - r.best_value / r.initial_value); digits = 2), "%",
            "  thời gian = ", round(r.time; digits = 1), "s",
            "  khả thi = ", check_feasible(inst, r.best))
    println("      nước được chọn: ", r.accepted, "   cải thiện best nhờ: ", r.improved_by)
end
