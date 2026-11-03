using Test, Random
include(joinpath(@__DIR__, "..", "src", "PDPCR.jl"))
using .PDPCR

const DATA = joinpath(@__DIR__, "..", "data")

@testset "PDP-CR heuristic" begin

    @testset "Ví dụ 1: nghiệm Hình 3 có makespan 21.5" begin
        inst = example1()
        S, T = source(inst), sink(inst)
        s = Solution(nnodes(inst))
        # robot 1: S -> 2 -> 1 -> T ; robot 2: S -> 3 -> 1 -> T
        for (a, b) in [(S, 2), (2, 1), (1, T), (S, 3), (3, 1), (1, T)]
            add_flow!(s, a, b, 1)
        end
        @test check_feasible(inst, s)
        @test makespan!(Workspace(nnodes(inst)), inst, s) ≈ 21.5
        g = greedy_solution(inst)
        @test check_feasible(inst, g)
        res = tabu_search(inst, TabuParams(:TSR, inst.n; max_iter = 200))
        @test check_feasible(inst, res.best)
        @test res.best_value <= 21.5 + 1e-9
    end

    @testset "Ví dụ 2-3: L̄, L* (Hình 6) và redistribute (Hình 9)" begin
        # node: A=1 B=2 C=3 D=4 E=5 S=6 S'=7, m = 7
        r = [5, 2, 2, 4, 3]; c = [1.0, 2, 4, 3, 1]
        t = fill(100.0, 7, 7)
        for ((a, b), v) in Dict((6, 1) => 4, (6, 2) => 7, (1, 4) => 10, (1, 3) => 5, (2, 3) => 3,
                                (2, 5) => 3, (3, 5) => 4, (4, 7) => 2, (5, 7) => 3,
                                (1, 5) => 4, (2, 4) => 2, (3, 4) => 1)   # Table 3
            t[a, b] = v
        end
        inst = Instance("ex2", 7, r, c, t)
        s = Solution(7)
        for (a, b, f) in [(6, 1, 5), (6, 2, 2), (1, 4, 4), (1, 3, 1), (2, 3, 1), (2, 5, 1),
                          (3, 5, 2), (4, 7, 4), (5, 7, 3)]
            add_flow!(s, a, b, f)
        end
        @test check_feasible(inst, s)
        ws = Workspace(7)
        cmax = makespan!(ws, inst, s); tails!(ws, inst, s)
        # Lưu ý: Hình 5 ghi t_AD = 6 nhưng Table 3 ghi t_AD = 10; ta dùng Table 3
        # nên L̄_D, L*_A khác hình, còn các giá trị dưới đây trùng Hình 6/8.
        @test [ws.L[i] + c[i] for i in (1, 2, 3)] ≈ [5, 9, 16]           # L̄_A, L̄_B, L̄_C
        @test [ws.Lstar[j] for j in (4, 5, 7)] ≈ [5, 4, 0]               # L*_D, L*_E, L*_S'
        @test cmax ≈ 24
        order = [6, 2, 1, 3, 4, 5, 7]                                   # S B A C | D E S'
        mv, newc = redistribute(inst, ws, s, order, 4, cmax, Random.Xoshiro(1))
        @test newc ≈ 22
        apply!(s, mv.changes)
        @test check_feasible(inst, s)
        @test makespan!(ws, inst, s) ≈ 22
        @test (flow(s, 1, 4), flow(s, 1, 5), flow(s, 3, 4), flow(s, 3, 5), flow(s, 2, 5)) == (2, 2, 2, 0, 1)
    end

    @testset "Các nước đi giữ tính khả thi + redis đánh giá đúng" begin
        inst = random_instance(30, 6; seed = 7)
        rng = Random.Xoshiro(3)
        s = greedy_solution(inst)
        @test check_feasible(inst, s)
        ws = Workspace(nnodes(inst)); ws2 = Workspace(nnodes(inst))
        for _ in 1:30
            cmax = makespan!(ws, inst, s); tails!(ws, inst, s)
            crit = critical_arcs(ws, inst, s, cmax)
            @test !isempty(crit)
            for mv in vcat(reroute_moves(inst, s, crit, 20, rng), reverse_moves(inst, ws, s, crit, 20, rng))
                apply!(s, mv.changes)
                v = makespan!(ws2, inst, s)
                isfinite(v) && @test check_feasible(inst, s)
                undo!(s, mv.changes)
            end
            order = topo_by_start(ws, inst)
            for ζ in 1:nnodes(inst)-1
                res = redistribute(inst, ws, s, order, ζ, cmax, rng)
                res === nothing && continue
                mv, v = res
                apply!(s, mv.changes)
                @test check_feasible(inst, s)
                @test makespan!(ws2, inst, s) ≈ v
                @test v <= cmax + 1e-9
                undo!(s, mv.changes)
            end
            # đi một bước ngẫu nhiên để thăm nghiệm khác
            mvs = reroute_moves(inst, s, crit, 5, rng)
            for mv in mvs
                apply!(s, mv.changes)
                isfinite(makespan!(ws2, inst, s)) && break
                undo!(s, mv.changes)
            end
        end
    end

    @testset "TSO / TSR trên instance Li & Lim" begin
        inst = read_lilim(joinpath(DATA, "lilim100", "lc101.txt"); m = 10, seed = 1)
        @test inst.n == 53
        for v in (:TSO, :TSR)
            res = tabu_search(inst, TabuParams(v, inst.n; max_iter = 300); rng = Random.Xoshiro(1))
            @test check_feasible(inst, res.best)
            @test res.best_value <= res.initial_value + 1e-9
            @test makespan!(Workspace(nnodes(inst)), inst, res.best) ≈ res.best_value
        end
    end

    @testset "Ghi / đọc instance" begin
        inst = random_instance(8, 3; seed = 2)
        p = tempname()
        write_instance(p, inst)
        inst2 = read_instance(p)
        @test inst2.r == inst.r && inst2.c ≈ inst.c && inst2.t ≈ inst.t
    end
end
