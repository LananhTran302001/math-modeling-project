# Pickup and Delivery Problem with Cooperative Robots (PDP-CR)

This repository contains implementations of exact and heuristic approaches for solving the **Pickup and Delivery Problem with Cooperative Robots (PDP-CR)**.

The project investigates two main approaches:

- **Mixed-Integer Linear Programming (MILP):** An exact optimization formulation based on robot flows and task scheduling.
- **Tabu Search:** A heuristic approach with three neighborhood operators:
  - **Reroute:** Reorganizes robot routes by exchanging successors.
  - **Reverse:** Reverses selected arcs to explore alternative task orders.
  - **Redistribute:** Reallocates robot flows while preserving task order.

Two Tabu Search variants are implemented:

- **TSO:** Baseline Tabu Search using Reroute and Reverse.
- **TSR:** Extended Tabu Search incorporating the proposed Redistribute neighborhood.

## Experiments

The algorithms are evaluated on benchmark instances adapted from Li & Lim and randomly generated instances. Experimental results compare solution quality and computational runtime.

## Technologies

- Julia
- JuMP
- HiGHS

## Cấu trúc

```
heuristic/
├── src/
│   ├── PDPCR.jl          # module chính
│   ├── instance.jl       # dữ liệu, Ví dụ 1, đọc Li & Lim, sinh ngẫu nhiên, ghi/đọc file
│   ├── solution.jl       # đồ thị luồng, Thuật toán 1-2 (L, L*), cung găng, kiểm tra khả thi, Gantt
│   ├── initial.jl        # lời giải ban đầu tham lam
│   ├── neighborhoods.jl  # reroute, reverse, redistribute (Thuật toán 3)
│   └── tabu.jl           # tabu search, tham số TSO/TSR (Table 4)
├── test/runtests.jl      # kiểm tra với Ví dụ 1, Ví dụ 2-3 và tính khả thi của mọi nước đi
├── scripts/
│   ├── run_example.jl    # chạy minh họa
│   ├── run_experiments.jl# so sánh TSO và TSR, xuất bảng kiểu Table 5 và dữ liệu Hình 12
│   └── plot_history.py   # vẽ đồ thị kiểu Hình 12 (matplotlib)
├── results/              # kết quả thực nghiệm
```

## Chạy

```bash
cd heuristic
julia test/runtests.jl                                   # kiểm thử
julia scripts/run_example.jl                             # minh họa
julia scripts/run_experiments.jl lilim 5000 600 10       # 10 instance Li & Lim, 5000 vòng lặp
julia scripts/run_experiments.jl random 5000 1800 10 100,200,300,400,500
```

## Dùng trong code

```julia
include("src/PDPCR.jl"); using .PDPCR, Random
inst = read_lilim("data/lilim100/lc101.txt"; m = 10, seed = 1)
res  = tabu_search(inst, TabuParams(:TSR, inst.n; max_iter = 5000, time_limit = 600);
                   rng = Random.Xoshiro(1), verbose = true)
res.best_value, check_feasible(inst, res.best)
print_schedule(inst, res.best)          # lộ trình và lịch của từng robot
```

## Team Contributions

| Contributor | Main Responsibilities | Report Sections |
|---|---|---|
| Tran Thi Lan Anh | Studied the problem and established its assumptions; defined the notation, data, decision variables, objective function, and constraints; implemented the MILP model; validated feasibility on small instances; conducted exact optimization experiments and recorded solutions, runtime, MIP gap, and node counts. | Introduction; Mathematical Model; Exact MILP Approach; Exact Method Results |
| Do Hoang Viet | Designed the initial solution and solution representation; implemented Tabu Search, TSO/TSR variants, and the Redistribute neighborhood based on the team's agreed definitions; validated solution feasibility and stopping criteria; conducted heuristic experiments and recorded solution quality and runtime. | Tabu Search; TSO and TSR; Redistribute Neighborhood; Heuristic Results |
| Joint Contributions | Agreed on benchmark datasets; established and coordinated the experimental setup; ran and compared both approaches on the same instances under consistent time limits; analyzed results and limitations; reviewed the report and prepared the presentation slides. | Numerical Experiments; Discussion; Conclusion |
