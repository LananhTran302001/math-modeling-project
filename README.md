# PDP-CR — Thuật toán heuristic (Julia)

Phần việc của **Người 2**: tabu search, các biến thể TSO/TSR, và lân cận *redistribute* cho bài toán Pickup and Delivery with Cooperative Robots (`Thao_ICAMCS_special_issue.pdf`).

Code chỉ dùng thư viện chuẩn của Julia (Random, Printf, Statistics, Test). Đã chạy thử với Julia 1.13.

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
