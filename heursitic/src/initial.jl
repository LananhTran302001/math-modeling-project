# ============================================================================
# Lời giải ban đầu: tham lam theo thời điểm bắt đầu sớm nhất (list scheduling)
#
# Mỗi robot k có (vị trí cuối last_k, thời điểm rảnh avail_k). Lặp n lần:
#   - với mỗi task j chưa làm: chọn r_j robot có thời điểm đến p_j sớm nhất
#     (avail_k + t[last_k, j]); thời điểm bắt đầu khả dĩ = max các thời điểm đến đó;
#   - chọn task có thời điểm bắt đầu nhỏ nhất (hoặc ngẫu nhiên trong top-`rcl`
#     ứng viên tốt nhất khi cần đa dạng hóa / restart);
#   - gán robot, thêm luồng f[last_k, j] += 1.
# Cuối cùng mọi robot về T. Đồ thị thu được không có chu trình vì luồng luôn
# đi theo chiều thời gian tăng.
# ============================================================================

function greedy_solution(inst::Instance; rng = Random.default_rng(), rcl::Int = 1)
    n, m = inst.n, inst.m
    S, T = source(inst), sink(inst)
    last = fill(S, m)
    avail = zeros(m)
    done = falses(n)
    sol = Solution(nnodes(inst))
    arrive = zeros(m)
    perm = collect(1:m)
    cand = Tuple{Float64,Int}[]
    for _ in 1:n
        empty!(cand)
        for j in 1:n
            done[j] && continue
            for k in 1:m
                arrive[k] = avail[k] + inst.t[last[k], j]
            end
            partialsortperm!(perm, arrive, inst.r[j])
            push!(cand, (arrive[perm[inst.r[j]]], j))
        end
        sort!(cand)
        _, j = cand[rand(rng, 1:min(rcl, length(cand)))]
        for k in 1:m
            arrive[k] = avail[k] + inst.t[last[k], j]
        end
        partialsortperm!(perm, arrive, inst.r[j])
        start = arrive[perm[inst.r[j]]]
        for q in 1:inst.r[j]
            k = perm[q]
            add_flow!(sol, last[k], j, 1)
            last[k] = j
            avail[k] = start + inst.c[j]
        end
        perm .= 1:m
        done[j] = true
    end
    for k in 1:m
        add_flow!(sol, last[k], T, 1)
    end
    return sol
end
