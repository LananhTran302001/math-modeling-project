# ============================================================================
# Các lân cận trên đồ thị luồng
#   - N_reroute^{max,ca,δ} (Poppenborg & Knust 2015, Mục 3.2 bài báo)
#   - N_reverse^{ca}       (Mục 3.2)
#   - N_redis              (Mục 3.1, Thuật toán 3) — đóng góp chính của bài báo
#
# Mỗi lân cận sinh ra các `Move`: danh sách thay đổi luồng + thông tin cho tabu.
# ============================================================================

struct Move
    kind::Symbol                     # :reroute, :reverse, :redis
    changes::Changes                 # (i, j, Δf)
    removed::Vector{Tuple{Int,Int}}  # cung bị giảm luồng (đưa vào TL_add)
    added::Vector{Tuple{Int,Int}}    # cung được tăng luồng (dùng để kiểm tra tabu)
    cut::Int                         # ζ (chỉ dùng cho :redis)
end

# ----------------------------------------------------------------------------
# Reroute: chọn cung găng (i,j) có t_ij lớn nhất và một cung (u,v) bất kỳ,
# chuyển q = min(f_ij, f_uv) đơn vị:  f_ij -= q, f_uv -= q, f_iv += q, f_uj += q.
# ----------------------------------------------------------------------------
function reroute_move(s::Solution, i, j, u, v)
    (u == i || v == j || i == v || u == j) && return nothing
    q = min(flow(s, i, j), flow(s, u, v))
    q > 0 || return nothing
    ch = Changes([(i, j, -q), (u, v, -q), (i, v, q), (u, j, q)])
    return Move(:reroute, ch, [(i, j), (u, v)], [(i, v), (u, j)], 0)
end

"""
    reroute_moves(inst, s, crit, k, rng)

N_reroute^{max,ca,δ}: cung găng (i,j) có t_ij lớn nhất; tối đa `k` cung (u,v) ngẫu nhiên.

Nếu `L !== nothing` (mặc định trong tabu search) chỉ lấy mẫu các (u,v) thỏa điều kiện
cần để không tạo chu trình: L_v ≥ L_i và L_j ≥ L_u (nếu đã có đường v ⇝ i thì
L_i ≥ L_v + c_v + t > L_v). Không có bộ lọc này, thực nghiệm cho thấy ~97% nước reroute
ngẫu nhiên tạo chu trình và bị loại.
"""
function reroute_moves(inst::Instance, s::Solution, crit, k::Int, rng; L = nothing,
                       all_crit::Bool = false)
    moves = Move[]
    (k <= 0 || isempty(crit)) && return moves
    # chỉ xét cung găng thật sự đổi được (không phải S->T); t lớn nhất trước, hòa thì ngẫu nhiên
    cands = shuffle(rng, filter(a -> !(a[1] == source(inst) && a[2] == sink(inst)), crit))
    isempty(cands) && return moves
    sort!(cands; by = a -> -inst.t[a[1], a[2]])
    # all_crit = false: đúng N^{max,ca,δ} (chỉ 1 cung găng);
    # all_crit = true : nếu chưa đủ k ứng viên thì xét tiếp các cung găng khác (N^{max,ca})
    for (i, j) in (all_crit ? cands : cands[1:1])
        all = [(u, v) for u in eachindex(s.succ) for v in keys(s.succ[u])
               if L === nothing || (L[v] >= L[i] && L[j] >= L[u])]
        shuffle!(rng, all)
        for (u, v) in all
            mv = reroute_move(s, i, j, u, v)
            mv === nothing || push!(moves, mv)
            length(moves) >= k && return moves
        end
    end
    return moves
end

# ----------------------------------------------------------------------------
# Reverse: đảo chiều cung găng (i,j) giữa hai task.
#   U = các tiền bối u của i (f_ui > 0), ưu tiên thời điểm robot đến i sớm,
#       lấy tối thiểu sao cho Σ f_ui ≥ f_ij; chuyển luồng u -> i thành u -> j.
#   V = các hậu bối v của j (f_jv > 0), ưu tiên t_jv nhỏ,
#       lấy tối thiểu sao cho Σ f_jv ≥ f_ij; chuyển luồng j -> v thành i -> v.
#   Sau đó f_ji = f_ij, f_ij = 0.
# ----------------------------------------------------------------------------
function reverse_move(inst::Instance, ws::Workspace, s::Solution, i::Int, j::Int)
    (i > inst.n || j > inst.n) && return nothing   # không đảo cung chạm depot
    fij = flow(s, i, j)
    fij > 0 || return nothing
    ch = Changes()
    push!(ch, (i, j, -fij)); push!(ch, (j, i, fij))
    # tập U
    us = sort!(collect(keys(s.pred[i])); by = u -> ws.L[u] + inst.c[u] + inst.t[u, i])
    rest = fij
    for u in us
        rest == 0 && break
        q = min(rest, s.pred[i][u])
        push!(ch, (u, i, -q)); push!(ch, (u, j, q))
        rest -= q
    end
    # tập V
    vs = sort!(collect(keys(s.succ[j])); by = v -> inst.t[j, v])
    rest = fij
    for v in vs
        rest == 0 && break
        q = min(rest, s.succ[j][v])
        push!(ch, (j, v, -q)); push!(ch, (i, v, q))
        rest -= q
    end
    # gộp các thay đổi trùng cung (vd. u == ... ) để apply/undo an toàn
    agg = Dict{Tuple{Int,Int},Int}()
    for (a, b, d) in ch
        agg[(a, b)] = get(agg, (a, b), 0) + d
    end
    any(((a, b),) -> a == b, keys(agg)) && return nothing
    chg = Changes([(a, b, d) for ((a, b), d) in agg if d != 0])
    sort!(chg; by = x -> x[3])        # giảm trước, tăng sau
    added = [(a, b) for (a, b, d) in chg if d > 0]
    return Move(:reverse, chg, [(i, j)], added, 0)
end

"N_reverse^{ca}: tối đa `k` cung găng (giữa hai task) được chọn ngẫu nhiên."
function reverse_moves(inst::Instance, ws::Workspace, s::Solution, crit, k::Int, rng)
    moves = Move[]
    k <= 0 && return moves
    for (i, j) in shuffle(rng, crit)
        mv = reverse_move(inst, ws, s, i, j)
        mv === nothing || push!(moves, mv)
        length(moves) >= k && break
    end
    return moves
end

# ----------------------------------------------------------------------------
# Redistribute (Thuật toán 3)
#
# Cắt thứ tự topo tại ζ: N_left = ζ node đầu (chứa S), N_right = phần còn lại (chứa T).
# Mọi đường S -> T cắt lát cắt đúng một lần, nên sau khi phân phối lại luồng
# qua lát cắt (giữ nguyên π^left_i, π^right_j) thì
#       makespan mới = max_{(i,j) cắt, f_ij > 0} Δ_ij,  Δ_ij = L̄_i + t_ij + L*_j,
# và đồ thị vẫn không có chu trình (thứ tự topo cũ vẫn hợp lệ).
# Tổng luồng qua lát cắt đúng bằng m, nên bài toán (RP) chỉ có ≤ m node mỗi bên.
# ----------------------------------------------------------------------------
"""
    redistribute(inst, ws, s, order, ζ, cmax, rng) -> (Move, Cmax_mới) hoặc nothing

Chạy Thuật toán 3 trên lát cắt ζ. Cần `ws.L`, `ws.Lstar` của nghiệm hiện tại.
"""
function redistribute(inst::Instance, ws::Workspace, s::Solution, order::Vector{Int},
                      ζ::Int, cmax::Float64, rng)
    N = nnodes(inst)
    inleft = falses(N)
    for k in 1:ζ
        inleft[order[k]] = true
    end
    # các cung cắt ngang và luồng hiện tại
    F = Dict{Tuple{Int,Int},Int}()
    for k in 1:ζ
        i = order[k]
        for (j, f) in s.succ[i]
            inleft[j] || (F[(i, j)] = f)
        end
    end
    F0 = copy(F)
    lefts = unique(first.(collect(keys(F))))           # π^left_i > 0
    rights = unique(last.(collect(keys(F))))           # π^right_j > 0
    Δ(i, j) = ws.L[i] + inst.c[i] + inst.t[i, j] + ws.Lstar[j]

    while true
        # dòng 4: cập nhật Cmax
        cur = maximum(Δ(i, j) for ((i, j), f) in F if f > 0)
        cmax = min(cmax, cur)
        forbidden(i, j) = Δ(i, j) >= cmax - EPS
        # dòng 6–7: chọn ngẫu nhiên một cung găng có luồng
        Amax = [(i, j) for ((i, j), f) in F if f > 0 && forbidden(i, j)]
        isempty(Amax) && break
        (istar, jstar) = Amax[rand(rng, 1:length(Amax))]
        # dòng 8: BFS tìm dãy Γ = (i*, j1, i2, j2, ..., iK, j*)
        #   - cung (i_k, j_k) thêm luồng: không bị cấm
        #   - cung (i_{k+1}, j_k) bớt luồng: f > 0
        parentR = Dict{Int,Int}()      # j -> i đã đi tới j
        parentL = Dict{Int,Int}(istar => 0)   # i -> j mà từ đó đi lùi về i
        queue = [istar]; head = 1; found = 0
        while head <= length(queue) && found == 0
            a = queue[head]; head += 1
            if a != istar && !forbidden(a, jstar)
                found = a; break
            end
            for j in rights
                (j == jstar || haskey(parentR, j) || forbidden(a, j)) && continue
                parentR[j] = a
                for b in lefts
                    if get(F, (b, j), 0) > 0 && !haskey(parentL, b)
                        parentL[b] = j
                        push!(queue, b)
                    end
                end
            end
        end
        found == 0 && break            # dòng 9–10: không tìm được Γ -> dừng
        # dựng lại đường: các cung tăng (+) và giảm (−)
        plus = [(found, jstar)]; minus = [(istar, jstar)]
        b = found
        while b != istar
            j = parentL[b]; a = parentR[j]
            push!(minus, (b, j)); push!(plus, (a, j))
            b = a
        end
        δ = minimum(F[e] for e in minus)
        for e in minus
            F[e] -= δ
        end
        for e in plus
            F[e] = get(F, e, 0) + δ
        end
    end
    # chuyển thành danh sách thay đổi so với nghiệm gốc
    ch = Changes()
    for (e, f) in F
        d = f - get(F0, e, 0)
        d != 0 && push!(ch, (e[1], e[2], d))
    end
    isempty(ch) && return nothing
    sort!(ch; by = x -> x[3])
    newc = maximum(Δ(i, j) for ((i, j), f) in F if f > 0)
    removed = [(a, b) for (a, b, d) in ch if d < 0]
    added = [(a, b) for (a, b, d) in ch if d > 0]
    return Move(:redis, ch, removed, added, ζ), newc
end
