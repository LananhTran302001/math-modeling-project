# ============================================================================
# Tabu search (Mục 4 bài báo) — hai biến thể:
#   TSO: chỉ reroute + reverse              (Table 4, dòng TSO)
#   TSR: reroute + reverse + redistribute   (Table 4, dòng TSR)
# ============================================================================

Base.@kwdef struct TabuParams
    tl_add::Int                 # kích thước TL_add
    tl_redis::Int               # kích thước TL_redis
    n_reroute::Int              # số ứng viên reroute mỗi vòng
    n_reverse::Int              # số ứng viên reverse mỗi vòng
    n_redis::Int                # số lát cắt ζ thử mỗi vòng
    β::Float64 = 0.5            # xác suất nhận ứng viên redis khi không cải thiện
    max_iter::Int = 5000
    time_limit::Float64 = Inf   # giây
    restart_best::Int = 300     # không cải thiện best sau bấy nhiêu vòng -> quay về best
    restart_new::Int = 600      # ... -> khởi động lại từ nghiệm tham lam ngẫu nhiên
    aspiration::Bool = true     # cho phép nước tabu nếu nó tốt hơn best
    rcl::Int = 3                # độ ngẫu nhiên của greedy khi restart
    reroute_filter::Bool = true # lọc (u,v) theo thời gian để tránh tạo chu trình (xem reroute_moves)
    reroute_all_crit::Bool = true  # thiếu ứng viên thì dùng thêm các cung găng khác
end

"Tham số theo Table 4 (làm tròn xuống, tối thiểu 1 để chạy được với n nhỏ)."
function TabuParams(variant::Symbol, n::Int; kw...)
    f(x) = max(1, floor(Int, x))
    if variant == :TSO
        return TabuParams(; tl_add = f(n / 8), tl_redis = 0, n_reroute = f(n / 3),
                          n_reverse = f(n / 9), n_redis = 0, kw...)
    elseif variant == :TSR
        return TabuParams(; tl_add = f(n / 8), tl_redis = f(n / 15), n_reroute = f(n / 5),
                          n_reverse = f(n / 15), n_redis = f(n / 15), kw...)
    else
        error("variant phải là :TSO hoặc :TSR")
    end
end

"Danh sách tabu FIFO có giới hạn kích thước."
struct TabuList{T}
    cap::Int
    items::Vector{T}
end
TabuList{T}(cap::Int) where {T} = TabuList{T}(cap, T[])
Base.in(x, tl::TabuList) = x in tl.items
function Base.push!(tl::TabuList, x)
    tl.cap <= 0 && return tl
    push!(tl.items, x)
    length(tl.items) > tl.cap && popfirst!(tl.items)
    return tl
end

"""
Quy ước tabu (thống nhất trong nhóm):
- reroute/reverse/redis ghi các cung bị *giảm* luồng vào TL_add;
- reroute tabu nếu một cung nó *tạo thêm* ((i,v) hoặc (u,j)) nằm trong TL_add
  (bài báo viết (i',u'), (u',j') — chúng tôi hiểu là các cung được thêm);
- reverse trên (i,j) tabu nếu (j,i) ∈ TL_add (đúng như bài báo);
- redis tại ζ tabu nếu ζ ∈ TL_redis.
"""
function is_tabu(mv::Move, tladd::TabuList, tlredis::TabuList)
    if mv.kind == :redis
        return mv.cut in tlredis
    elseif mv.kind == :reverse
        (i, j) = mv.removed[1]
        return (j, i) in tladd
    else
        return any(a -> a in tladd, mv.added)
    end
end

function update_tabu!(mv::Move, tladd::TabuList, tlredis::TabuList)
    if mv.kind == :redis
        push!(tlredis, mv.cut)
    else
        foreach(a -> push!(tladd, a), mv.removed)
    end
end

struct TabuResult
    best::Solution
    best_value::Float64
    initial_value::Float64
    history::Vector{Float64}       # best sau mỗi vòng (để vẽ như Hình 12)
    iterations::Int
    time::Float64
    accepted::Dict{Symbol,Int}     # số nước được chọn theo loại
    improved_by::Dict{Symbol,Int}  # số lần cải thiện best theo loại
end

"""
    tabu_search(inst, params; init = nothing, rng, verbose = false)

Tabu search trên đồ thị luồng. `init` là nghiệm ban đầu (mặc định: greedy).
"""
function tabu_search(inst::Instance, P::TabuParams; init::Union{Nothing,Solution} = nothing,
                     rng = Random.Xoshiro(1), verbose::Bool = false, log_every::Int = 500)
    t0 = time()
    ws = Workspace(nnodes(inst))
    ws2 = Workspace(nnodes(inst))
    cur = init === nothing ? greedy_solution(inst) : copy(init)
    curval = makespan!(ws, inst, cur)
    isfinite(curval) || error("nghiệm ban đầu không khả thi")
    best, bestval = copy(cur), curval
    initval = curval
    tladd = TabuList{Tuple{Int,Int}}(P.tl_add)
    tlredis = TabuList{Int}(P.tl_redis)
    history = Float64[]
    accepted = Dict(:reroute => 0, :reverse => 0, :redis => 0)
    improved_by = Dict(:reroute => 0, :reverse => 0, :redis => 0)
    since_best = 0
    nnode = nnodes(inst)
    iter = 0
    while iter < P.max_iter && time() - t0 < P.time_limit
        iter += 1
        # --- thông tin của nghiệm hiện tại -----------------------------------
        makespan!(ws, inst, cur)
        tails!(ws, inst, cur)
        crit = critical_arcs(ws, inst, cur, curval)

        # --- sinh và đánh giá ứng viên ---------------------------------------
        best_other = nothing; val_other = Inf
        best_redis = nothing; val_redis = Inf
        consider(mv, v) = (is_tabu(mv, tladd, tlredis) && !(P.aspiration && v < bestval - EPS))

        for mv in Iterators.flatten((reroute_moves(inst, cur, crit, P.n_reroute, rng;
                                                   L = P.reroute_filter ? ws.L : nothing,
                                                   all_crit = P.reroute_all_crit),
                                     reverse_moves(inst, ws, cur, crit, P.n_reverse, rng)))
            apply!(cur, mv.changes)
            v = makespan!(ws2, inst, cur)          # workspace riêng, giữ L, L* của cur
            undo!(cur, mv.changes)
            consider(mv, v) && continue
            if v < val_other
                best_other, val_other = mv, v
            end
        end
        if P.n_redis > 0
            order = topo_by_start(ws, inst)
            cuts = shuffle(rng, 1:nnode-1)
            tried = 0
            for ζ in cuts
                tried >= P.n_redis && break
                (ζ in tlredis) && continue
                tried += 1
                res = redistribute(inst, ws, cur, order, ζ, curval, rng)
                res === nothing && continue
                mv, v = res
                if v < val_redis
                    best_redis, val_redis = mv, v
                end
            end
        end

        # --- chọn nước đi ----------------------------------------------------
        chosen = nothing; chosen_val = Inf
        if min(val_other, val_redis) < curval - EPS
            # có ứng viên cải thiện: nhận ngay ứng viên tốt nhất
            chosen, chosen_val = val_redis <= val_other ? (best_redis, val_redis) : (best_other, val_other)
        elseif val_redis <= val_other && best_redis !== nothing
            # tốt nhất đến từ redis: nhận với xác suất β, ngược lại nhận reroute/reverse
            if best_other === nothing || rand(rng) < P.β
                chosen, chosen_val = best_redis, val_redis
            else
                chosen, chosen_val = best_other, val_other
            end
        elseif best_other !== nothing
            chosen, chosen_val = best_other, val_other
        end

        if chosen !== nothing && isfinite(chosen_val)
            apply!(cur, chosen.changes)
            curval = chosen_val
            update_tabu!(chosen, tladd, tlredis)
            accepted[chosen.kind] += 1
        end

        # --- cập nhật best và cơ chế restart ---------------------------------
        if curval < bestval - EPS
            best, bestval = copy(cur), curval
            since_best = 0
            chosen === nothing || (improved_by[chosen.kind] += 1)
        else
            since_best += 1
            if since_best == P.restart_best
                cur, curval = copy(best), bestval
            elseif since_best >= P.restart_new
                cur = greedy_solution(inst; rng = rng, rcl = P.rcl)
                curval = makespan!(ws, inst, cur)
                since_best = 0
                if curval < bestval - EPS
                    best, bestval = copy(cur), curval
                end
            end
        end
        push!(history, bestval)
        if verbose && iter % log_every == 0
            @printf("  iter %5d  cur = %10.3f  best = %10.3f  (%.1fs)\n", iter, curval, bestval, time() - t0)
        end
    end
    return TabuResult(best, bestval, initval, history, iter, time() - t0, accepted, improved_by)
end
