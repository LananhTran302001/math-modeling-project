# ============================================================================
# Biểu diễn nghiệm: đồ thị luồng F (resource-flow graph, Poppenborg & Knust)
#
#   f[i][j] > 0  <=>  có f_ij robot làm task j ngay sau task i
#   Ràng buộc:  sum_j f_ij = r_i,  sum_i f_ij = r_j  (bảo toàn luồng),
#               f_ij nguyên >= 0, F không có chu trình.
#
# Lưu dạng danh sách kề thưa (≤ ~5n cung) để copy và đánh giá nhanh.
# ============================================================================

mutable struct Solution
    succ::Vector{Dict{Int,Int}}   # succ[i][j] = f_ij
    pred::Vector{Dict{Int,Int}}   # pred[j][i] = f_ij
end

Solution(N::Int) = Solution([Dict{Int,Int}() for _ in 1:N], [Dict{Int,Int}() for _ in 1:N])

Base.copy(s::Solution) = Solution(copy.(s.succ), copy.(s.pred))

flow(s::Solution, i, j) = get(s.succ[i], j, 0)

"Cộng `d` đơn vị luồng vào cung (i,j); xóa cung khi luồng về 0."
function add_flow!(s::Solution, i::Int, j::Int, d::Int)
    v = get(s.succ[i], j, 0) + d
    v < 0 && error("luồng âm trên cung ($i,$j)")
    if v == 0
        delete!(s.succ[i], j); delete!(s.pred[j], i)
    else
        s.succ[i][j] = v; s.pred[j][i] = v
    end
    return s
end

"Một bước biến đổi luồng: danh sách (i, j, Δf)."
const Changes = Vector{Tuple{Int,Int,Int}}

apply!(s::Solution, ch::Changes) = (for (i, j, d) in ch; add_flow!(s, i, j, d); end; s)
undo!(s::Solution, ch::Changes) = (for k in length(ch):-1:1; (i, j, d) = ch[k]; add_flow!(s, i, j, -d); end; s)

arcs(s::Solution) = [(i, j, f) for i in eachindex(s.succ) for (j, f) in s.succ[i]]

# ----------------------------------------------------------------------------
# Workspace để tính đường dài nhất mà không cấp phát lại bộ nhớ
# ----------------------------------------------------------------------------
struct Workspace
    indeg::Vector{Int}
    order::Vector{Int}
    L::Vector{Float64}
    Lstar::Vector{Float64}
end
Workspace(N::Int) = Workspace(zeros(Int, N), zeros(Int, N), zeros(N), zeros(N))

"""
    topo_order!(ws, inst, s) -> Bool

Thứ tự topo (Kahn). Trả về `false` nếu F có chu trình (nghiệm không khả thi).
"""
function topo_order!(ws::Workspace, inst::Instance, s::Solution)
    N = nnodes(inst)
    indeg, order = ws.indeg, ws.order
    for v in 1:N
        indeg[v] = length(s.pred[v])
    end
    head = 0; tail = 0
    for v in 1:N
        if indeg[v] == 0
            tail += 1; order[tail] = v
        end
    end
    while head < tail
        head += 1
        u = order[head]
        for v in keys(s.succ[u])
            indeg[v] -= 1
            if indeg[v] == 0
                tail += 1; order[tail] = v
            end
        end
    end
    return tail == N
end

"""
    makespan!(ws, inst, s) -> Float64

Thuật toán 1: L_j = max_{i pred j} (L_i + c_i + t_ij), L_S = 0. Makespan = L_T.
Trả về `Inf` nếu đồ thị có chu trình.
"""
function makespan!(ws::Workspace, inst::Instance, s::Solution)
    topo_order!(ws, inst, s) || return Inf
    L, c, t = ws.L, inst.c, inst.t
    fill!(L, 0.0)
    @inbounds for u in ws.order
        Lu = L[u] + c[u]
        for v in keys(s.succ[u])
            x = Lu + t[u, v]
            x > L[v] && (L[v] = x)
        end
    end
    return L[sink(inst)]
end

"""
    tails!(ws, inst, s)

Thuật toán 2: L*_i = max_{j succ i} (c_i + t_ij + L*_j), L*_T = 0
(độ dài đường dài nhất từ điểm nhận của i tới depot cuối). Cần gọi `makespan!` trước.
"""
function tails!(ws::Workspace, inst::Instance, s::Solution)
    Ls, c, t = ws.Lstar, inst.c, inst.t
    fill!(Ls, 0.0)
    @inbounds for k in nnodes(inst):-1:1
        u = ws.order[k]
        best = 0.0
        for v in keys(s.succ[u])
            x = c[u] + t[u, v] + Ls[v]
            x > best && (best = x)
        end
        Ls[u] = best
    end
    return Ls
end

const EPS = 1e-6

"""
    critical_arcs(ws, inst, s, cmax)

Các cung (i,j) có f_ij > 0 nằm trên một đường găng: L_i + c_i + t_ij + L*_j = Cmax.
"""
function critical_arcs(ws::Workspace, inst::Instance, s::Solution, cmax::Float64)
    res = Tuple{Int,Int}[]
    for i in 1:nnodes(inst), j in keys(s.succ[i])
        if ws.L[i] + inst.c[i] + inst.t[i, j] + ws.Lstar[j] >= cmax - EPS
            push!(res, (i, j))
        end
    end
    return res
end

"""
    topo_by_start(ws, inst)

Thứ tự topo dùng cho phép cắt ζ: sắp theo thời điểm bắt đầu sớm nhất L_i,
hòa thì theo vị trí trong thứ tự Kahn (luôn là thứ tự topo hợp lệ vì
f_ij > 0 ⇒ L_j ≥ L_i + c_i + t_ij ≥ L_i).
"""
function topo_by_start(ws::Workspace, inst::Instance)
    N = nnodes(inst)
    pos = zeros(Int, N)
    for k in 1:N
        pos[ws.order[k]] = k
    end
    return sort(collect(1:N); by = v -> (ws.L[v], pos[v]))
end

# ----------------------------------------------------------------------------
# Kiểm tra tính khả thi đầy đủ (dùng trong test và sau khi chạy heuristic)
# ----------------------------------------------------------------------------
function check_feasible(inst::Instance, s::Solution; verbose = false)
    S, T, N = source(inst), sink(inst), nnodes(inst)
    msgs = String[]
    for i in 1:N, (j, f) in s.succ[i]
        f > 0 || push!(msgs, "cung ($i,$j) có luồng $f ≤ 0")
        get(s.pred[j], i, -1) == f || push!(msgs, "succ/pred lệch tại ($i,$j)")
        i == T && push!(msgs, "có cung đi ra từ depot cuối")
        j == S && push!(msgs, "có cung đi vào depot đầu")
        i == j && push!(msgs, "khuyên tại $i")
    end
    for v in 1:N
        out = sum(values(s.succ[v]); init = 0)
        inn = sum(values(s.pred[v]); init = 0)
        v != T && out != inst.r[v] && push!(msgs, "node $v: luồng ra $out ≠ r = $(inst.r[v])")
        v != S && inn != inst.r[v] && push!(msgs, "node $v: luồng vào $inn ≠ r = $(inst.r[v])")
    end
    ws = Workspace(N)
    topo_order!(ws, inst, s) || push!(msgs, "đồ thị có chu trình")
    verbose && foreach(println, msgs)
    return isempty(msgs)
end

# ----------------------------------------------------------------------------
# Phân rã luồng thành lộ trình từng robot + lịch (để vẽ Gantt / báo cáo)
# ----------------------------------------------------------------------------
"""
    robot_routes(inst, s) -> Vector{Vector{Int}}

Phân rã luồng: mỗi robot đi một đường S -> ... -> T. Robot đến task j sớm sẽ chờ
tới thời điểm bắt đầu L_j.
"""
function robot_routes(inst::Instance, s::Solution)
    S, T = source(inst), sink(inst)
    rem = copy(s)
    routes = Vector{Vector{Int}}()
    ws = Workspace(nnodes(inst)); makespan!(ws, inst, s)
    for _ in 1:inst.m
        route = [S]; u = S
        while u != T
            # đi theo cung có thời điểm đến sớm nhất (phân rã "tự nhiên" theo thời gian)
            v = argmin(j -> ws.L[j], collect(keys(rem.succ[u])))
            add_flow!(rem, u, v, -1)
            push!(route, v); u = v
        end
        push!(routes, route)
    end
    return routes
end

function print_schedule(io::IO, inst::Instance, s::Solution)
    ws = Workspace(nnodes(inst))
    cmax = makespan!(ws, inst, s)
    S, T = source(inst), sink(inst)
    name(v) = v == S ? "S" : v == T ? "T" : string(v)
    println(io, "Makespan = ", round(cmax; digits = 3))
    for (k, route) in enumerate(robot_routes(inst, s))
        parts = String[]
        for (a, b) in zip(route[1:end-1], route[2:end])
            arrive = ws.L[a] + inst.c[a] + inst.t[a, b]
            start = b == T ? arrive : ws.L[b]
            wait = start - arrive
            push!(parts, "→ $(name(b))[đến $(round(arrive; digits=2))" *
                         (wait > EPS ? ", chờ $(round(wait; digits=2))" : "") *
                         (b == T ? "" : ", làm $(round(start; digits=2))–$(round(start + inst.c[b]; digits=2))") * "]")
        end
        println(io, "  Robot $k: S ", join(parts, " "))
    end
end
print_schedule(inst::Instance, s::Solution) = print_schedule(stdout, inst, s)
