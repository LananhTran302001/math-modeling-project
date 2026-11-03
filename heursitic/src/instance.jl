# ============================================================================
# Dữ liệu bài toán PDP-CR
#
# Quy ước chỉ số (Julia đánh chỉ số từ 1):
#   1..n   : các task
#   S = n+1: depot xuất phát (task 0 trong bài báo)
#   T = n+2: depot kết thúc (task n+1 trong bài báo)
# ============================================================================

struct Instance
    name::String
    n::Int                      # số task
    m::Int                      # số robot
    r::Vector{Int}              # r[i]: số robot cần cho task i (r[S] = r[T] = m)
    c::Vector{Float64}          # c[i]: thời gian xử lý task i (c[S] = c[T] = 0)
    t::Matrix{Float64}          # t[i,j]: thời gian đi từ đích d_i tới điểm nhận p_j
    pickup::Vector{NTuple{2,Float64}}    # tọa độ p_i (có thể rỗng nếu dữ liệu cho sẵn t)
    delivery::Vector{NTuple{2,Float64}}  # tọa độ d_i
end

source(inst::Instance) = inst.n + 1
sink(inst::Instance) = inst.n + 2
nnodes(inst::Instance) = inst.n + 2

euclid(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)

"""
    Instance(name, m, depot, pickup, delivery, r; c = nothing)

Xây dựng instance từ tọa độ. Thời gian xử lý mặc định là khoảng cách Euclid p_i -> d_i
(như Mục 5.1 của bài báo); thời gian di chuyển t[i,j] = ||d_i - p_j||.
"""
function Instance(name::AbstractString, m::Int, depot::NTuple{2,Float64},
                  pickup::Vector{NTuple{2,Float64}}, delivery::Vector{NTuple{2,Float64}},
                  r::Vector{Int}; c = nothing)
    n = length(pickup)
    @assert length(delivery) == n && length(r) == n
    @assert all(1 .<= r .<= m) "mỗi task phải cần từ 1 tới m robot"
    P = vcat(pickup, [depot, depot])      # p của S và T là depot
    D = vcat(delivery, [depot, depot])
    cc = c === nothing ? [euclid(pickup[i], delivery[i]) for i in 1:n] : collect(Float64, c)
    cc = vcat(cc, [0.0, 0.0])
    rr = vcat(r, [m, m])
    N = n + 2
    t = [euclid(D[i], P[j]) for i in 1:N, j in 1:N]
    t[n+1, n+2] = 0.0                      # robot không dùng: đi thẳng S -> T
    return Instance(String(name), n, m, rr, cc, t, P, D)
end

"""
    Instance(name, m, r, c, t)

Xây dựng instance khi đã có sẵn ma trận thời gian di chuyển `t` cỡ (n+2)x(n+2)
(theo quy ước chỉ số ở trên). Dùng cho Ví dụ 1 trong bài báo.
"""
function Instance(name::AbstractString, m::Int, r::Vector{Int}, c::Vector{Float64},
                  t::Matrix{Float64})
    n = length(r)
    @assert size(t) == (n + 2, n + 2)
    return Instance(String(name), n, m, vcat(r, [m, m]), vcat(c, [0.0, 0.0]), t,
                    NTuple{2,Float64}[], NTuple{2,Float64}[])
end

# ----------------------------------------------------------------------------
# Ví dụ 1 trong bài báo (n = 3, m = 2), dùng nguyên bảng thời gian Table 2.
# Nghiệm trong Hình 3 có makespan 21.5.
# ----------------------------------------------------------------------------
function example1()
    INF = 1e6  # cặp không xuất hiện trong bảng (i == j)
    # thứ tự chỉ số: 1,2,3 = task; 4 = S (task 0); 5 = T (task 4)
    t = fill(INF, 5, 5)
    tab = Dict((0, 1) => 6, (0, 2) => 7, (0, 3) => 10, (0, 4) => 0,
               (1, 2) => 2, (1, 3) => 7, (1, 4) => 3,
               (2, 1) => 4, (2, 3) => 6, (2, 4) => 6,
               (3, 1) => 3, (3, 2) => 4, (3, 4) => 9)
    idx(k) = k == 0 ? 4 : (k == 4 ? 5 : k)
    for ((a, b), v) in tab
        t[idx(a), idx(b)] = v
    end
    return Instance("example1", 2, [2, 1, 1], [3.0, 2.5, 2.5], t)
end

# ----------------------------------------------------------------------------
# Đọc benchmark Li & Lim (PDPTW) và chuyển sang PDP-CR theo Mục 5.1:
#   - mỗi cặp (pickup, delivery) thành 1 task
#   - c_i = khoảng cách Euclid pickup -> delivery
#   - r_i ~ Uniform{1..rmax} (sinh ngẫu nhiên với seed cố định để cả nhóm dùng chung)
#   - bỏ qua cửa sổ thời gian, tải trọng
# ----------------------------------------------------------------------------
function read_lilim(path::AbstractString; m::Int = 10, rmax::Int = 5, seed::Int = 1,
                    ntasks::Union{Nothing,Int} = nothing)
    lines = filter(!isempty, strip.(readlines(path)))
    rows = [parse.(Float64, split(l)) for l in lines[2:end]]
    depot = (rows[1][2], rows[1][3])
    byid = Dict(Int(row[1]) => row for row in rows)
    pickup = NTuple{2,Float64}[]
    delivery = NTuple{2,Float64}[]
    for row in rows[2:end]
        Int(row[8]) == 0 || continue           # chỉ xét dòng pickup (trường pickup-index = 0)
        dl = byid[Int(row[9])]
        push!(pickup, (row[2], row[3]))
        push!(delivery, (dl[2], dl[3]))
    end
    if ntasks !== nothing
        pickup, delivery = pickup[1:ntasks], delivery[1:ntasks]
    end
    rng = Random.Xoshiro(seed)
    rmax = min(rmax, m)
    r = rand(rng, 1:rmax, length(pickup))
    name = splitext(basename(path))[1]
    return Instance(name, m, depot, pickup, delivery, r)
end

"""
    random_instance(n, m; seed, side = 100.0, rmax = 5)

Sinh instance ngẫu nhiên kiểu Li & Lim (tọa độ đều trên hình vuông `side` x `side`,
depot ở tâm). Dùng khi cần cỡ lớn (200-500 task) hoặc cỡ nhỏ để so với MILP.
"""
function random_instance(n::Int, m::Int; seed::Int = 1, side::Float64 = 100.0, rmax::Int = 5)
    rng = Random.Xoshiro(seed)
    pt() = (round(side * rand(rng); digits = 1), round(side * rand(rng); digits = 1))
    pickup = [pt() for _ in 1:n]
    delivery = [pt() for _ in 1:n]
    r = rand(rng, 1:min(rmax, m), n)
    return Instance("rand_n$(n)_m$(m)_s$(seed)", m, (side / 2, side / 2), pickup, delivery, r)
end

"""
    write_instance(path, inst)

Ghi instance ra file văn bản để người làm MILP đọc lại đúng cùng dữ liệu.
Định dạng:
    dòng 1: n m depot_x depot_y
    n dòng: i r_i c_i px py dx dy          (task 1..n)
    sau đó: ma trận t cỡ (n+2)x(n+2), chỉ số 1..n task, n+1 = depot đầu, n+2 = depot cuối
"""
function write_instance(path::AbstractString, inst::Instance)
    open(path, "w") do io
        hasxy = !isempty(inst.pickup)
        dep = hasxy ? inst.pickup[source(inst)] : (NaN, NaN)
        println(io, join((inst.n, inst.m, dep[1], dep[2]), " "))
        for i in 1:inst.n
            p = hasxy ? inst.pickup[i] : (NaN, NaN)
            d = hasxy ? inst.delivery[i] : (NaN, NaN)
            println(io, join((i, inst.r[i], inst.c[i], p[1], p[2], d[1], d[2]), " "))
        end
        for i in 1:nnodes(inst)
            println(io, join(inst.t[i, :], " "))
        end
    end
    return path
end

function read_instance(path::AbstractString)
    lines = readlines(path)
    h = split(lines[1])
    n, m = parse(Int, h[1]), parse(Int, h[2])
    dep = (parse(Float64, h[3]), parse(Float64, h[4]))
    r = Int[]; c = Float64[]; P = NTuple{2,Float64}[]; D = NTuple{2,Float64}[]
    for k in 1:n
        v = split(lines[1+k])
        push!(r, parse(Int, v[2])); push!(c, parse(Float64, v[3]))
        push!(P, (parse(Float64, v[4]), parse(Float64, v[5])))
        push!(D, (parse(Float64, v[6]), parse(Float64, v[7])))
    end
    t = reduce(vcat, [permutedims(parse.(Float64, split(lines[1+n+k]))) for k in 1:n+2])
    name = splitext(basename(path))[1]
    isnan(dep[1]) && return Instance(name, m, r, c, t)
    return Instance(name, n, m, vcat(r, [m, m]), vcat(c, [0.0, 0.0]), t,
                    vcat(P, [dep, dep]), vcat(D, [dep, dep]))
end
