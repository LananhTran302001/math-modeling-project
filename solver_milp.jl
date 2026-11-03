using JuMP
using HiGHS

# ============================================================
# 1. ĐỌC INSTANCE 
# ============================================================

function read_instance_txt(file_path::String)
    lines = readlines(file_path)

    # Loại bỏ các dòng trống
    lines = filter(line -> !isempty(strip(line)), lines)

    isempty(lines) && error("File đầu vào rỗng: $file_path")

    # Dòng đầu: n m depot_x depot_y
    header = split(strip(lines[1]))

    if length(header) < 4
        error(
            "Header không hợp lệ. Nội dung thực tế: $(repr(lines[1])). " *
            "Định dạng yêu cầu: n m depot_x depot_y"
        )
    end

    n = parse(Int, header[1])
    m = parse(Int, header[2])
    depot_x = parse(Float64, header[3])
    depot_y = parse(Float64, header[4])

    n > 0 || error("Số task n phải lớn hơn 0.")
    m > 0 || error("Số robot m phải lớn hơn 0.")

    # Kiểm tra số dòng tối thiểu:
    # 1 dòng header + n dòng task + (n+2) dòng ma trận
    expected_lines = 1 + n + (n + 2)

    length(lines) == expected_lines ||
        error(
            "Số dòng không hợp lệ trong $file_path. " *
            "Mong đợi $expected_lines dòng, nhận được $(length(lines))."
        )

    # Đọc thông tin task
    # Mỗi dòng: task_id r_i c_i pickup_x pickup_y delivery_x delivery_y

    r = Dict{Int, Int}()
    c = Dict{Int, Float64}()

    pickup = Dict{Int, Tuple{Float64, Float64}}()
    delivery = Dict{Int, Tuple{Float64, Float64}}()

    for k in 1:n
        row = split(strip(lines[1 + k]))

        length(row) == 7 ||
            error("Dòng task thứ $k phải có đúng 7 giá trị.")

        task_id = parse(Int, row[1])
        r_i = parse(Int, row[2])
        c_i = parse(Float64, row[3])

        pickup_x = parse(Float64, row[4])
        pickup_y = parse(Float64, row[5])

        delivery_x = parse(Float64, row[6])
        delivery_y = parse(Float64, row[7])

        1 <= task_id <= n ||
            error("Task ID không hợp lệ: $task_id")

        1 <= r_i <= m ||
            error("r[$task_id] = $r_i không nằm trong khoảng 1..$m")

        c_i >= 0 ||
            error("Thời gian xử lý c[$task_id] phải không âm.")

        r[task_id] = r_i
        c[task_id] = c_i

        pickup[task_id] = (pickup_x, pickup_y)
        delivery[task_id] = (delivery_x, delivery_y)
    end

    # Depot xuất phát: 0
    # Tasks: 1..n
    # Depot kết thúc: n+1

    S = 0
    T = n + 1

    # Thời gian xử lý tại hai depot bằng 0
    c[S] = 0.0
    c[T] = 0.0

    # Đọc ma trận thời gian di chuyển t_ij
    # Ma trận có (n+2) hàng và (n+2) cột

    t_matrix = zeros(Float64, n + 2, n + 2)

    matrix_start = n + 2

    for i in 1:(n + 2)
        row = split(strip(lines[matrix_start + i - 1]))

        length(row) == n + 2 ||
            error(
                "Hàng ma trận $i phải có $(n+2) giá trị, " *
                "nhận được $(length(row))."
            )

        for j in 1:(n + 2)
            t_matrix[i, j] = parse(Float64, row[j])

            t_matrix[i, j] >= 0 ||
                error("t[$i,$j] phải không âm.")
        end
    end

    println("Đã đọc instance: $(basename(file_path))")
    println("  Số task: $n")
    println("  Số robot: $m")
    println("  Depot: ($depot_x, $depot_y)")
    println("  Kích thước ma trận t: $(size(t_matrix))")

    return n, m, r, c, t_matrix
end


# ============================================================
# 2. GIẢI MÔ HÌNH FLOW-BASED MILP
# ============================================================

function solve_flowbased_milp(txt_file_path::String)

    # Đọc dữ liệu TXT
    n, m, r, c, t_matrix = read_instance_txt(txt_file_path)

    # Tập chỉ số:
    # S = 0       : depot xuất phát
    # N = 1..n    : các task
    # T = n+1     : depot kết thúc

    S = 0
    T = n + 1

    N = 1:n
    N_plus = S:n
    N_minus = 1:T
    N_all = S:T

    # Big-M
    max_t = maximum(t_matrix)

    M = sum(values(c)) + (n + 1) * max_t + 1000.0

    # Khởi tạo mô hình
    model = Model(HiGHS.Optimizer)

    set_attribute(model, "time_limit", 3600.0)
    set_attribute(model, "mip_rel_gap", 0.05)

    # --------------------------------------------------------
    # 2.1. Biến quyết định
    # --------------------------------------------------------

    # x[i,j] = 1 nếu có tuyến di chuyển từ i đến j
    @variable(
        model,
        x[i in N_plus, j in N_minus; i != j],
        Bin
    )

    # f[i,j]: số robot di chuyển từ i đến j
    @variable(
        model,
        f[i in N_plus, j in N_minus; i != j] >= 0
    )

    # z[i]: thời điểm bắt đầu task i
    @variable(model, z[i in N_all] >= 0)

    # --------------------------------------------------------
    # 2.2. Hàm mục tiêu
    # --------------------------------------------------------

    # Minimize makespan
    @objective(model, Min, z[T])

    # --------------------------------------------------------
    # 2.3. Các ràng buộc
    # --------------------------------------------------------

    # (1) Tổng luồng đi ra khỏi depot xuất phát bằng m
    @constraint(
        model,
        sum(f[S, j] for j in N_minus) == m
    )

    # (2) Tổng luồng đi vào depot kết thúc bằng m
    @constraint(
        model,
        sum(f[i, T] for i in N_plus) == m
    )

    # (3) và (4) Bảo toàn luồng tại từng task
    for j in N
        @constraint(
            model,
            sum(f[i, j] for i in N_plus if i != j) == r[j]
        )

        @constraint(
            model,
            sum(f[j, k] for k in N_minus if k != j) == r[j]
        )
    end

    # (5) Thời điểm bắt đầu tại depot xuất phát
    @constraint(model, z[S] == 0)

    # (6) Ràng buộc thời gian và loại bỏ subtour
    #
    # z[j] >= z[i] + c[i] + t[i,j] - M(1-x[i,j])
    #
    # t_matrix được đánh chỉ số từ 1 trong Julia,
    # nên phải truy cập t_matrix[i+1, j+1].

    for i in N_plus, j in N_minus
        if i != j
            t_ij = t_matrix[i + 1, j + 1]

            @constraint(
                model,
                z[j] + M * (1 - x[i, j]) >=
                z[i] + c[i] + t_ij
            )
        end
    end

    # (7) và (8) Liên kết luồng f với tuyến x
    for i in N_plus, j in N_minus
        if i != j

            # Có tuyến thì phải có ít nhất 1 robot
            @constraint(model, f[i, j] >= x[i, j])

            # Giới hạn luồng trên mỗi tuyến
            cap = (i in N) ? r[i] : m

            @constraint(
                model,
                f[i, j] <= cap * x[i, j]
            )
        end
    end

    # --------------------------------------------------------
    # 2.4. Giải mô hình
    # --------------------------------------------------------

    println("\nĐang giải mô hình MILP...")
    optimize!(model)

    # --------------------------------------------------------
    # 2.5. In kết quả
    # --------------------------------------------------------

    status = termination_status(model)

    println("\n================ KẾT QUẢ ================")
    println("Instance: $txt_file_path")
    println("Trạng thái: $status")

    if !has_values(model)
        println("Không có nghiệm khả thi để hiển thị.")
        return
    end

    makespan = objective_value(model)

    println("Makespan: ", round(makespan, digits=4))

    println("\nThời điểm bắt đầu:")
    for i in N_all
        println(
            "  z[$i] = ",
            round(value(z[i]), digits=4)
        )
    end

    println("\nCác tuyến có luồng robot dương:")
    for i in N_plus, j in N_minus
        if i != j
            flow_value = value(f[i, j])

            if flow_value > 1e-3
                println(
                    "  $i -> $j: ",
                    round(flow_value, digits=4),
                    " robot"
                )
            end
        end
    end

    println("\nGiá trị mục tiêu: ", objective_value(model))
    println("Trạng thái kết thúc: $status")

end



# ============================================================
# 3. CHẠY 
# ============================================================

input_dir = ""
output_dir = ""

mkpath(output_dir)

txt_files = filter(
    f -> endswith(lowercase(f), ".txt"),
    readdir(input_dir)
)

for file in sort(txt_files)

    txt_path = joinpath(input_dir, file)

    filename = splitext(file)[1]
    output_path = joinpath(
        output_dir,
        "$(filename)_result.txt"
    )

    println("Đang xử lý: $file")

    try
        open(output_path, "w") do io
            redirect_stdout(io) do
                solve_flowbased_milp(txt_path)
            end
        end

        println("Đã lưu kết quả: $output_path")

    catch e
        println("Lỗi khi xử lý $file: ", sprint(showerror, e))
    end
end
