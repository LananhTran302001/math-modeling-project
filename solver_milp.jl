using JuMP
using HiGHS
using JSON

function solve_flowbased_milp(json_file_path::String)
    # 1. Đọc dữ liệu JSON
    data = JSON.parsefile(json_file_path)
    
    m = data["metadata"]["m_robots"]       # Số lượng robot m
    n = data["metadata"]["n_tasks"]        # Số lượng nhiệm vụ n
    tasks_data = data["tasks"]             # Bảng 1: Thông tin các nhiệm vụ
    t_matrix = data["traveling_time_matrix_t_ij"] # Bảng 2: Ma trận t_ij
    
    # Tập chỉ số:
    # N: 1..n
    # N_plus: 0..n
    # N_minus: 1..(n+1)
    # N_all: 0..(n+1)
    N = 1:n
    N_plus = 0:n
    N_minus = 1:(n+1)
    N_all = 0:(n+1)
    
    # Tạo Dictionary lưu c_i và r_i
    c = Dict(task["i"] => Float64(task["c_i"]) for task in tasks_data)
    r = Dict(task["i"] => Int(task["r_i"]) for task in tasks_data)
    
    # Tính Big-M động an toàn dựa trên tổng thời gian c_i và t_ij tối đa 
    max_t = maximum([maximum(row) for row in t_matrix]) 
    M = sum(values(c)) + (n + 1) * max_t + 1000.0
    
    # 2. Khởi tạo mô hình JuMP với bộ giải HiGHS
    model = Model(HiGHS.Optimizer)
    
    # 3. Khai báo Biến quyết định
    # x[i, j] nhị phân: bằng 1 nếu có robot di chuyển trực tiếp từ i sang j
    @variable(model, x[i in N_plus, j in N_minus; i != j], Bin)
    
    # f[i, j] liên tục >= 0: số lượng robot (luồng) di chuyển từ i sang j
    @variable(model, f[i in N_plus, j in N_minus; i != j] >= 0)
    
    # z[i] liên tục >= 0: thời điểm bắt đầu thực hiện nhiệm vụ i
    @variable(model, z[i in N_all] >= 0)
    
    # 4. Hàm Mục tiêu (1): Tối thiểu hóa tổng thời gian hoàn thành (Makespan = z_{n+1})
    @objective(model, Min, z[n+1])
    
    # 5. Khai báo các Ràng buộc
    
    # (2) Luồng xuất phát từ Depot 0 bằng m robot
    @constraint(model, sum(f[0, j] for j in N_minus) == m)
    
    # (3) Luồng trở về Depot n+1 bằng m robot
    @constraint(model, sum(f[i, n+1] for i in N_plus) == m)
    
    # (4) & (5) Bảo toàn luồng tại từng nhiệm vụ j
    for j in N
        @constraint(model, sum(f[i, j] for i in N_plus if i != j) == r[j])
        @constraint(model, sum(f[j, k] for k in N_minus if k != j) == r[j])
    end
    
    # (6) Thời điểm xuất phát tại Depot 0 bằng 0
    @constraint(model, z[0] == 0)
    
    # (7) Ràng buộc thứ tự thời gian & Loại bỏ chu trình (Subtour Elimination)
    # Lưu ý: Ma trận t_matrix trong Julia là 1-indexed nên truy cập bằng [i+1][j+1]
    for i in N_plus, j in N_minus
        if i != j
            t_ij = Float64(t_matrix[i + 1][j + 1])
            @constraint(model, z[j] + M * (1 - x[i, j]) >= z[i] + c[i] + t_ij)
        end
    end
    
    # (8) & (9) Liên kết giữa luồng f[i, j] và tuyến x[i, j]
    for i in N_plus, j in N_minus
        if i != j
            @constraint(model, f[i, j] >= x[i, j])
            cap = (i in N) ? r[i] : m
            @constraint(model, cap * x[i, j] >= f[i, j])
        end
    end
    
    # 6. Giải mô hình
    println("⏳ Đang giải mô hình MILP bằng JuMP & HiGHS...")
    optimize!(model)
    
    # 7. Hiển thị kết quả
    status = termination_status(model)
    println("\n================ KẾT QUẢ TỐI ƯU ($status) ================")
    if status == MOI.OPTIMAL
        makespan = objective_value(model)
        println("🎯 Tổng thời gian hoàn thành (Makespan): ", round(makespan, digits=2))
        
        println("\n Lịch trình bắt đầu từng nhiệm vụ (z_i):")
        for i in N_all
            println("   - Task $i: bắt đầu tại t = ", round(value(z[i]), digits=2))
        end
        
        println("\n Các luồng di chuyển của robot (f_ij > 0):")
        for i in N_plus, j in N_minus
            if i != j
                val = value(f[i, j])
                if val > 0.001
                    println("   - Từ Task $i -> Task $j: ", round(Int, val), " robot")
                end
            end
        end
    else
        println("Không tìm thấy lời giải tối ưu.")
    end
end


input_dir = "./instances/pdp_100"
output_dir = "./logs/pdp_100"

# Tạo thư mục output nếu chưa tồn tại
mkpath(output_dir)

# Lấy danh sách tất cả file JSON trong thư mục input
json_files = filter(
    f -> endswith(lowercase(f), ".json"),
    readdir(input_dir)
)

# Chạy lần lượt từng file
for file in sort(json_files)
    json_path = joinpath(input_dir, file)

    # Ví dụ: lc103.json -> lc103_result.txt
    filename = splitext(file)[1]
    output_path = joinpath(output_dir, "$(filename)_result.txt")

    println("Đang xử lý: $file")

    open(output_path, "w") do io
        redirect_stdout(io) do
            solve_flowbased_milp(json_path)
        end
    end

    println("Đã lưu kết quả: $output_path")
end
