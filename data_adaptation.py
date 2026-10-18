import os
import glob
import math
import random
import json
from pathlib import Path

def convert_single_file_to_pdp_cr(input_file_path, fleet_size_m=10, seed=42):
    """
    Chuyển đổi 1 file benchmark Li & Lim gốc sang định dạng đầy đủ thông tin của PDP-CR
    tương ứng với Bảng 1 và Bảng 2 trong bài báo.
    """
    random.seed(seed)
    
    nodes = {}
    pickup_nodes = []
    
    with open(input_file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    if not lines:
        return None
        
    # Đọc thông tin các nút (bỏ qua dòng 0 chứa metadata gốc của Li & Lim)
    for line in lines[1:]:
        parts = line.strip().split()
        if not parts:
            continue
            
        node_id = int(parts[0])
        x, y = float(parts[1]), float(parts[2])
        demand = float(parts[3])
        p_id, d_id = int(parts[-2]), int(parts[-1])
        
        nodes[node_id] = {'x': x, 'y': y, 'demand': demand, 'p_id': p_id, 'd_id': d_id}
        if d_id > 0:  # Nút Pickup
            pickup_nodes.append(node_id)
            
    if 0 not in nodes:
        return None
        
    # Tọa độ Depot (Tọa độ 0 trong file gốc)
    depot_x, depot_y = nodes[0]['x'], nodes[0]['y']
    depot_coord = [depot_x, depot_y]
    
    # Tạo danh sách các nhiệm vụ (i = 1 ... n)
    tasks = []
    for i, p_node_id in enumerate(pickup_nodes, start=1):
        d_node_id = nodes[p_node_id]['d_id']
        
        p_coord = [nodes[p_node_id]['x'], nodes[p_node_id]['y']]
        d_coord = [nodes[d_node_id]['x'], nodes[d_node_id]['y']]
        
        # c_i = Khoảng cách Euclide từ p_i đến d_i
        c_i = math.dist(p_coord, d_coord)
        
        # r_i = ngẫu nhiên trong khoảng
        r_i = random.randint(1, 5)
        
        tasks.append({
            "i": i,
            "p_i": p_coord,
            "d_i": d_coord,
            "r_i": r_i,
            "c_i": round(c_i, 2)
        })
        
    n = len(tasks)
    
    # --- Tạo bảng danh sách tất cả các điểm All Nodes (0, 1..n, n+1) chuẩn theo Bảng 1 ---
    depot_start = {
        "i": 0,
        "p_i": depot_coord,
        "d_i": depot_coord,
        "r_i": 0,
        "c_i": 0.0
    }
    
    depot_end = {
        "i": n + 1,
        "p_i": depot_coord,
        "d_i": depot_coord,
        "r_i": 0,
        "c_i": 0.0
    }
    
    all_nodes_info = [depot_start] + tasks + [depot_end]
    
    # --- Tính ma trận t_ij kích thước (n+2) x (n+2) chuẩn theo Bảng 2 ---
    # t_ij: Thời gian di chuyển từ d_i đến p_j
    t_matrix = []
    for i_node in all_nodes_info:
        row = []
        i_idx = i_node["i"]
        d_i_coord = i_node["d_i"]
        for j_node in all_nodes_info:
            j_idx = j_node["i"]
            p_j_coord = j_node["p_i"]

            if i_idx == j_idx or (i_idx == 0 and j_idx == n + 1):
                dist = 0.0
            else:
                dist = math.dist(d_i_coord, p_j_coord)
            row.append(round(dist, 2))
        t_matrix.append(row)
        
    # Đóng gói kết quả với đầy đủ cấu trúc bài báo
    data_out = {
        "metadata": {
            "m_robots": fleet_size_m,
            "n_tasks": n
        },
        "tasks": all_nodes_info,
        "traveling_time_matrix_t_ij": t_matrix
    }
    
    return data_out

def batch_adapt_dataset(input_dir, output_dir, fleet_size_m=10, seed=42, file_extension="*.txt"):
    """
    Chuyển đổi toàn bộ thư mục input_dir và lưu file JSON đầy đủ vào output_dir.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    file_list = [f for f in glob.glob(os.path.join(input_dir, file_extension)) if os.path.isfile(f)]
    
    success_count = 0
    for file_path in file_list:
        file_name = Path(file_path).stem
        output_file_name = f"{file_name}.json"
        output_file_path = os.path.join(output_dir, output_file_name)
        
        converted_data = convert_single_file_to_pdp_cr(
            file_path, 
            fleet_size_m=fleet_size_m, 
            seed=seed
        )
            
        if converted_data:
            with open(output_file_path, 'w', encoding='utf-8') as f_out:
                json.dump(converted_data, f_out, indent=4, ensure_ascii=False)
            
            print(f" [Thành công] {Path(file_path).name} -> {output_file_name}")
            success_count += 1
            
    print(f"\n HOÀN THÀNH: Đã lưu {success_count}/{len(file_list)} file chuẩn hóa vào '{output_dir}'.")

if __name__ == "__main__":
    INPUT_FOLDER = ""      # Đường dẫn thư mục chứa file gốc
    OUTPUT_FOLDER = "" # Thư mục lưu kết quả JSON
    
    batch_adapt_dataset(
        input_dir=INPUT_FOLDER,
        output_dir=OUTPUT_FOLDER,
        fleet_size_m=10,   # Số lượng robot m = 10 theo bài báo
        seed=42            # Khởi tạo seed ngẫu nhiên
    )