# G chọn toàn bộ lịch trước khi chạy

Giữ nguyên engine Morse/iCC commit `259d4bf`, F2, H0–H1, max-simplex 2,
S² đơn vị, alpha 0.85, cutoff 2. Không phát triển CONTINUE hoặc thay luật PV.
`selector.py` chỉ đọc tọa độ, dự đoán thời gian của 10 lịch, rồi `execv` engine.
Không dựng checkpoint ngoài lịch. Bộ chọn chỉ dùng thư viện chuẩn Python.

Run hiện hành xem CURRENT-RUN.txt. INDEX.md trong run dẫn tới báo cáo, bảng,
mô hình khóa, dữ liệu và log gốc. Thời gian G bao gồm khởi động Python,
đọc dữ liệu, O(n² × 11) đặc trưng khoảng cách, dự đoán, ghi lựa chọn,
thay tiến trình bằng engine, PH, xuất kết quả và giải phóng.

`offline.py --run RUN prepare|build|survey|fit|final` là các bước nghiên cứu;
`summarize.py --run RUN` tổng hợp mà không chạy PH. Dùng Python có NumPy và
Ripser 0.6.15 cho offline.py; selector.py không cần hai thư viện này.
Mã/mô hình khóa và kết quả chính được giữ trong models/, tables/, RESULTS.md.
Các bước thực nghiệm giữ budget.json đã ghi từ đầu, kể cả lúc gặp lỗi.

Trong gói, `replay_handoff.py --bundle THU_MUC_GOI --out THU_MUC_MOI` dựng
engine từ source kèm CLI11 rồi chạy lại A–E trên đúng sáu input cuối. Không
huấn luyện lại hoặc đổi lịch. REPRODUCE.md có lệnh và giới hạn phép đo.
Kiểm tra gói trong lượt này dùng `--budget-from` trỏ ngân sách ban đầu;
một lượt tái hiện độc lập về sau có ngân sách riêng, giữ nguyên hồ sơ gốc.
