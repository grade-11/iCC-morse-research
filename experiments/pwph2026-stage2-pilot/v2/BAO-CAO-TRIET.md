# Báo cáo cho Triết — PwPH v2, 16/10/2026

**Chưa sửa được mất thông tin PH nền. Đã xây G điều khiển quotient/defer thật. Final G: max-dim 2 (H0–H1) PASS 8/8; max-dim 3 (H0–H2) PASS 8/8.** Chưa có lợi ích tốc độ so fixed tốt: fixed/G=0.953 ở H0–H1 và 0.962 ở H0–H2 (<1 bất lợi cho G). Không suy ra bảo toàn PH tổng quát.

Có 40 dataset mới ({'dev': 24, 'final': 8, 'validation': 8}), n32/64/100/200, S2 đều, S2 tập trung cực và cụm đa thang đo. Seed07–11 chỉ là regression; seed final30–37 chỉ sinh sau lock. Không chỉnh thiết kế sau final. Đến tổng hợp đã dùng 112.2 phút/240 phút (monotonic), gồm truy nguyên, build, calibration và đo; deadline cũ giữ nguyên khi resume.

## Độ đúng đã kiểm, lỗi còn lại

Hai phản ví dụ cũ tái lập ở 1/4 luồng. Oracle độc lập khớp từng simplex/weight witness và reduction F2 của C++; nhưng phức sau quotient không còn đủ thông tin VR: ca 11 điểm còn ba nhãn trước H2 mới; ca 10 điểm còn hai nhãn trước hai H1 và một H2 mới. PV vẫn qua các bộ lọc đang cài đặt. Không có căn cứ vá parser/matching hoặc đổi alpha để ép PASS. Sáu lỗi n100 đầy đủ vẫn giữ nguyên; chưa chứng minh cùng nguyên nhân với ca nhỏ. Max-dim 2 là cấu hình kiểm riêng, không che lỗi H1 đã biết ở max-dim 3.

Bổ sung riêng bộ kiểm H0 từ MSF thật từng cửa sổ, kiểm conservation và toàn multiset với baseline/Ripser. Giả thiết: khoảng cách dương, alpha∈(0,1), luật disjoint/frozen PV hiện tại. Sáu unit tests và audit 576 cấu hình đạt; 18 H0 chưa xác minh cũ được bổ sung, evidence gốc giữ nguyên. Đây là kiểm chứng H0 có điều kiện, không sửa H1/H2. CTest, n100 chuẩn, regression cũ và nhánh toàn DEFER được kiểm. Source+CLI11 đóng gói build lại offline đạt 3/3 CTest, binary SHA256 trùng binary benchmark. Endpoints/multiplicities giữa iCC so đầy đủ; Ripser theo float32 đã chốt, không nới tolerance. Fixture cutoff/điểm trùng vẫn KNOWN_DIFFERENCE.

## Đo final: thời gian toàn tiến trình và RSS

Mỗi nhóm đủ PH có 1 warmup +5 vòng xoay/đảo thứ tự, tuần tự 4 luồng, cùng binary, không verbose. Bảng max-dim 2: giây trung vị; tỷ số paired fixed/G >1 lợi cho G; RSS MiB. — là không đủ điều kiện, không phải 0.

| Final (hình học/n/seed) | Fixed / G (s) | Fixed/G [min,max] | RSS fixed / G | G PH |
|---|---:|---:|---:|---|
| uniform_s2/32/30 | 0.011 / 0.014 | 0.86 [0.67,0.90] | 5.137 / 5.438 | PASS |
| uniform_s2/64/31 | 0.018 / 0.021 | 0.85 [0.60,1.11] | 5.828 / 6.086 | PASS |
| uniform_s2/100/32 | 0.033 / 0.034 | 1.02 [0.80,1.08] | 7.059 / 7.266 | PASS |
| uniform_s2/200/33 | 0.123 / 0.127 | 0.99 [0.93,1.08] | 20.523 / 20.734 | PASS |
| polar_s2/64/34 | 0.021 / 0.020 | 1.04 [0.97,1.06] | 5.766 / 5.957 | PASS |
| polar_s2/100/35 | 0.046 / 0.046 | 1.00 [0.90,1.26] | 9.520 / 9.750 | PASS |
| multiscale/64/36 | 0.020 / 0.021 | 0.95 [0.92,1.03] | 5.805 / 5.949 | PASS |
| multiscale/100/37 | 0.036 / 0.039 | 0.94 [0.84,1.00] | 7.320 / 7.496 | PASS |

Đã khóa theo chiều và nhóm n≤64/n≥100: `{'2-small': {'fixed': 'sparse', 'G': 'gs5'}, '2-large': {'fixed': 'sparse', 'G': 'gs5'}, '3-small': {'fixed': 'sparse', 'G': 'gs5'}, '3-large': {'fixed': 'sparse', 'G': 'gs5'}}`. Chọn bằng coverage PH, rồi H0, rồi thời gian development/validation. Bảng đầy đủ gồm fixed gốc, fixed tốt, v1, G, H2 và min–max ở performance.csv/paired-ratios.csv. Các tỷ số tổng hợp là geometric mean của median paired ratio từng dataset: v1/G=1.319 (D2), 0.797 (D3); fixed gốc/G=1.106/0.829. G chọn cùng hành động với fixed tốt ở 7/8 ca mỗi chiều. Final chỉ có một checkpoint nên không đo lợi ích EMA cho quyết định sau. Cùng lịch/hành động mà lệch nhỏ chưa chứng minh giá trị scheduling; tỷ số có điều kiện trên ca đủ đúng, không bao gồm ca sai.

## G và chất lượng dự báo

`G = Ccontinue − (Cselect+quotient + Crebuild + Cfuture_afterQ)` từ cùng trạng thái tới cùng h, đơn vị giây. h là checkpoint kế tiếp; G chọn có quotient hay trì hoãn, không chọn PV theo chi phí. Ridge đơn giản theo nhãn/cofacets/PV/support/độ rộng, cập nhật scale từ quá khứ. Quotient khi G>max(0.5ms; margin×Ccontinue); RAM KiB là ràng buộc riêng với 30% headroom. Ngoài miền dự báo fallback DEFER; không gọi là an toàn PH tổng quát. Toàn telemetry/controller/checkpoint/logging tính vào runtime.

Ví dụ `uniform_s2-n32-s20261030` D2, ε=1.2→h=2, A=32, C=368: continue dự báo 2.86ms; Q/rebuild/future=0.01/2.48/1.89ms; G=-1.52ms, ngưỡng 0.50ms → **DEFER** (G-below-margin). Work cửa sổ sau quan sát 1.41ms; final chỉ đo nhánh đã chọn.

Ridge .01 ban đầu: D2 32 cặp, MAE 4.8ms, đúng dấu 78.1% (luôn dương 84.4%); D3 32 cặp, MAE 169.9ms, đúng dấu 90.6% (luôn dương 100.0%).
Bổ sung chân trời cuối: D2 40 cặp, MAE 9.7ms, đúng dấu 77.5% (luôn dương 87.5%); D3 40 cặp, MAE 657.8ms, đúng dấu 75.0% (luôn dương 100.0%).
Bỏ nhóm PV/support làm MAE validation xấu hơn nên giữ feature đầy đủ; đúng dấu riêng có thể tốt hơn, không giấu kết quả đó. Nhãn G_actual là hiệu chi phí các phase không chồng lấn, chưa phải chi phí phản thực whole-process chính xác: phí hủy support-container Q và phần controller/log chưa tách. Biên 0,5ms chưa là cận nhiễu được chứng minh. Chi phí replay/fitting offline ghi riêng; online không chạy thử hai nhánh.

## Kiểm kê và bàn giao

Thu nhỏ lỗi mới: uniform_s2-n64-s20261013 D3: 64→11 điểm; polar_s2-n100-s20261018 D3: 100→14 điểm. Input/trace giữ nguyên, không tuyên bố tối tiểu.
RAM study: COMPLETED, đổi hành động=True, trạng thái tương đương=True; đây là phản ứng dự báo, chưa chứng minh tiết kiệm RAM thực.
Hai cặp validation khó được lặp để đo nhiễu: D2 G thực min–max 72.0–73.5ms; D3 G thực min–max 3842.8–3947.7ms.
704 cấu hình logic (704 thực thi riêng), trạng thái {'PASS': 668, 'PH_MISMATCH': 36}; 64 warmup +320 đo final, trạng thái {'PASS': 384}. Alias không tính thêm process/dataset. TIMEOUT/MEMORY_LIMIT/ERROR giữ riêng trong inventory/raw logs. Lỗi tooling snapshot/replay/orchestration và lần auto-review commit timeout được lưu; retry không xóa evidence.

Chưa có điều kiện đủ bảo toàn PH mọi quotient; estimator còn sai số, mẫu final nhỏ và chỉ một máy. Chưa làm giai đoạn 3 hoặc hoàn tất giai đoạn 2. Đọc PH-ROOT-CAUSE.md, H0-LEDGER.md, CONTROLLER-v2.md; inventory/correctness/performance/G-predictions/AUDIT trong run. README có resume giữ deadline. Một archive kèm manifest/SHA256 chứa source, patch, tests, input, model và raw logs, không build/venv/cache. Commit cục bộ, không push.
