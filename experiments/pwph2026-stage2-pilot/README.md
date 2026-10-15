# Pilot lịch PwPH giai đoạn 2

Đọc `BAO-CAO-TRIET.md` để lấy kết luận. `CURRENT-RUN.txt` chỉ thư mục evidence.
Không dùng số đo cũ của n100 để suy ra tăng tốc adaptive. Tham số heuristic trong
`CONTROLLER.md`; `LOCK.json` trong run lưu thời điểm khóa đối chứng và hash.

Các chương trình chạy từ root repo, bằng Python đã có Ripser 0.6.15:

```bash
/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python experiments/pwph2026-stage2-pilot/resume.py
```

Lệnh tự lưu stdout/stderr/context/exit-status trong thư mục `resume-TIMESTAMP`.
Chỉ tiếp tục trong ngân sách 6 giờ ban đầu, dành 10 phút đóng gói; không tự chạy
định kỳ. Không retry ca lỗi, không ghi đè evidence, không chạy lại ca hoàn tất
cùng identity. Nếu bị ngắt giữa ca (context chưa có exit_status), dừng để kiểm
tra ca đó, giữ lại và ghi lý do trước khi retry với tên riêng. Trước resume cần
đảm bảo không còn runner khác hoạt động. Không chạy hai benchmark đồng thời.

Identity mỗi process gồm argv, binary/data/script hash; `cases/*/context.json`
giữ môi trường, UTC timestamps, monotonic time, GNU Time peak RSS và exit code.
`index.jsonl` cập nhật qua replace nguyên tử, `PROGRESS.md` cập nhật từng ca.
`ph-check.json` giữ trạng thái PH từng lần đo, round0 là warmup; round1..5 là
năm lần đo. Thứ tự phương pháp xoay vòng và đảo chiều, tuần tự 4 luồng/process.
TIMEOUT/MEMORY_LIMIT/ERROR tách riêng với exit code gốc; PH_MISMATCH có thể có
exit code iCC=0. H0_NOT_VERIFIED không đủ điều kiện vào bảng PASS H0–H2.

Parser tổng quát hóa quy tắc trong benchmark_icc.py: Counter Decimal giữ
multiplicities, chỉ lấy finite death thuộc cửa sổ và survivors ở cửa sổ cuối.
Không nối/khử trùng survivor theo birth. Verbose được kiểm tra riêng bằng parser
compare_icc_ph.py trên actual schedule và so compact; raw log luôn giữ nguyên.
MSF H0 dùng lại validate_n100.msf; chỉ PASS khi tất cả hợp nhất đã xong trước
quotient đầu và khớp baseline/Ripser. Hai known differences vẫn ở regression.

Tạo lại bảng (không chạy PH) bằng `analyze.py --run RUN --binary RUN/build/iCC`.
`--lock` chỉ dùng một lần trước held-out; không được ghi đè lock để chọn lại.
`shrink.py` giữ hai điều tra lịch cố định, tối đa 20 phút/ca; không tuyên bố
tối tiểu toàn cục. Không sửa kernel hoặc quy ước cutoff/zero-distance.

CLI có auto-generated schedule (`--eps-interval-count`, `--eps-interval-scale`),
nhưng nó kết thúc tại nextafter(max pair distance), checkpoint đầu theo rank
khoảng cách. Không thể ép đồng thời checkpoint đầu 1.2 và cutoff 2.0 qua các
flag đó. Vì thế không đưa vào bảng đối chứng cùng điều kiện của pilot này.

Build mới: Release, BUILD_TESTING=ON, CLI11 lấy từ build/_deps/cli11-src;
commands chính xác trong configure/context.json. Không dùng build/iCC cũ.
Archive không chứa binary/build/venv/cache; muốn tái lập cần build source đã gói,
đối chiếu compiler/flags trong môi trường run rồi tạo batch mới khi hash đổi.
