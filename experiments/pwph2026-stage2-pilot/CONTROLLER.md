# Quy tắc heuristic v1, viết trước đánh giá adaptive

Đầu tiên xử lý 0→1.2, giữ nguyên chọn support/PV của PwPH. Sau mỗi cửa sổ
không cuối, quotient (thu gọn vùng thành nhãn đại diện) như cũ, quan sát:

`r = 1 - active_after / active_before` (tỷ lệ không đơn vị).
`width = clamp(0.5 * r / target, 0.2, 0.5)` (đơn vị epsilon), target=0.35.

Chọn điểm xa nhất trên lưới {1.2,1.4,1.5,1.6,1.7,1.8,2.0} không vượt
`current+width`; nếu không có thì lấy điểm kế tiếp. Bắt buộc tiến lên tới 2.
CLI cho đổi `--adaptive-target` trong (0,1]; pilot khóa 0.35.

Lý do: quotient giảm nhiều nhãn cho phép đi xa hơn; khi ít giảm, xử lý một
khoảng ngắn để có dịp quotient tiếp. 0.5 là bước của lịch n100, 0.2 là bước
ngắn trên lưới, target 35% là mốc thử đơn giản gần mức giảm n100 đã biết.
Không fit G, không dùng PH hay thời gian tương lai, không chọn ngưỡng theo held-out.
Các seed phát triển có thể dùng để đánh giá thất bại; không buộc sửa heuristic
cho tới khi vượt đối chứng. Giữ cả kết quả âm.

Active counts hoạt động cả khi tắt verbose. Log cuối tiến trình ghi quyết định,
counts, top cofacets, matching/rebuild/window time và controller microseconds;
record cửa sổ tiếp chính là giá trị thực sau quyết định trước. Không dự báo
cofacets/thời gian. Không đọc RAM để điều khiển v1; giám sát RSS cây tiến trình
bên ngoài dùng cùng giới hạn cho mọi phương pháp. Đây là giới hạn của heuristic.
Logging hoãn tới sau loop; runtime toàn tiến trình gồm telemetry, logging,
checkpoint, quotient, rebuild và controller. Verbose chạy riêng để lấy MSF H0.

`WindowState` giữ active labels, PV representatives và witness semantics qua
cửa sổ. Chỉ `DimensionFrame` và cache geometry kết thúc ở ranh giới, đúng luồng
PwPH có sẵn. Không chạy lại point cloud từ đầu ở mỗi quyết định.

Khảo sát lịch cố định ban đầu trên binary n100 phát hiện H2 sai ở một số lịch
1.4/1.5. Không thay kernel hay sửa persistence bookkeeping trong nhiệm vụ này;
adaptive cũng phải qua kiểm tra PH riêng. Đây là pilot sau checkpoint đầu 1.2,
chưa hiện thực hàm lợi ích G, chưa làm giai đoạn 3.
