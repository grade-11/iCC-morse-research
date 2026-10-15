# Báo cáo pilot giai đoạn 2 cho Triết — 16/10/2026 (+07)

**Đã code và chạy heuristic online sau checkpoint đầu 1.2; chưa hoàn tất toàn bộ giai đoạn 2, chưa làm giai đoạn 3.** Đối chứng chính là PwPH lịch cố định đã khóa `[1.2,1.7,2.0]`. Cần đọc kết quả theo từng seed; không lấy mức nhanh hơn iCC không thu gọn làm đóng góp mới của adaptive.

Đã khảo sát năm lịch cố định trên ba seed phát triển n100; khóa lịch tốt đủ điều kiện PH bằng trung bình nhân của trung vị thời gian từng seed. Heuristic target=0.35 được giữ nguyên. Sau đó đánh giá hai seed độc lập, rồi n200 cùng nhóm seed, cùng binary. Không mở rộng n400/dữ liệu thực tế; không push/merge/công bố.

Baseline = cách đối chứng; checkpoint = mốc xử lý; quotient = thu gọn vùng thành nhãn đại diện; rebuild = dựng lại phức; peak RSS = bộ nhớ cư trú lớn nhất; regression = kiểm tra sửa code có làm hỏng kết quả đã đúng; ablation = thử riêng từng thành phần.

## Độ đúng

Binary mới có **52 cấu hình PASS H0–H2**, **6 PH_MISMATCH** (các trạng thái khác: không). Mỗi cấu hình được kiểm tra riêng; mọi vòng đo còn so H1/H2 với baseline và Ripser. H0 dùng MSF thật của iCC trong diagnostic riêng, chỉ nhận PASS khi graph gốc đã liên thông trước quotient đầu, khớp baseline chính xác và Ripser sau float32. Không cộng local H0 của các quotient.

So giữa các iCC bằng multiset Decimal của endpoints đầy đủ đã in, giữ multiplicities và survivor cuối. Với Ripser 0.6.15 chiếu endpoints iCC sang float32, không tolerance. Chỉ bỏ survivor tạm ở các cửa sổ chưa cuối, theo actual schedule; parser verbose và compact khớp. Không tuyên bố đúng float64 với Ripser hoặc chứng minh PH tổng quát.

Regression n100 và CTest 3/3 đạt; kiểm tra harness 3/3 đạt, gồm timeout, lỗi exit7, memory limit giả lập, resume và parser nhiều cửa sổ. Hai fixture cạnh đúng cutoff/điểm trùng vẫn là **KNOWN_DIFFERENCE**: iCC strict `<`, Ripser gồm equality; piecewise bỏ cạnh 0 và Ripser không xuất zero-length H0. Không đổi quy ước.

## Thời gian và RAM so đối chứng chính

Mỗi cấu hình đủ PH: một warmup + năm vòng đo, thứ tự xoay/đảo, từng process tuần tự, 4 luồng. Whole-run không verbose dùng monotonic clock và GNU Time; gồm mọi chi phí controller/telemetry/logging/quotient/rebuild. Bảng dưới là seed độc lập; mỗi ô là **trung vị [min–max]**. Tỷ số theo cặp = fixed/adaptive; >1 có lợi cho adaptive.

| n / seed | Fixed (s) | Adaptive (s) | Tỷ số thời gian theo cặp | Fixed RSS (MiB) | Adaptive RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| 100 / 20261010 | 0.1623 [0.0858–0.1665] | 0.1616 [0.0770–0.1686] | 1.007 [0.929–1.114] | 10.50 [10.46–10.59] | 10.39 [10.25–10.62] |
| 100 / 20261011 | 0.2007 [0.1843–0.2084] | 0.2113 [0.2008–0.2681] | 0.957 [0.741–1.000] | 13.84 [13.48–14.12] | 12.11 [11.99–12.27] |
| 200 / 20261010 | 1.7930 [1.7376–1.8359] | 1.7745 [1.6721–1.8042] | 1.018 [0.979–1.061] | 83.11 [82.95–83.40] | 83.16 [82.72–83.36] |
| 200 / 20261011 | 1.9754 [1.8932–2.1060] | 2.0140 [1.9155–2.0451] | 1.007 [0.938–1.031] | 83.83 [83.61–83.87] | 83.40 [83.12–83.45] |

Trong 10 ca adaptive, 9 ca chọn cùng lịch cố định đã khóa. Chênh lệch thời gian ở các ca trùng lịch không thể quy cho thay đổi mốc. Riêng n100 seed độc lập 20261011, adaptive chậm hơn khoảng 5,3% theo tỷ số hai trung vị nhưng giảm peak RSS khoảng 12,5%; đây là đánh đổi RAM/thời gian, chưa là cải thiện cả hai. Pilot này chưa đủ bằng chứng để khẳng định heuristic thắng lịch PwPH cố định cấu hình tốt.

Đủ số đo từng seed phát triển/độc lập, min–max và tỷ số RAM theo cặp trong `summary.csv`, `paired-ratios.csv`; từng vòng và warmup trong `measurements.csv`. Không trộn với benchmark n100 cũ. Lịch đối chứng chuyển nguyên từ phần phát triển n100 sang n200, chưa chọn lại lịch tốt nhất riêng cho n200. Lịch tốt nhất theo từng ca nhìn lại chỉ là tham khảo hậu nghiệm; không thay lịch khóa bằng lịch thắng trên held-out.

## Controller đã chọn thế nào

`r = 1 − active_after/active_before`; `width = clamp(0.5*r/0.35, 0.2, 0.5)` epsilon. Lấy mốc xa nhất trong lưới tăng không vượt current+width, hoặc mốc kế tiếp để bảo đảm tiến lên, kết thúc 2.0. Nhãn giảm nhiều thì cho khoảng sau rộng hơn; ít giảm thì xử lý gần hơn để có dịp quotient tiếp. Chỉ dùng trạng thái hiện tại, không dùng PH/thời gian tương lai; không restart từ đầu.

Ví dụ seed độc lập n100/20261011 (log nonverbose, đơn vị ms và µs):

```text
[adaptive] eps_lo=0 eps_hi=1.2 next=1.6 reason=reduction-scaled-width active_before=100 active_after=68 reduction=0.32 target=0.35 top_cofacets=54077 matching_ms=17.1763 rebuild_ms=12.1042 elapsed_ms=29.7628 controller_us=0.304
[adaptive] eps_lo=1.2 eps_hi=1.6 next=1.8 reason=reduction-scaled-width active_before=68 active_after=52 reduction=0.235294 target=0.35 top_cofacets=103003 matching_ms=29.2151 rebuild_ms=29.1238 elapsed_ms=58.6972 controller_us=0.261
[adaptive] eps_lo=1.6 eps_hi=1.8 next=2 reason=reduction-scaled-width active_before=52 active_after=36 reduction=0.307692 target=0.35 top_cofacets=116441 matching_ms=18.9761 rebuild_ms=41.1189 elapsed_ms=60.5106 controller_us=0.249
[adaptive] eps_lo=1.8 eps_hi=2 next=-1 reason=final active_before=36 active_after=36 reduction=0 target=0.35 top_cofacets=58905 matching_ms=12.021 rebuild_ms=29.527 elapsed_ms=41.5587 controller_us=0
```

Đây là **heuristic dựa trên mức giảm nhãn**, chưa hiện thực hàm lợi ích G trong AC_PwPH. Cofacets và timing được ghi để giải thích/đánh giá, chưa tham gia quy tắc; RAM là ràng buộc giám sát bên ngoài. Không bật FFI/verbose để lấy tín hiệu điều khiển. Giữ nguyên kernel Morse, witness weights, điều kiện/thứ tự PV và quotient của từng cửa sổ; chỉ thay mốc tiếp theo. Log quyết định hoãn tới sau vòng xử lý, record cửa sổ sau là kết quả quan sát thực sau quyết định trước; không có dự báo.

## Kết quả âm và lỗi

Ba lịch có checkpoint 1.4/1.5 sai H2 ở n100 seed20261007/20261008 trên cả binary gốc và mới. Ví dụ `[1.2,1.5,2.0]`, seed07: birth H2 0.97696823381778153 thành 0.93668824329311717, death giữ 1.6602624057539093. Sai khác vẫn tồn tại sau float32. Loại khỏi bảng hiệu năng hợp lệ, không ép PASS.

- Ca 20261007/f15: tái lập, giảm 100→11 điểm, 114 lần thử, 15.8 s; giữ `shrunk/20261007-f15/retained.csv` và toàn bộ input/log. Đây là ca đã thu nhỏ, không tuyên bố tối tiểu toàn cục.
- Ca 20261008/f4: tái lập, giảm 100→10 điểm, 179 lần thử, 17.8 s; giữ `shrunk/20261008-f4/retained.csv` và toàn bộ input/log. Đây là ca đã thu nhỏ, không tuyên bố tối tiểu toàn cục.

Ripser mới trên cả hai ca thu nhỏ khớp baseline H1/H2; PwPH mất H2 (11 điểm), hoặc mất cả H1/H2 (10 điểm), không phải sai số float32. H0 các ca thu nhỏ không dùng để kết luận hiệu năng.

Audit chỉ mục lần đầu exit1 vì kiểm hash cấu hình là duy nhất, quên các vòng lặp chủ ý. Đã sửa khóa thành tên ca/vòng + hash cấu hình, chạy lại audit đọc dữ liệu; giữ log lỗi và source v1. Không chạy lại phép đo PH/timing vì lỗi audit này.

Không sửa lỗi PH nền trong lượt thử lịch này. Bộ đối chứng n100 cũ không bị reset/ghi đè. Lần tạo nhánh đầu exit128 do .git chỉ đọc; thực hiện lại qua quyền sandbox đã được phép và thành công. Probe môi trường phát hiện không có psutil: harness dùng /proc, không cài package. Matplotlib từng cảnh báo cache không ghi được; chuyển MPLCONFIGDIR sang /tmp cho các đợt sau. TIMEOUT/MEMORY_LIMIT giả lập chỉ thuộc unit test, không gọi là lỗi thực nghiệm.

## Giới hạn và điều cần trao đổi với thầy

Pilot chỉ thích nghi sau checkpoint đầu **1.2**, cutoff **2.0**, alpha **0.85**, simplex max-dim **3**, F2; S² đơn vị trong R³ tự sinh, không phải dữ liệu gốc tác giả. Generator `dSphere_gen.dSphere_uniformRandom(dimension=3, r1=r2=1, trials=0, selection=min)`, NumPy seed đã ghi, hash CSV/source và phiên bản package trong datasets/environment/context. Trials=0 khớp chính xác tọa độ n100 đã review; không thay phân phối hay sửa bộ sinh.

Cần trao đổi về điều kiện bảo toàn H2 khi thay lịch, H0 nối qua quotient khi chưa liên thông, và thiết kế estimator G/ablation trước khi tuyên bố đóng góp. Chưa đánh giá nhiều phần cứng, dữ liệu khó, ngân sách RAM khác, hoặc ý nghĩa thống kê của chênh lệch nhỏ. Số test nhiều không là chứng minh tổng quát.

Giới hạn RAM cây tiến trình 6,37 GiB (min(8 GiB,70% MemAvailable lúc bắt đầu)), giám sát /proc khoảng 50 ms, **gần đúng**, không hard cgroup limit; mỗi process thử tối đa 300 s. Không thay swap/WSL hoặc kill process ngoài nhiệm vụ. Peak RSS đo cho toàn process C++ bằng GNU Time, khác tổng RSS lấy mẫu để giám sát. Không cộng các timer lồng nhau.

CLI auto-generated lịch đã được đọc: dùng `--eps-interval-count`/`--eps-interval-scale`, nhưng không ép được cả mốc đầu 1.2 và cutoff 2.0; không đưa vào bảng cùng điều kiện. Ripser chỉ kiểm PH, không benchmark tốc độ wrapper.

## Bàn giao và tiếp tục

Run: `experiments/pwph2026-stage2-pilot/runs/stage2-20261015T190630Z`. `PROGRESS.md`, `index.jsonl`, `LOCK.json`, `correctness-overview.json` và raw `cases/` là điểm tiếp tục; không cần đọc từng log để lấy kết luận. Source/patch, dữ liệu, scripts, log, checksum được gói trong **một** archive `stage2-20261015T190630Z-handoff.tar.gz` cạnh báo cáo này; kiểm SHA256 cạnh archive và manifest bên trong. Không chứa build/venv/cache; build context và hash binary vẫn được lưu. Commit cục bộ được ghi trong `git-final.json`.

Lệnh resume có tự lưu log, kiểm hash và bỏ qua các ca hoàn tất:

```bash
/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python experiments/pwph2026-stage2-pilot/resume.py
```

Resume chỉ trong ngân sách 6 giờ ban đầu và phải không còn runner khác. Khi run đã đủ ca, không tự benchmark lại. Sau khi hết ngân sách, nghiên cứu tiếp cần lượt mới, giữ nguyên evidence này.
