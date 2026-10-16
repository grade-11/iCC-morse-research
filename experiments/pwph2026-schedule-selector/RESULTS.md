# Báo cáo G chọn lịch trước khi chạy — 16/10/2026

Đã triển khai, biên dịch, khảo sát hữu hạn, khóa mô hình trước khi sinh dữ liệu cuối và đo riêng toàn bộ chi phí G. Tập cuối có 6 bộ độc lập, H0–H1 (đồng điều bền vững chiều 0 và 1), trường F2 gồm hai phần tử; G đạt kiểm PH 6/6 bộ. G chọn [1.0,2.0] cho cả sáu bộ, nhưng chậm hơn lịch cố định tốt trung bình nhân 3.42 lần khi tính đủ chi phí. Đây là kết quả bất lợi, chưa có bảo đảm PH cho mọi dữ liệu.

Ngân sách cố định bắt đầu **2026-10-16T06:47:00+00:00**, dừng mở thí nghiệm **2026-10-16T08:27:00+00:00**, hạn **2026-10-16T08:47:00+00:00**; không đặt lại khi lỗi. Thời điểm lập bảng 2026-10-16T07:10:41.948397+00:00. Repo/HEAD đầu lượt trong initial-context.json; commit bàn giao trong git-final.json. Không push, không reset/clean.

1. **G chọn lịch dựa trên gì?** Quotient là phép lấy phức thương bằng các đỉnh giả PV. Bộ chọn đọc tọa độ trước PH, lấy n và mật độ cạnh ở các ngưỡng k/10 (k nguyên từ 8 đến 18), tính vector 10 đặc trưng: hệ số chặn, log(n/100), mật độ cạnh tại 0.8, số mốc quotient, có quotient hay không, mốc quotient đầu/2, tổng mật độ cạnh ở các mốc quotient, hai tương tác với log(n) và mật độ 0.8, cùng proxy biên Σdε(1−dε)(2−ε). Proxy là đại lượng thay thế từ khoảng cách gốc; không phải số witness, PV hay cofacet tương lai đã biết. Cofacet là đơn hình kề trên một mặt; witness là tổ hợp đỉnh gốc đại diện dùng xác định đơn hình và trọng số.

Hồi quy ridge (bình phương tối thiểu có phạt hệ số) dùng 10 hệ số cho log thời gian toàn tiến trình của engine (bộ tính toán Morse/iCC). Lấy hàm mũ để dự đoán giây cho từng lịch, G(x,L)=T̂(noQ)−T̂(L). Các dự đoán này cùng phạm vi đọc dữ liệu, Morse/iCC tới 2, xuất PH và giải phóng của engine; phần Python chọn lịch là chi phí chung cộng thêm của phương án D, đo trọn bên ngoài và báo riêng, chưa được hồi quy thành dự đoán lợi ích ròng. Đây là giới hạn của G hiện tại: G dự đoán dương chưa đủ kết luận tổng thời gian D thắng A. Chọn lịch có thời gian dự đoán thấp nhất trong tập đủ kiểm tra phát triển/chọn cấu hình. Không gọi đây là controller (bộ điều khiển) thích nghi sau mỗi quotient. Không dựng các checkpoint (mốc kiểm tra) ngoài lịch để quyết định bỏ qua.

Có 12 bộ development (phát triển), 6 bộ selection (chọn cấu hình), 6 bộ final (đánh giá cuối); seed (giá trị khởi tạo ngẫu nhiên) khác nhau và sinh final sau LOCK.json. Chỉ thử ridge 0.01 và 1, chọn theo MAE (sai số tuyệt đối trung bình) trên selection. Lambda được chọn 0.01; MAE selection 0.008090 giây. Không huấn luyện bằng final. Bộ lọc lịch chỉ là bằng chứng trên mẫu đã thấy; không phải chứng chỉ topo. Lịch đủ điều kiện: ['noQ', 'q10', 'q12', 'q14', 'q16', 'author', 'q11_16', 'q13_18']. Ngoài miền n/mật độ đã khóa, điểm trùng hoặc sai hình học đơn vị: dự phòng noQ. Giới hạn miền trong model.json.

2. **Giữ DMT và luật PV?** Có. DMT là lý thuyết Morse rời rạc; PV là đỉnh giả. Engine được xuất nguyên mã commit (phiên bản ghi nhận trong Git) 259d4bf, biên dịch Release (cấu hình tối ưu) với GCC 13.3.0, C++23, -O3 -DNDEBUG -fopenmp, 4 luồng. Mã băm SHA-256 của binary (tệp thực thi) `a950a64ed2c4f9434bc589ad00dbe8bffbc923e599db9a8b0b6e3e64579eb566` trùng bản đóng băng trước lượt này. Không sửa matching (ghép cặp Morse), quotient, chọn PV, hỗ trợ witness hoặc cách dựng phức. Bộ chọn chỉ chuyển lịch sang CLI (giao diện dòng lệnh). Không dùng f2-stream/f2-rebuild; CONTINUE giữ trạng thái tạm dừng. Hồ sơ cũ được bảo toàn.

Đã đọc đủ 7 trang hi.pdf (hash/text trong hi-proposal.*). Các tỷ lệ tăng tốc, chứng chỉ an toàn, tham số và mô tả Gcurr=|A|−1 trong đó chưa được xác minh trên code này. G v2 cũ đã dự đoán chênh lệch thời gian. Code PV hiện tại đã thu nhiều vùng rồi cập nhật chung; không nhận gom lô là cải tiến mới. Không áp dụng chính sách đổi PV trong tài liệu. Chi phí đặc trưng hiện tại O(n²×11), không phải O(1), và được tính trong D.

3. **G chọn lịch nào?** Trong mọi danh sách, 2.0 là cutoff (ngưỡng kết thúc), không quotient ở đó. Mốc được tạo từ k nguyên chia 10. Có đúng 10 ứng viên, không quét tổ hợp khác.

| Bộ cuối | Lịch G | Trong miền | T̂ lịch (s) | G dự đoán (s) |
|---|---|---|---:|---:|
| final-uniform-n64-s2026101900 | [1.0, 2.0] | Có | 0.01615 | 0.00541 |
| final-clustered-n64-s2026101901 | [1.0, 2.0] | Có | 0.01647 | 0.00509 |
| final-uniform-n100-s2026101902 | [1.0, 2.0] | Có | 0.03769 | 0.02187 |
| final-clustered-n100-s2026101903 | [1.0, 2.0] | Có | 0.04245 | 0.01851 |
| final-uniform-n200-s2026101904 | [1.0, 2.0] | Có | 0.13946 | 0.14872 |
| final-clustered-n200-s2026101905 | [1.0, 2.0] | Có | 0.14749 | 0.14303 |

4. **So sánh tốc độ và RAM.** Lịch cố định (fixed) tốt nhất chọn trước final là `q12` = [1.2, 2.0]. A: không quotient; B: tác giả [1.2,1.7,2]; C: lịch cố định đó; D: G gồm dự đoán; E: chạy lại đúng lịch G bỏ dự đoán. F: Ripser chỉ kiểm PH, không chọn lịch và không công bố benchmark (đo đối chứng hiệu năng) tốc độ/RAM với Ripser.

Mỗi phương án có 1 warm-up (làm nóng) và 5 lượt đo, tuần tự, xoay và đảo thứ tự giữa lượt. Bảng dưới là median (trung vị) toàn tiến trình, giây; khoảng [min,max] là nhỏ nhất–lớn nhất. RSS (bộ nhớ thường trú) là peak (mức cao nhất) GNU time, KiB; không dùng peak theo cửa sổ verbose (log chi tiết) làm số đo toàn tiến trình. G dùng execv: Python nhường cùng tiến trình cho engine, không giữ thêm tiến trình Python chờ PH. Phép đo bắt đầu trước khởi động Python, gồm đọc, đặc trưng, dự đoán, ghi lựa chọn, PH, xuất và giải phóng. Diagnostic (lượt chẩn đoán) bật verbose riêng, không nằm trong bảng thời gian.

| Bộ | Phương án | PH | Trung vị [min,max] s | RSS trung vị [min,max] KiB |
|---|---|---|---:|---:|
| final-uniform-n64-s2026101900 | A | PASS | 0.02196 [0.02037,0.02227] | 6604 [6572,6732] |
| final-uniform-n64-s2026101900 | B | PASS | 0.01896 [0.01863,0.02798] | 5720 [5564,5876] |
| final-uniform-n64-s2026101900 | C | PASS | 0.01702 [0.01568,0.01902] | 5684 [5556,5724] |
| final-uniform-n64-s2026101900 | D | PASS | 0.09260 [0.08236,0.09696] | 16820 [16808,16888] |
| final-uniform-n64-s2026101900 | E | PASS | 0.01850 [0.01699,0.01920] | 5968 [5840,5996] |
| final-clustered-n64-s2026101901 | A | PASS | 0.02268 [0.01952,0.02410] | 6776 [6492,6852] |
| final-clustered-n64-s2026101901 | B | PASS | 0.01772 [0.01703,0.02239] | 5572 [5420,5640] |
| final-clustered-n64-s2026101901 | C | PASS | 0.01719 [0.01600,0.02313] | 5968 [5876,6016] |
| final-clustered-n64-s2026101901 | D | PASS | 0.09106 [0.08260,0.09181] | 16820 [16780,16888] |
| final-clustered-n64-s2026101901 | E | PASS | 0.01570 [0.01544,0.01736] | 5868 [5736,6036] |
| final-uniform-n100-s2026101902 | A | PASS | 0.04806 [0.04610,0.05500] | 10108 [10028,10140] |
| final-uniform-n100-s2026101902 | B | PASS | 0.03615 [0.03432,0.03978] | 7320 [7156,7356] |
| final-uniform-n100-s2026101902 | C | PASS | 0.03058 [0.02752,0.03235] | 6992 [6884,7140] |
| final-uniform-n100-s2026101902 | D | PASS | 0.11766 [0.11414,0.12779] | 16832 [16816,16944] |
| final-uniform-n100-s2026101902 | E | PASS | 0.03686 [0.03324,0.04067] | 8008 [7812,8088] |
| final-clustered-n100-s2026101903 | A | PASS | 0.04675 [0.04392,0.05644] | 10488 [9740,10616] |
| final-clustered-n100-s2026101903 | B | PASS | 0.03149 [0.03063,0.03688] | 6720 [6688,6792] |
| final-clustered-n100-s2026101903 | C | PASS | 0.03071 [0.02641,0.03578] | 7796 [7424,7908] |
| final-clustered-n100-s2026101903 | D | PASS | 0.11355 [0.10850,0.12835] | 16808 [16804,16844] |
| final-clustered-n100-s2026101903 | E | PASS | 0.03497 [0.03415,0.03633] | 7892 [7840,7972] |
| final-uniform-n200-s2026101904 | A | PASS | 0.27673 [0.26907,0.30506] | 43684 [43316,44100] |
| final-uniform-n200-s2026101904 | B | PASS | 0.15780 [0.14810,0.16999] | 19484 [19240,19640] |
| final-uniform-n200-s2026101904 | C | PASS | 0.13133 [0.11389,0.14888] | 20768 [20608,20980] |
| final-uniform-n200-s2026101904 | D | PASS | 0.26234 [0.25040,0.28082] | 17264 [16976,17348] |
| final-uniform-n200-s2026101904 | E | PASS | 0.11744 [0.11540,0.12886] | 17184 [17096,17408] |
| final-clustered-n200-s2026101905 | A | PASS | 0.29057 [0.28729,0.30419] | 44104 [43988,44232] |
| final-clustered-n200-s2026101905 | B | PASS | 0.17973 [0.17571,0.20641] | 25744 [25688,26200] |
| final-clustered-n200-s2026101905 | C | PASS | 0.13377 [0.12715,0.13787] | 19292 [19100,19492] |
| final-clustered-n200-s2026101905 | D | PASS | 0.26280 [0.24460,0.27259] | 22596 [22300,22740] |
| final-clustered-n200-s2026101905 | E | PASS | 0.13109 [0.12484,0.13143] | 22500 [22456,22708] |

Tỷ số thời gian X/D >1 nghĩa là G nhanh hơn X. Chỉ dùng giao các bộ mà cả X và D đều PASS; tổng hợp bằng trung bình nhân của tỷ số trung vị từng bộ. Không lấy ca sai PH làm thành tích phương pháp PH chính xác.

- A/D: 0.4757 trên 6/6 bộ; loại 0 bộ: [].

- B/D: 0.3344 trên 6/6 bộ; loại 0 bộ: [].

- C/D: 0.2921 trên 6/6 bộ; loại 0 bộ: [].

- E/D: 0.3009 trên 6/6 bộ; loại 0 bộ: [].

5. **Đúng, sai, chưa xác minh.** PASS = khớp kiểm chứng; PH_MISMATCH = sai PH; UNVERIFIED = chưa xác minh; TIMEOUT/MEMORY_LIMIT = vượt thời gian/RAM. So multiset (tập giữ số lần lặp) của mọi khoảng; C++ xuất endpoint (đầu mút) 17 chữ số, so Decimal chính xác giữa các lịch. Ripser 0.6.15, F2, khoảng cách Euclid, maxdim=1, thresh=2, không chọn tập điểm đại diện, không xuất đối chu trình; đối chiếu sau chuyển endpoint iCC sang float32 (số thực 32 bit), không thêm tolerance (sai số chấp nhận). H0 lấy từ MSF (rừng khung nhỏ nhất) thật do clearing (loại cạnh đã giải H0) của iCC tạo ở từng cửa sổ; nối endpoint bằng ledger (sổ theo dõi), kiểm số thành phần liên thông và so chính xác với phương án tham chiếu không quotient và Ripser. H0 được kiểm trong lượt chẩn đoán riêng của đúng lịch; mọi lượt đo H1 đều so với lượt chẩn đoán đó. Đây là kiểm H0 có điều kiện khoảng cách dương, alpha<1; không suy ra bảo toàn H1 tổng quát.

iCC lọc khoảng cách<epsilon, Ripser gồm equality. Survivor (khoảng còn sống) cuối có death=-1 ở iCC, +∞ ở Ripser; chỉ nói sống tới cutoff, không khẳng định sống mãi. Bỏ survivor tạm ở các cửa sổ chưa cuối. Không tự bỏ khoảng dài 0 ở parser: giữ mọi khoảng engine thực xuất; Ripser thường không xuất khoảng 0. Kernel (nhân tính toán) Morse bỏ cặp H1 có birth=death, parser (bộ đọc kết quả) không lọc thêm khoảng nào; số cặp dài 0 được kiểm riêng trong AUDIT.json. Đặc trưng kiểm điểm phân biệt, không có cạnh đúng cutoff 2 trong dữ liệu ngẫu nhiên này; các khác biệt cutoff/điểm trùng trong hồ sơ cũ vẫn tồn tại. Không nới sai số để đổi FAIL thành PASS.

Trên final, 30 cấu hình phương án/dataset (các alias vẫn được đo riêng); trạng thái từng phương án: {'A': {'PASS': 6}, 'B': {'PASS': 6}, 'C': {'PASS': 6}, 'D': {'PASS': 6}, 'E': {'PASS': 6}}. Tổng 180 lượt gồm làm nóng: {'PASS': 180}. Khảo sát 180 cấu hình lịch/dataset: {'PASS': 178, 'PH_MISMATCH': 2}. Bảng correctness.csv và survey-correctness.* giữ toàn bộ ca sai/chưa kiểm; final-correctness.json giữ chẩn đoán và H0.

Ca theo dõi 8 điểm tách khỏi tập huấn luyện/đánh giá cuối: [('noQ', 'PASS'), ('author', 'PASS'), ('[0.5,2]', 'PH_MISMATCH')]. Không phân tích lại nguyên nhân/thu nhỏ. Ca [0.5,2] nằm ngoài 10 ứng viên, dùng tái hiện lỗi đã biết. Các ca sai không làm dừng khảo sát các lịch khác, và không đổi lựa chọn final sau khi xem Ripser.

6. **Chi phí controller có bù lợi ích?** Chi phí D−E dưới đây là chênh lệch hai trung vị toàn tiến trình, có nhiễu lịch hệ điều hành; paired (ghép theo vòng) min/max nằm trong controller-overhead.csv. pre_exec_s chỉ đo phần bên trong Python trước exec, thiếu khởi động interpreter và bước ghi cuối, nên không thay thế D−E.

| Bộ | D−E (s) | Phần trước exec (s) |
|---|---:|---:|
| final-uniform-n64-s2026101900 | 0.074101 | 0.007206 |
| final-clustered-n64-s2026101901 | 0.075352 | 0.007030 |
| final-uniform-n100-s2026101902 | 0.080799 | 0.017674 |
| final-clustered-n100-s2026101903 | 0.078580 | 0.016071 |
| final-uniform-n200-s2026101904 | 0.144892 | 0.074855 |
| final-clustered-n200-s2026101905 | 0.131710 | 0.066162 |

Trung bình chênh lệch trung vị D−E: 0.097572s. C/E (giá trị lịch trước phí G): 0.970818760791853; C/D (sau đủ phí G): 0.29213224945133903. G chưa thắng lịch cố định tốt trên giao PASS; không có bằng chứng chi phí chọn lịch được bù so với C.

MAE dự đoán thời gian lịch đã chọn trên 6 bộ cuối, so với E: 0.008310s. Đây là phân tích sau khóa, không dùng đổi lịch.

Học và khảo sát không miễn phí: từ thời điểm bắt đầu bảo thủ tới khóa mô hình mất 756.4s, gồm đọc ngữ cảnh, triển khai, lỗi thiết lập, biên dịch, khảo sát và học. Các số phase (giai đoạn) dưới đây đo trong harness (chương trình điều phối) sau khi nạp thư viện Python, chưa gồm khởi động/nạp thư viện của chính harness. Mọi ca PH/Ripser bên trong vẫn có số đo toàn tiến trình riêng. Tổng thời gian các tiến trình đã ghi: 193.492s (gồm cả build/kiểm chứng/đo, không cộng nhầm vào thời gian G).

| Giai đoạn chuẩn bị | Thời gian harness (s) |
|---|---:|
| prepare | 0.045389 |
| build | 45.452941 |
| survey | 155.205046 |
| fit | 0.732547 |

7. **Điều đáng báo thầy và giới hạn.** Có một bộ chọn toàn lịch chạy được, mã nhỏ, dữ liệu độc lập trước/sau khóa, giữ nguyên thuật toán Morse/iCC và PV, đối chiếu H0–H1 riêng, đo đầy đủ phí điều khiển và giữ ca sai. Kết quả tốc độ và coverage (tỷ lệ bộ vượt kiểm PH) ở trên là kết quả thực, gồm cả bất lợi. Mô hình có thể chọn cùng một lịch cho nhiều hình học; khi đó chưa chứng minh lợi ích chọn lịch theo dữ liệu. Chưa có chứng chỉ an toàn quotient, chưa khắc phục upstream, chưa đánh giá H2 (đồng điều chiều 2), chưa chứng minh khái quát ngoài hai họ hình học/n64–200 và máy này. Không dùng thành tích phát triển thay cho final, không gọi đây là tái lập toàn bộ bài báo.

Tài nguyên: giới hạn process-group (nhóm tiến trình riêng của ca) 5.00 GiB, timeout (hạn ca) tối đa 300s hoặc thời gian còn lại trước mốc đóng thí nghiệm. Monitor (bộ theo dõi) khoảng 50ms, không phải giới hạn cgroup cứng. Không đổi swap/WSL, không chạy đồng thời các ca nặng. Không dừng tiến trình người dùng.

Lỗi công cụ có 2 traceback lưu riêng, cộng initial-tooling-errors.json. Mọi raw command (lệnh gốc), stdout/stderr, exit code (mã kết thúc), mã băm dữ liệu đầu vào trước/sau và RSS ở cases/. Không sửa log cũ. INDEX.md, REPRODUCE.md và gói kèm manifest SHA-256 phục vụ bàn giao.

![So sánh thời gian và RAM trên sáu bộ cuối](runs/selector-20261016T064700Z/performance.png)

Tất cả timing chính là lượt đo mới, cùng binary, chỉ H0–H1; không trộn H0–H2 của hồ sơ cũ. Tập cuối có một seed cho mỗi tổ hợp n/hình học, chưa đủ suy rộng. Mật độ được lấy từ khoảng cách gốc; các file PLAN/model/LOCK dùng số nguyên k, epsilon=k/10; xem candidates.csv để đọc lịch theo epsilon.

Hai cấu hình khảo sát sai đều trên `dev-clustered-n100-s2026101709`: `[0.8,2.0]` và `[0.8,1.2,1.6,2.0]`. Cùng bar H1 có birth=1.0620512302574743, death đúng=1.0899597931988616, death quotient=1.0898437955693845. Không sửa tolerance, chỉ dùng nhãn sai trên development để loại lịch trước khóa; không đổi lịch bằng đáp án final.

Kiểm gói dùng cùng ngân sách ban đầu: dựng lại source/CLI11 từ gói trùng SHA binary, 3/3 CTest, chạy lại 30 lượt A–E trên một input cuối đều PASS. Đây là tái hiện input cũ, không tăng số mẫu final. `measurement-check.json` kiểm GNU peak RSS giữ qua execv (75716 KiB trước/sau). AUDIT.json kiểm 13 điều kiện toàn vẹn và tách tập đều PASS.

Metadata epsilon của header mặc định chỉ in 6 chữ số được phục hồi từ chính argv float64 của lệnh để xác định cửa sổ. Không sửa endpoint PH 17 chữ số hoặc raw stdout; parser-input.log là bản dẫn xuất. Không thay tolerance.
