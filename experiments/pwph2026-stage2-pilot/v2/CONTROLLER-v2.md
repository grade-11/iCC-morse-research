# Thiết kế G v2 (viết trước fit và đánh giá)

Trạng thái s sau matching ở epsilon, trước chọn PV: active labels/PVs và
representatives, dependency supports/protected labels, số top cofacets C,
epsilon và chiều simplex D. Hai hành động dùng **cùng trạng thái, cùng h**:

- CONTINUE: không chọn/commit PV, giữ nguyên WindowState, dựng lại cùng nhãn
  từ 0 tới h theo kiến trúc PwPH; không giữ interface matching cũ. Chi phí dựng
  lại được đo đầy đủ, không giả rằng đây là incremental matching miễn phí.
- QUOTIENT: chọn/commit PV theo đúng luật cũ, rồi dựng lại tới cùng h.

h là mốc kế tiếp trên lưới hữu hạn đã công bố. Những checkpoint bị DEFER vẫn
phải xử lý matching và xuất PH; G điều khiển **có quotient tại mốc đó hay không**.
Trì hoãn không đổi nhãn/PV, nhưng bookkeeping vẫn chỉ xuất finite death ở cửa
sổ hiện tại và giữ survivor cuối. Chạy mask toàn DEFER phải khớp baseline.

`G_hat = (R_hat(s,0,h)+F_hat(s,0,h)) - [Q_hat(s)+R_hat(s,1,h)+F_hat(s,1,h)]`.
Mọi số đo/prediction là **giây**:

- Q: processPwph + ensureBinomialCapacity, gồm selection và commit, chỉ khi Q.
- R: prepareWindow, ctor DimensionFrame và advance của cửa sổ tới h.
- F: thời gian cửa sổ tới h trước quyết định tiếp - R. Gồm toàn interface,
  bookkeeping, giải phóng và xuất PH; không gọi là kernel Morse thuần.
- Controller/telemetry/decision logging nằm trong whole-run; không cộng timer
  witness/interface lồng trong R/F. Chi phí checkpoint hiện tại đã xảy ra như
  nhau trước cả hai hành động nên không là lợi ích tiết kiệm tương lai.

Thu mẫu offline bằng replay tiền tố, thay đúng một hành động Q/DEFER rồi dừng
sau h, trước quyết định ở h. Xác nhận state hash và trace tương đương ở epsilon;
ghi chi phí chạy cả hai nhánh riêng, không giấu vào benchmark online. Đánh giá
sai số G và đúng dấu trên validation có nhãn chi phí hai nhánh thực sự.

Estimator thử nghiệm: ridge tuyến tính cho R,F và peak RSS, dùng feature rẻ
chuẩn hóa theo đơn vị đã chốt (C, A, epsilon→h, tỷ lệ PV, support count và tương
tác action). Q là ridge theo số supports/entries. Fit trên development, kiểm
validation; không fit hay chọn threshold bằng final. Hệ số xuất file text có
hash. Online chỉ hiệu chỉnh một hệ số thời gian chung bằng quan sát quá khứ
(EMA), không chạy nhánh phản thực. Tối đa vài biến thể, không quét lớn.

QUOTIENT nếu G lớn hơn `max(0.0005 s, margin*Ccontinue)`; mặc định margin 5%.
DEFER khi G âm/không đủ biên. Đặc trưng ngoài miền calibration dùng fallback
DEFER được công bố, không có nghĩa đảm bảo PH nếu đã quotient trước đó.
RAM là ràng buộc riêng: RSS dự báo phải dưới 70% ngân sách để chừa headroom;
nếu chỉ Q khả thi có thể Q vì RAM dù G âm, log rõ lý do; nếu cả hai không khả
thi thì dừng ca có trạng thái riêng. Không cộng byte vào giây.

Đây là controller G thử nghiệm, không định lý tối ưu, không phải sửa PH, không
phát triển chọn PV theo chi phí. Các cấu hình vẫn cần PH validation trước đo.
Thiết kế/feature cụ thể và lock sẽ được bổ sung bằng số liệu development/validation.

## Feature và ranh giới số đo đang dùng

Đặt `w=C/100000`, `z=log(h/epsilon)`, `a=A/100`, `p=PV/A`, action q∈{0,1}.
Vector tám phần tử của R/F/RSS là
`[1,w,w*z,a^(D+1)*z,p*w,q*w,q*raw_supports/100,q]`.
Q dùng `[1,raw_supports/100,support_entries/1000]`.
Ridge chuẩn hóa bằng mean/std của development; intercept không phạt. Khi xuất
file text, biến đổi coefficients về feature thô để C++ không dùng thư viện ML.
Hai ridge thử .01 và 1, chọn bằng normalized G MAE trên validation, báo cả
MAE giây và đúng dấu (toàn bộ/|G thực|>0.5 ms). Không bỏ nhãn phản thực bất lợi.

Trong nhánh DEFER, chi phí reset support hiện tại được cộng vào target
F_continue; Q target là timer bao processPwph + ensureBinomialCapacity, không
đếm lại các công đoạn bên trong. Phí hủy phần support-container còn lại ở cuối
vòng Q chưa tách được trong timer này; cùng telemetry/controller/logging, nó
vẫn nằm trong whole-run. Nhãn G_actual trong bảng calibration vì thế là hiệu
chi phí các phase đã đo, không phải nhãn counterfactual whole-process chính xác.
Đây là giới hạn của estimator, nhất là ca rất ngắn; biên 0.5ms là dự phòng
thử nghiệm, không phải upper bound đã chứng minh. RSS target là VmHWM của tiền tố tới h, không phải riêng RSS của
cửa sổ. Prediction RSS lấy max với HWM đã quan sát. Đây là ràng buộc dự báo,
còn giới hạn RSS cây process của harness là cơ chế bảo vệ thực, gần đúng 50 ms.

Online scale ban đầu 1, cập nhật nhân theo EMA .2 từ observed/predicted cost
của cửa sổ vừa xong, clamp scale trong [.25,4]. Không cập nhật model bằng
nhánh chưa chạy. Domain fallback khi A hoặc C >1.5 lần cực đại calibration.
Log `[g]` ghi dự báo và hành động; `[state]` ở cửa sổ sau ghi cost thực của nhánh
đã chọn. Counterfactual của final không được gọi là đã đo.

G có thể chọn khác khi bật verbose do timing/FFI. Vì vậy validation chạy G
không verbose trước, rồi **replay đúng các hành động đã quan sát** để lấy MSF
và PH diagnostic; kiểm hash trạng thái và multiset với lượt không verbose.
Không lấy PH của một trajectory verbose khác để xác nhận trajectory benchmark.

## Các biến thể và đánh giá trước final

Bốn chính sách đã chạy: g0/g5 trên [1.2,1.5,1.7,2] với margin 0/5%;
gd5 trên [1.2,1.4,1.6,1.8,2] margin5%; gs5 trên [1.2,2] margin5%.
Margin0 vẫn giữ ngưỡng tuyệt đối 0.5ms. Biến thể gs5 có giả thuyết giảm phí
checkpoint, thu thêm nhãn hai nhánh tới cùng h=2 thay vì chỉ ngoại suy từ h
ngắn. PLAN-AMENDMENT-2 được ghi trước batch này, trước mọi dataset final.

Đã có 320 cặp replay tương đương: development 192+48 và validation64+16,
tức 640 process C++ thu nhãn chi phí offline. Không tính chúng là PH PASS.
Mỗi chiều có 40 cặp validation khi tính cả chân trời cuối. Ridge .01 được
chọn trước final so với 1; ablation bỏ tương tác PV và q*raw làm MAE G xấu hơn
ở model ban đầu, dù độ đúng dấu của ablation có thể tốt hơn. Giữ feature đầy
đủ theo tiêu chí MAE đã khai báo; không lựa bằng final.

Model gs5 fit development cộng cặp terminal, vẫn lambda .01 đã chọn trước.
So **cùng** 40 nhãn validation: D2 model ban đầu MAE8.859ms, global9.694ms;
D3 lần lượt673.003ms và657.837ms. Riêng terminal8 cặp: D2 25.042→21.298ms,
D3 2685.372→2290.567ms. Có đánh đổi giữa các chân trời, không gọi đây là cải
thiện dự báo đồng đều. Bảng estimator-common-validation.json giữ từng mẫu.

Đúng dấu cũng chưa tốt hơn baseline rất đơn giản: cả40 cặp D3 có G thực dương,
luôn dự báo dương đạt100%, trong khi model ban đầu90%, global75%. D2 baseline
luôn dương87.5%, hai model80%/77.5%. Confusion/class balance trong
G-prediction-summary.json. Những số này không là bằng chứng dự báo online tốt.
Hai cặp có sai số lớn nhất được lặp thêm hai lần; G thực D2 71.98–73.47ms,
D3 3842.84–3947.74ms. Sai số model lớn hơn nhiều mức nhiễu quan sát này.

RAM được thử thực ở cùng state hash trên multiscale n32 D3: ngân sách dự báo
39.313MiB làm DEFER (memory-reject-quotient), 512MiB làm QUOTIENT (positive-G).
Hai process exit0, dừng sau cửa sổ kế; không dùng làm kết quả PH toàn cutoff.
RSS thực chỉ vài MiB, nên đây là kiểm chức năng phản ứng với **dự báo** RAM,
không chứng minh tiết kiệm bộ nhớ hoặc nhu cầu áp lực RAM thực. Monitor cây
process áp cùng hai cap, xấp xỉ50ms; GNU time giữ peak RSS đầy đủ hơn sample.
Kết quả tìm nhánh âm không có candidate phù hợp trước đó vẫn được giữ riêng.

Fixed tốt và G được chọn theo chiều/nhóm n≤64/n≥100 từ development+validation,
ưu tiên coverage PH, rồi H0 đã kiểm, rồi geometric mean thời gian. Lựa chọn
này ưu tiên độ đúng quan sát, không nói quotient của một ca mới là an toàn.
LOCK-v2.json khóa timestamp, binary, code chạy, model, map model, lưới và final
plan trước khi sinh seed30–37. Final không thay estimator/threshold/lịch.

Benchmark G luôn không verbose. Nếu state/action khác trajectory đã kiểm vì
thời gian cập nhật scale, replay mask thực và kiểm state identity/MSF ledger
cho trajectory đó ngoài lượt đo. Chi phí diagnostic này ghi riêng; online
không được nhìn nhánh kia. Nếu PH đo sai, giữ lỗi và dừng đo cấu hình đó.

Lựa chọn khóa thực tế là sparse/gs5 ở cả bốn nhóm. Vì gs5 chỉ có một quyết
định trước final (epsilon1.2, h2), đánh giá final này không đo lợi ích của EMA
cho quyết định tiếp theo. EMA đã hoạt động trong các ứng viên nhiều checkpoint
ở development/validation; không gán lợi ích final cho cơ chế cập nhật đó.
