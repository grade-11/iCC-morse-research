# Truy nguyên PH — kết quả đã kiểm chứng, chưa sửa thuật toán

Hai phản ví dụ cũ đã tái hiện trên binary dd92867 ở 1 và 4 luồng: ca 11 điểm
mất một H2; ca 10 điểm mất hai H1 và một H2. Baseline khớp Ripser float32.
Max-dim 2 được kiểm riêng: H1 hai ca khớp, nhưng ca 10 điểm H0_NOT_VERIFIED
vì graph chưa liên thông trước quotient đầu. Không suy ra max-dim 3 đúng từ đó.

Trace mới chỉ thêm diagnostics có flag, không sửa matching/selection/weights.
`root-cause-oracle.json` và `oracle.py` kiểm tra độc lập từng cửa sổ bằng:

1. Duyệt toàn bộ joint witness assignments trên các PV thực từ trace; so **từng
   simplex và từng trọng số float64** với enumerator C++. Không có sai khác ở
   cả hai ca, cả 1/4 luồng.
2. Rút gọn ma trận biên F2 tường minh trên chính phức đã trace, rồi so các bar
   thuộc cửa sổ hiện tại với raw PH C++. Tất cả H1/H2 ở mọi cửa sổ đều khớp.
3. So phức witness đó với phức VR gốc trên cùng khoảng; sai khác PH vẫn còn.

Ca 11 điểm: chưa quotient được ở 1.2; sau 1.5 có PV
`{3,4,5,6,10}` và `{1,2,7,8,9}`, chỉ còn ba active labels.
Tại cửa sổ cuối 2.0, phức trên ba nhãn không thể chứa H2; VR gốc có bar
H2 `(1.7508756119083391,1.7757032130821528)`.

Ca 10 điểm: sau 1.2 PV `{0,1,3,4,5}`; sau 1.4 thêm PV `{2,6,7,8,9}`,
chỉ còn hai nhãn. Tại 2.0 VR gốc có hai H1 và một H2 mới (birth >1.7),
còn phức trên hai nhãn không thể biểu diễn chúng. Các cửa sổ tới 1.7 chưa có
các bar tương lai này, nên chỉ nhìn live features lúc quotient không phát hiện.

Đây là bằng chứng mất thông tin ở phép biểu diễn bằng nhãn/quotient và phức
witness của **hai ca thu nhỏ**, không phải lỗi parser, race 1/4 luồng, hay một
lỗi số học trong enumerator/matching đã được chứng minh. Chưa chứng minh rằng
sáu ca n100 đầy đủ cùng nguyên nhân; không suy diễn từ một chỉ số FFI.

Bản thảo PwPH2026 §§V–VI nêu rõ quotient CW khác quotient simplicial, đường
kính <=alpha*epsilon làm vùng cục bộ contractible tại lúc tạo, và bảo vệ critical
vertices nhằm giữ các birth đang sống. §VII thừa nhận ảnh hưởng interval chiều
thấp; kết luận top-dimensional preservation cần lựa chọn alpha phù hợp, không
cung cấp một điều kiện kiểm tra đủ cho mọi dữ liệu/lịch ở alpha=0.85.
DMT2025 §IV.C/Algorithm4 xây lại trên virtual labels; §VI cũng thừa nhận có
thể mất short features liên hệ hình học với vùng quotient. Không dùng định lý
quotient CW để kết luận phép chiếu simplicial đang dùng bảo toàn toàn bộ module
persistence tương lai. Hai PDF được đọc cục bộ, hash và trang trong `papers/`.

**Chưa có bản sửa PH có chứng minh.** Không đổi alpha, kernel, protection,
expected bars hoặc parser để ép PASS; cũng không gọi bỏ quotient hoàn toàn là
đã sửa PwPH. Lượt này tiếp tục controller trên từng cấu hình được đối chiếu,
giữ PH_MISMATCH/H0_NOT_VERIFIED riêng. Cần thiết kế phép biểu diễn giữ các
identification/chain maps hoặc điều kiện đủ mạnh hơn, vượt một bản vá hẹp.


Bổ sung H0 (commit riêng): `h0_ledger.py` nối endpoint từ MSF **thật** qua
các cửa sổ theo lập luận connected-region quotient, kiểm conservation, so
multiset baseline/Ripser. Sáu unit tests đạt; audit 576 cấu hình ban đầu đạt
576/576, gồm 18 H0_NOT_VERIFIED cũ. Evidence cũ không đổi, bổ sung ở
`screen-h0-enhanced-*.json`. Phạm vi/giả thiết trong H0-LEDGER.md; đây là
bổ sung kiểm chứng H0 có điều kiện, không sửa biểu diễn H1/H2 nền.

## Hai phản ví dụ mới từ development

Không thấy sai H1 trong batch mới đã chạy ở D2 hoặc D3; H1 cũ ở ca10 điểm
vẫn còn và vẫn là giới hạn bắt buộc ghi nhận. Chọn hai lỗi H2 fixed khác hình
học, trước final; ddmin có mục tiêu, không tuyên bố tối tiểu:

- uniform S2 n64 seed20261013, lịch[.9,1.3,1.7,2]: còn11 điểm,106 ứng viên,
  23.8s. Trước cửa sổ1.7→2 chỉ còn3 nhãn; thiếu H2
  (1.7138536425200717,1.7915925075328596).
- polar S2 n100 seed20261018, lịch[1.2,1.7,2]: còn14 điểm,140 ứng viên,
  34.6s. Trước1.7→2 còn6 nhãn; thiếu H2
  (1.6049548916725145,1.7530400614293653). Sáu nhãn riêng nó chưa chứng minh
  bất khả biểu diễn H2, nên không dùng lập luận cardinality cho ca này.

Cả hai được kiểm lại baseline/Ripser, lưu trace và input riêng. Oracle brute
joint-witness kiểm **mọi simplex/weight** của các phức nhỏ đúng bằng C++; reduction
F2 độc lập cũng khớp mọi PH C++ từng cửa sổ. Phức quotient khác VR gốc ở H2.
Như vậy trên hai input **đã thu nhỏ**, sai khác nằm ở biểu diễn quotient/witness,
không phải một lỗi matching/weight đã phát hiện. Chưa có kết luận nguyên nhân
cho toàn bộ input64/100 điểm từ việc thu nhỏ này. new-root-cause-oracle.json
lưu đầy đủ; không đổi expected endpoints hay tolerance.
