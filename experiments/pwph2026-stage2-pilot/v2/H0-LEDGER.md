# H0 qua quotient: phạm vi có điều kiện, từ MSF thật

Đây là bổ sung kiểm chứng/bookkeeping endpoint H0, không sửa H1/H2 hoặc quy tắc
chọn PV. Không áp dụng mù phép cộng các local bars.

Điều kiện: mọi khoảng cách giữa hai điểm khác nhau dương; alpha∈(0,1); dùng
đúng quy tắc disjoint/frozen PV đang cài đặt. Mỗi support được chọn chỉ chứa
original active labels, diameter<=alpha*epsilon<epsilon; do đó PV nằm trong
một clique liên thông tại lúc hình thành, và tiếp tục liên thông về sau.

Với đồ thị 1-skeleton tại t>=epsilon, thay một tập đỉnh liên thông bằng một
nhãn không đổi số component. Cạnh giữa các nhãn có trọng số min khoảng cách
giữa witness sets, nên có cạnh đúng khi có ít nhất một cạnh gốc nối hai tập.
Như vậy component và thời điểm hợp nhất component được giữ qua bước quotient.
Kết luận này chỉ nói H0; không suy ra simplicial witness complex giữ H1/H2.

MSF của cửa sổ i được lấy từ chính phần bù active-edge mask mà iCC clearing
đã tạo, không từ một union-find thay thế. Tất cả vertices sinh ở 0. Cửa sổ đầu
giữ finite deaths thật; cửa sổ sau chỉ giữ finite deaths >=eps_lo. Với strict
cutoff `<`, death bằng eps_lo chưa xảy ra ở cửa sổ trước nên thuộc cửa sổ sau.
Không giữ lại các local MSF bars cũ có death<eps_lo. Survivor cuối có multiplicity
bằng components cuối, không cộng survivors tạm.

Mỗi bước kiểm `new_mergers = components_previous - components_current`, bước
đầu là n_original-components_first. Kiểm tổng finite+survivors=n_original;
so toàn multiset chính xác với MSF baseline và Ripser float32. Nếu điều kiện,
conservation hoặc reference không đạt: không ghi PASS. Khi baseline có death0,
giữ KNOWN_DIFFERENCE và không dùng chứng minh khoảng cách dương. Alpha=1 hoặc
thay chính sách quotient cần audit riêng; không coi là đã hỗ trợ tổng quát.

`test_h0_ledger.py` kiểm multiplicities, nhiều survivors, connected-region
quotient khi graph gốc còn nhiều components, bỏ old local bars, và các đường
reject. `h0-ledger-audit.json` so với raw MSF/Ripser trên toàn batch có sẵn;
mọi H0_NOT_VERIFIED cũ vẫn được giữ trong evidence ban đầu, kết quả ledger mới
là lớp kiểm chứng bổ sung riêng có source/hash. Không sửa CSV, log hoặc expected
bars để đạt PASS. Bộ fixture cutoff equality/điểm trùng của n100 giữ nguyên.
