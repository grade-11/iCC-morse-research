# CONTINUE giữ trạng thái — thiết kế trước thử nghiệm

Phân biệt ba đường: native iCC hiện tại dựng lại 0→h; tham chiếu F2 thực sự giữ
reduced columns qua checkpoint; và native CONTINUE cần thiết kế thêm. Không gọi
tham chiếu F2 là native Morse continuation, không dùng tốc độ Python để benchmark.
Không sửa quotient trong lượt thay đổi này.

Hiện trạng native: `PwphPipeline::run` tạo `DimensionFrame` trong block mỗi
window; frame chết trước `WindowState::invalidateCurrentWindow`. Chỉ active labels
và PV representatives sống qua window. Ngay trong frame,
`releaseInterfaceWorkspace` reset mate arrays, full cofacet lookup và raw support;
`releaseActiveFacetIndex` xóa mask; `advance` giải phóng facet list cũ, chuyển
cofacets sang chiều sau. Khi kết thúc max-dim, không còn interface của chiều thấp.
`FullCofacetIndex` mượn `cofacets.data()`; vector tăng capacity sẽ làm con trỏ
đó hết hiệu lực nếu giữ nguyên hash functor. Chỉ giữ matchList không giữ đủ state.

Thứ tự: mỗi list dùng (weight,bindex) tăng; matching duyệt facets giảm. Paper2025
Algorithm3 mô tả incremental traversal bên trong một interface; không đồng nghĩa
append filtration windows. Paper2026 §VI(e) nói rõ rebuild từ0. New weights ở
strict checkpoint là >=lo; tie groups không tách bởi <epsilon, nhưng cần giữ
quy tắc bindex và face/coface order. Sau quotient, labels bị loại/thêm và simplex
không còn cùng định danh: cần reset đầy đủ hoặc chain-map/vineyard migration có
kiểm chứng, không tái dùng vị trí rank cũ.

Một native CONTINUE cần giữ theo mọi chiều: facet/cofacet storage có stable
handles; complete lookup, active/clearing/compression masks; mates và reduction
paths/support provenance; danh sách live births; MSF union-find và edge prefix;
frontier enumeration. H0 prefix-Kruskal nối tiếp được khi label set không đổi.
H1/H2 cần cập nhật cả old live facets khi new cofaces tới. Apparent pairs và
compression của interface dưới ảnh hưởng active rows của interface trên; không
chỉ append arrays hoặc chạy matching cho new facets. Mọi checkpoint chỉ là
live_at_checkpoint, không phát một class là vô hạn vĩnh viễn trước cutoff cuối.

Phản ví dụ tối thiểu cho việc đóng băng survivors: square bốn điểm, cạnh1,
diagonal sqrt(2). H1 sinh1, còn sống tại1.1, chết sqrt(2) trước1.5. Giữ kết quả
unmatched cũ mà không cho new coface ghép với old birth sẽ sai. Test phải còn
xác nhận rank/pivot migration, clearing qua interface và tie groups; square đơn
lẻ không chứng minh một continuation implementation đúng tổng quát.

Prototype `stream_reference.py` dùng standard explicit F2 boundary reduction:
append simplex theo (weight,dimension,bindex), faces trước cofaces khi tie; giữ
reduced columns, pivot→owner, live births. Chỉ reduce cột mới, dùng cột đã reduce
để cập nhật, không recompute các cột cũ. Tuy nhiên nó pre-enumerate tới cutoff
cuối (ghi rõ chi phí), không dùng native Morse kernel/clearing và không phục vụ
G. Quotient gọi vào sẽ bị từ chối rõ. Đây là reference/proof-of-concept cho
no-quotient, chưa phải tính năng CONTINUE production của iCC.

Điểm dừng native: chưa có API append an toàn, chưa chứng minh cách extension
đồng thời các interface của matcher reverse-facet và cách remap sau quotient.
Không vá hàng loạt lifespan/rank/matching/quotient cùng lúc rồi dùng PASS nhỏ
để suy ra an toàn. Chỉ mở mode native sau khi có thiết kế state migration và
regressions tương ứng; không loại phí rebuild khỏi thời gian để gọi incremental.

## Mode C++ đã hiện thực và kiểm chứng trong lượt này

`--continue-engine f2-stream` hiện giữ state thật cho H0–H1, với engine F2 khác
native Morse. `f2-rebuild` là đối chứng cùng engine reset mỗi checkpoint.
Default `native-rebuild` không đổi thuật toán; so hồi quy riêng với binary259d.
Chỉ max-simplex2, nhiều checkpoint phải có mask toàn0. Quotient, G/v1, H2,
trace và early-stop bị từ chối rõ trước chạy pipeline. Chưa có migration native
sau quotient; không gọi đây là CONTINUE đầy đủ của controller G.

Khác Python reference pre-enumeration, C++ chỉ sort các cạnh dưới cutoff một
lần mỗi engine. Adjacency được mở rộng theo nhóm edge-weight; new triangle được
sinh khi cạnh cuối của nó tới. Cả nhóm edges xử lý trước triangles cùng trọng
số; triangle tie dùng colex/bindex. IDs của vertices/edges ổn định theo thứ tự
filtration, vectors di chuyển không làm đổi IDs; pivots lưu reduced column
bằng ID, không mượn con trỏ vào buffer cofacets. Chỉ columns mới được reduce,
old reduced columns được đọc để XOR. Old live birth được đóng khi một new
column có pivot đó; survivor intermediate không được giữ là infinity vĩnh viễn.

Reduction C++ hiện serial; không gán speedup với native4-thread chỉ cho việc giữ
state. So f2-stream với f2-rebuild cùng kernel mới là đối chứng trực tiếp cho
reuse. Cả việc sort edges, rebuild control, output và teardown đều tính vào
whole-process runtime. Các counter new/total edge/triangle chứng minh khối
lượng prefix không bị reduce lại; full prefix state hash khớp engine fresh.

Kiểm gồm square/rectangle sống qua checkpoint, nhiều partitions kể cảk/10,
zero distances, equidistant octahedron, graph disconnected và cutoff equality.
H0 của mode mới so multiset với actual MSF native PH, H1 so native/Ripser;
H0 zero bars được giữ trong output thực, Ripser bỏ zero nên ghi KNOWN riêng.
H2 C++ bị từ chối, dù Python oracle hỗ trợ để đọc phản ví dụ nhỏ.

Metadata epsilon mặc định của native chỉ in6 chữ số, endpoints17 chữ số. Với
argv .3, binary64 thực là .29999999999999999. Parser mới tái dựng **chỉ bounds**
bằng format17g của argv; raw logs/endpoints giữ nguyên, không tolerance.
Rectangle .18×.24 có death đúng binary64 .3 là regression cho trường hợp này.
So lại cả340 cấu hình nguồn gốc: multiset không đổi. Upstream U nguyên trạng
in6 chữ số vẫn chỉ được so ở printed precision, không gọi là full64-bit PASS.
