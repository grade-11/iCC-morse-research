# Truy nguyên upstream và thử CONTINUE H0–H1

Run hiện tại nằm trong `CURRENT-RUN.txt`; bắt đầu đọc `INDEX.md` rồi `REPORT.md`
ở run đó. Đây là workflow độc lập với hai lượt pilot trước. Ngân sách4h và
20phút cuối dành cho đóng gói được lưu một lần trong `budget.json`; không tạo lại
ngân sách khi resume. Tất cả experiment giữ stdout/stderr, argv, input/binary
SHA, exit code, whole-process monotonic time, GNU time RSS và trạng thái giới hạn.
Không push; các file chưa commit của người dùng được giữ nguyên.

`launch.py STAGE` chạy một script trong thư mục này và ghi outer-capture. Các
stage đã hoàn tất được nhận diện bằng kết quả lưu; case có context nhưng thiếu
exit code phải được kiểm tra trước khi retry, không ghi đè. Không chạy đồng thời
benchmark với công việc nặng khác. `resume.py` kiểm tra deadline/boot, khóa runner
và marker; không tự chạy tiếp khi đã hoàn tất hoặc qua budget.

Các stage chính: `prepare_versions`, `instrumentation`, `h0_checks`,
`stream_checks`, `counterexample8`, `survey`, `n200_dossier`,
`provenance_summary`, `grid_foundation`, `g_boundary`, `n300_h2`,
`bounds_h0_audit`, `evidence_audit`. Main matrix nằm trong `matrix.py`.
Không chạy lại stage đã hoàn tất chỉ để tạo thêm log. N200 cũ đã được kiểm tra
hash/manifest và copy chọn lọc, không gán timing của D cho S.

New C++ mode: `--continue-engine f2-stream` dùng explicit F2 reduction khác
native Morse; chỉ H0/H1, giữ pivots/live births thật giữa checkpoints. Đối chứng
`f2-rebuild` dùng cùng kernel và reset mỗi checkpoint. Nhiều checkpoint phải có
`--quotient-mask` toàn0. H2, quotient, G/v1 và trace bị từ chối. Đọc
`CONTINUE-DESIGN.md`; đây chưa phải native Morse CONTINUE/migration sau quotient.

Tái hiện từ gói bàn giao (đặt `RUN` là thư mục vừa giải nén):

```bash
cmake -S "$RUN/source" -B /tmp/icc-continue-repro -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DCMAKE_DISABLE_FIND_PACKAGE_CLI11=ON \
  -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
  -DFETCHCONTENT_SOURCE_DIR_CLI11="$RUN/third-party/CLI11"
cmake --build /tmp/icc-continue-repro -j4
ctest --test-dir /tmp/icc-continue-repro --output-on-failure
/tmp/icc-continue-repro/iCC --tool ph --file-name "$RUN/inputs/two-squares-s2-8.csv" \
  --max-dim 2 --max-eps 2 --threads 4
/tmp/icc-continue-repro/iCC --tool piecewise --file-name "$RUN/inputs/two-squares-s2-8.csv" \
  --max-dim 2 --eps-breaks 0.5 2 --pv-cap-scale 0.85 --threads 4
/tmp/icc-continue-repro/iCC --tool piecewise --file-name "$RUN/inputs/two-squares-s2-8.csv" \
  --max-dim 2 --eps-breaks 0.5 2 --pv-cap-scale 0.85 --threads 4 \
  --continue-engine f2-stream --quotient-mask 0 --verbose
```

Để tái hiện upstream nguyên trạng, đổi `-S` thành `$RUN/versions/U` và dùng
build directory khác; U17 có patch **chỉ precision output** lưu riêng. Không đưa
các flag mới của nhóm vào U. Các `cases/*/command.txt` là lệnh gốc chính xác;
đường dẫn tuyệt đối chỉ tới môi trường lúc đo. Khi chuyển máy, đổi input/binary
prefix theo `versions.json`, giữ nguyên nội dung input/config; với n200 cũ dùng
`n200-reused/copy-map.json`. Reference cache có config đầy đủ và checksum.
Không dùng thời gian tái hiện có verbose/sanitizer thay cho benchmark đã lưu.
