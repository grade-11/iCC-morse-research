# Hồ sơ S² n100 — tái lập và giới hạn

Phạm vi: 100 điểm trên S² đơn vị trong R³, seed20261007, dữ liệu tương đương tự sinh;
không phải CSV gốc của tác giả. Không kết luận hoàn tất Table I/giai đoạn 1.
Chạy các lệnh dưới đây từ `/home/alex_1709/projects/iCC-team`.

## Môi trường và đầu vào

Python đã dùng: `/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python`,
Python 3.12.3; ripser 0.6.15, numpy 2.5.4, scipy 1.18.1, scikit-learn 1.9.1,
persim 0.3.8. Freeze đầy đủ nằm trong reference cũ và `build-context/ripser-freeze.txt`
của dossier. Không sửa Python hệ thống hay `.venv-plot`.

C++: GCC 13.3.0, C++23, CMake 3.28.3, OpenMP, CLI11 v2.6.2;
build Release `-O3 -DNDEBUG`, tối đa 4 worker. Kiểm tra version/flags thực tế
trong `build-context/`; hash binary phụ thuộc môi trường build.

```bash
set -e
PY=/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python
TASK=experiments/pwph2026-s2-n100
printf '%s\n' \
 'e39f71c26b3f71c96860c04316b9eae7ec192891ebdd39c80bbebb2029f76b48  experiments/pwph2026-s2-n100/data/sphere_S2_R3_n100_seed20261007.csv' \
 '730a2aea8b946d2d04995e189647188115146bdcbf78f4024e8e6d219d275e72  experiments/pwph2026-s2-n100/benchmark_icc.py' \
 '25f1937999c41bf0780626006827ab36ddf7cc8e1574f655f1b45dbdaf43bbf4  experiments/pwph2026-s2-n100/runs/whole-run-monotonic-ma523yck.tar.gz' \
 'a949d7ad5bcc22062a32c6784b6851530414494159a9c326fb386adcef5d131a  experiments/pwph2026-s2-n100/runs/ripser-reference-noa_bg6w.tar.gz' \
 | sha256sum --check || exit 2
OUT=$(mktemp -d "$TASK/runs/n100-recheck-XXXXXX")
```

`evidence_run.py` luôn tạo thư mục con mới và lưu `command.txt`, exact argv/cwd,
UTC timestamps, monotonic elapsed, versions, hashes trước/sau trong `context.json`,
`stdout.log`, `stderr.log`, `exit-status.txt`, `report.json`. Nó từ chối ghi đè.
Exit status của lệnh con được giữ; timeout=124, lỗi capture/input=125.
`report.json` là JSON stdout của chương trình con, hoặc bản tóm tắt capture;
luôn kiểm tra thêm `exit-status.txt` và stderr, không đọc riêng chữ PASS.

## Đối chiếu hai archive cũ, không chạy PH

```bash
"$PY" "$TASK/evidence_run.py" --output "$OUT/archive-comparison" \
 --input "$TASK/compare_ripser_n100.py" --input "$TASK/benchmark_icc.py" \
 --input "$TASK/runs/whole-run-monotonic-ma523yck.tar.gz" \
 --input "$TASK/runs/ripser-reference-noa_bg6w.tar.gz" -- \
 "$PY" "$TASK/compare_ripser_n100.py" \
 --icc-archive "$TASK/runs/whole-run-monotonic-ma523yck.tar.gz" \
 --ripser-archive "$TASK/runs/ripser-reference-noa_bg6w.tar.gz" \
 --benchmark-script "$TASK/benchmark_icc.py"
```

Comparator exit 0=H1/H2 PASS, 1=mismatch, 2=evidence không hợp lệ.
Giữ nguyên pin parser SHA256; nếu khác, review/provenance trước khi cập nhật pin.
H0 trong comparator archive cũ vẫn `NOT_CHECKED`, vì log cũ không xuất H0.
Các hash binary/HEAD cũ là metadata được ghi trong context, không có binary
trong archive để xác thực độc lập. Không gán HEAD hiện tại cho benchmark cũ.

## Kiểm tra workflow và build riêng

```bash
"$PY" "$TASK/evidence_run.py" --output "$OUT/workflow-capture" -- \
 "$PY" "$TASK/test_workflow_n100.py" --output "$OUT/workflow-cases"
"$PY" "$TASK/evidence_run.py" --output "$OUT/configure" --timeout 120 -- \
 cmake -S . -B "$OUT/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
 -DFETCHCONTENT_SOURCE_DIR_CLI11="$PWD/build/_deps/cli11-src" \
 -DFETCHCONTENT_FULLY_DISCONNECTED=ON
"$PY" "$TASK/evidence_run.py" --output "$OUT/build-capture" --timeout 600 -- \
 cmake --build "$OUT/build" --parallel 4
"$PY" "$TASK/evidence_run.py" --output "$OUT/ctest" -- \
 ctest --test-dir "$OUT/build" --output-on-failure
```

Dừng khi một lệnh build/configure lỗi; không dùng nhầm `build/iCC` cũ.
CLI11 source trên máy này đã có ở `build/_deps/cli11-src`; nếu thiếu, cần cung cấp
CLI11 2.6.2 đã xác minh hoặc cài dependency trong phạm vi tác vụ, rồi ghi lại nguồn.
Không chạy build/test đồng thời với phép đo.

## Xác nhận binary mới (tuần tự, không benchmark lặp)

```bash
"$PY" "$TASK/evidence_run.py" --output "$OUT/validation-capture" --timeout 600 \
 --input "$OUT/build/iCC" --input "$TASK/validate_n100.py" -- \
 "$PY" "$TASK/validate_n100.py" --binary "$OUT/build/iCC" --output "$OUT/validation"
```

Lệnh chạy baseline và piecewise n100, mỗi mode verbose/nonverbose, tối đa 4 luồng;
fixture H0 dùng 1 luồng. Timeout từng n100 120 s, fixture 30 s.
Nó kiểm tra Decimal multisets H1/H2 cũ/mới, float32 với Ripser, invariant
PV/active/cofacets/labels/FFI, bảo toàn counters và H0 từ đường MSF thực của iCC.
Fixture kiểm chứng lặp death/multiplicity, nhiều survivors, cạnh đúng cutoff,
điểm trùng và khác biệt quy ước; raw Ripser fixture NPZ được giữ.
Đây là kiểm tra hồi quy n100 cố định; không phải bộ kiểm thử PH tổng quát.

## H0, dimension và counters

`max-dim=3` là chiều simplex; output dimension=2/3 tương ứng H1/H2.
`DimensionFrame::buildInitialActiveFacetIndex` gọi MSF clearing trong
`SimplexUtility::getActiveEdgeIndexHashTable`. Cạnh được union bị loại khỏi active
edge mask và chính là các death H0. Chẩn đoán verbose `[h0-msf]`/`[h0-msf-pair]`
đọc phần bù mask đó, đếm components trên active labels, không chạy MST độc lập.
Mỗi cửa sổ báo H0 graph hiện tại, không nối các bar sau quotient với bar gốc.
Không thêm H0 vào hợp đồng output nonverbose.

Riêng n100: 99 death H0 tối đa 0.47270135747272135 < 1.2, nên đã xác định toàn bộ
finite H0 trước quotient đầu; giữ 1 survivor theo quy ước unreduced H0. H0 baseline
và cửa sổ piecewise đầu khớp Decimal, khớp Ripser sau float32. Cửa sổ sau vẫn 1
component. Không suy rộng sang dữ liệu còn nhiều components khi quotient bắt đầu.

Hai khác biệt có sẵn: iCC dùng cạnh có distance < epsilon, Ripser gồm equality;
piecewise bỏ cạnh distance=0, baseline giữ chúng. Ripser bỏ zero-length H0 bars,
trong khi MSF diagnostic giữ chúng nếu xuất hiện. n100 không có các trường hợp
biên này. Không bỏ multiplicities/zero-length để cưỡng ép PASS fixture.

`[pv-filters]`: mỗi đơn vị là một support được duyệt theo thứ tự khám phá đảo.
Thứ tự first-rejection: chứa PV → chứa nhãn protected → overlap với PV đã chọn
trong lượt → cardinality (rỗng hoặc >=64) → diameter > alpha*eps_hi → separation.
`label_clean` chỉ qua hai bộ lọc đầu; `eligible_at_visit` qua tất cả điều kiện tại
thời điểm xét, và được chọn ngay. Không phải tập ứng viên hợp lệ độc lập thứ tự.
Các rejection counters là loại đầu tiên gặp; không mô tả số supports thỏa từng
predicate độc lập. `raw = rejected_pv + rejected_protected + label_clean`;
`label_clean = rejected_overlap + rejected_cardinality + rejected_diameter +
rejected_separation + selected`. Cửa sổ cuối không quotient: supports/filters N/A.

## Tổng hợp profiling và đóng gói

Dossier đã có `tables/whole-run-paired.csv`, `whole-run-summary.csv`,
`windows-old-verbose.csv`, `windows-old-summary.csv`, `windows-new-validation.csv`,
`timings-long.csv`. Các bảng mới đọc lại bytes archive và so với thư mục raw,
không chỉ tin tên batch. Có thể tái tạo bảng bằng:

```bash
# summarize_n100 cần thư mục raw cũ cùng tên archive và kết quả validation-v2.
# Với OUT mới, tạo symlink validation-v2 trỏ tới validation của chính OUT.
ln -s validation "$OUT/validation-v2"
"$PY" "$TASK/evidence_run.py" --output "$OUT/profiling-capture" -- \
 "$PY" "$TASK/summarize_n100.py" --dossier "$OUT"
```

Không cần chạy lại benchmark/reference để dùng hồ sơ đã có. Nếu thật sự cần tạo
reference mới, bọc `run_ripser_n100.py` bằng `evidence_run.py`; runner tự tạo thư mục
và archive mới. `benchmark_icc.py run` giữ cấu hình binary cũ `./build/iCC`; không
trỏ ngầm nó sang binary khác hoặc trộn số đo với batch cũ.

Thống kê loại warmup round00, giữ mọi round01..10 kể cả realtime bất thường;
dùng monotonic elapsed. Worker elapsed sums không phải CPU/wall time và nằm trong
timer cofacet rebuild, không cộng vào timer ngoài. RSS theo cửa sổ là VmHWM đã
reset `/proc/self/clear_refs=5`; chỉ dùng khi memory_status=ok. Nó khác whole-run
peak từ GNU Time, và verbose có thêm FFI/stdout/diagnostics.

Archive cuối gồm báo cáo, nhật ký nối tiếp, scripts/source/patch, raw archives nhỏ,
logs, tests và build context; loại build/venv/cache/bản thảo/Zone.Identifier.
Kiểm tra mọi entry trong `MANIFEST.sha256.json` và hash `.tar.gz.sha256` cạnh archive.
Không tự đồng bộ Library/ChatGPT; nhật ký cập nhật nằm trong gói bàn giao.
