# PwPH tiếp nối v2 — lượt tự chủ 4 giờ

Đọc `BAO-CAO-TRIET.md`; truy nguyên chi tiết trong `PH-ROOT-CAUSE.md`, thiết kế
và giới hạn G trong `CONTROLLER-v2.md`. `CURRENT-RUN.txt` trỏ run mới. Pilot v1,
CSV, archive và báo cáo cũ không được sửa. Đây chưa là hoàn tất giai đoạn 2/3.

`budget.json` giữ start/deadline monotonic và boot ID; không đặt lại 4 giờ khi
resume. Sau experiment-stop chỉ đóng gói trong 20 phút còn lại. Mỗi C++/Ripser
process tối đa 300 s, 4 threads, tuần tự. RSS cây process theo /proc là giám sát
gần đúng, không hard cgroup. Chỉ kill process group do tác vụ tạo.

```bash
MPLCONFIGDIR=/tmp/icc-v2-mpl OPENBLAS_NUM_THREADS=4 OMP_THREAD_LIMIT=4 \
 /home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python \
 experiments/pwph2026-stage2-pilot/v2/resume.py
```

Resume tự lưu stdout/stderr/argv/timestamps/exit code, chỉ tiếp tục stage thiếu;
mỗi stage lấy lock riêng cho run. Không chạy hai runner. Khi binary/source/model
đã khóa đổi thì không trộn với batch cũ. Evidence trường hợp interrupted chưa
có exit status phải được giữ và xem xét trước retry có tên/lý do riêng.

`lab.py` tái dùng logger, MSF và comparators đã review của v1/n100. Parser chỉ
tổng quát hóa dimension/cutoff, không đổi endpoints/tolerance/multiplicities.
G diagnostic replay actual mask và kiểm trạng thái. H0 dùng ledger từ MSF thật từng cửa sổ, có kiểm conservation và so toàn bộ
multiset với baseline/Ripser. Phạm vi có điều kiện: khoảng cách dương, alpha<1,
luật disjoint/frozen PV hiện tại; xem H0-LEDGER.md. Audit bổ sung giữ nguyên mọi
H0_NOT_VERIFIED ban đầu. Không suy ra bảo toàn H1/H2.

`PLAN.json` khai báo dataset trước batch; `PLAN-AMENDMENT*.json` lưu việc thu gọn
controller candidates và biến thể chân trời cuối trước final. `LOCK-v2.json` khóa model/params/grid/fixed.
Seed07–11 đều là dữ liệu cũ. Final seed30–37 chỉ sinh sau lock; seed38–39 dự trữ.

`study.py` có stages plan/calibrate/fit/screen/lock/final. Calibration replay
hai tiền tố có hash/trạng thái quan sát bằng nhau, thay một hành động rồi đo tới
cùng h; đây là chi phí offline riêng. `ablate.py`, `stability.py`, `shrink_new.py`
là các kiểm tra có mục tiêu. Không benchmark Ripser hoặc Python wrapper speed.

`summarize.py --run RUN` chỉ đọc evidence và tạo bảng. `performance.csv` chỉ
nhận nhóm đủ năm vòng; warmup riêng trong measurements. `correctness.csv` giữ
mọi ca sai/chưa xác minh. Alias fixed-original trùng best-fixed dùng chung số
đo, không gọi là thêm process/dataset. So reference/G theo round, >1 lợi cho G.

Archive cuối gồm source/patch/tests/scripts, data, models, raw logs, papers text
và checksum; không build/venv/cache. Build Release/CLI11/configure commands trong
run. Muốn build lại từ source gói, tạo thư mục build mới và batch mới nếu hash
binary khác; không gán số đo cũ cho binary mới.

`audit.py --run RUN` kiểm provenance/lock/completeness. Sau commit cục bộ,
`package.py --run RUN` tạo một archive và SHA256, tự kiểm từng member với manifest.
Archive chứa thêm gói pilot cũ để truy xuất input/dependencies, không đổi gói đó.

Gói có source CLI11 đã dùng ở `third-party/CLI11`, kèm commit/hash và license;
không cần fetch mạng để build C++. Từ thư mục RUN đã giải nén, có thể chạy:

```bash
cmake -S "$PWD/source" -B /tmp/icc-v2-rebuild -DCMAKE_BUILD_TYPE=Release \
 -DBUILD_TESTING=ON -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
 -DFETCHCONTENT_SOURCE_DIR_CLI11="$PWD/third-party/CLI11"
cmake --build /tmp/icc-v2-rebuild -j4
ctest --test-dir /tmp/icc-v2-rebuild --output-on-failure
/tmp/icc-v2-rebuild/iCC --tool piecewise \
 --file-name data/uniform_s2-n32-s20261030.csv --max-dim 2 --threads 4 \
 --eps-breaks 1.2 2.0 --pv-cap-scale .85 --g-model models/global-d2.txt --g-margin .05
```

Python/numpy/ripser/scipy/matplotlib versions trong environment-context.json;
logger ghi cả phiên bản từng reference. Các đường dẫn raw evidence giữ đúng
cwd gốc. Resume là cho repo/run gốc trên cùng boot; trên máy khác cần tạo run
và budget mới rõ ràng, không giả tiếp tục ngân sách monotonic cũ.

Build offline từ source/CLI11 đóng gói đã được kiểm: 3/3 CTest đạt và SHA256
binary trùng binary benchmark (source-reproduction.json). Resume đã smoke-test,
exit0, không thêm process thực nghiệm; xem resume-smoke.json. CLOCK-CHECK.json
ghi chênh lệch giữa UTC descriptive và elapsed monotonic; ngân sách dùng mốc
monotonic cũ nghiêm hơn, không đặt lại hoặc kéo dài deadline.
