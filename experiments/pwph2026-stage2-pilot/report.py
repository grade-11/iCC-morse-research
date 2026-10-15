#!/usr/bin/env python3
"""Build the Vietnamese handoff from saved evidence (no experiments)."""
import argparse
import csv
import json
from pathlib import Path
from collections import Counter
from pilot import BASE,sha,write

def report(run):
    lock=json.loads((run/'LOCK.json').read_text());binary_hash=lock['binary_sha256']
    rows=list(csv.DictReader((run/'summary.csv').open()))
    ratios=list(csv.DictReader((run/'paired-ratios.csv').open()))
    validations=[json.loads(p.read_text()) for p in run.glob('validate-*.json')]
    validations=[v for v in validations if v['binary_sha256']==binary_hash]
    counts=Counter(v['status'] for v in validations)
    write(run/'correctness-overview.json',dict(counts=dict(counts),configurations=[dict(
        case=v['case'],status=v['status'],H0=v.get('H0'),H1_H2_exact=v.get('H1_H2_exact'),
        H1_H2_ripser_float32=v.get('H1_H2_ripser_float32'),windows=v.get('windows')) for v in validations]))
    def row(n,seed,m):return next(x for x in rows if int(x['n'])==n and int(x['seed'])==seed and x['method']==m)
    def fmt(r):return f"{float(r['median_s']):.4f} [{float(r['min_s']):.4f}–{float(r['max_s']):.4f}]"
    def mem(r):return f"{float(r['median_rss_kib'])/1024:.2f} [{float(r['min_rss_kib'])/1024:.2f}–{float(r['max_rss_kib'])/1024:.2f}]"
    lines=['# Báo cáo pilot giai đoạn 2 cho Triết — 16/10/2026 (+07)', '',
        '**Đã code và chạy heuristic online sau checkpoint đầu 1.2; chưa hoàn tất toàn bộ giai đoạn 2, chưa làm giai đoạn 3.** '
        'Đối chứng chính là PwPH lịch cố định đã khóa `[1.2,1.7,2.0]`. '
        'Cần đọc kết quả theo từng seed; không lấy mức nhanh hơn iCC không thu gọn làm đóng góp mới của adaptive.', '',
        'Đã khảo sát năm lịch cố định trên ba seed phát triển n100; khóa lịch tốt đủ điều kiện PH '
        'bằng trung bình nhân của trung vị thời gian từng seed. Heuristic target=0.35 được giữ nguyên. '
        'Sau đó đánh giá hai seed độc lập, rồi n200 cùng nhóm seed, cùng binary. '
        'Không mở rộng n400/dữ liệu thực tế; không push/merge/công bố.', '',
        'Baseline = cách đối chứng; checkpoint = mốc xử lý; quotient = thu gọn vùng thành nhãn đại diện; '
        'rebuild = dựng lại phức; peak RSS = bộ nhớ cư trú lớn nhất; regression = kiểm tra sửa code có làm hỏng '
        'kết quả đã đúng; ablation = thử riêng từng thành phần.', '',
        '## Độ đúng', '',
        f"Binary mới có **{counts.get('PASS',0)} cấu hình PASS H0–H2**, **{counts.get('PH_MISMATCH',0)} PH_MISMATCH** "
        f"(các trạng thái khác: {dict((k,v) for k,v in counts.items() if k not in ('PASS','PH_MISMATCH')) or 'không'}). "
        'Mỗi cấu hình được kiểm tra riêng; mọi vòng đo còn so H1/H2 với baseline và Ripser. '
        'H0 dùng MSF thật của iCC trong diagnostic riêng, chỉ nhận PASS khi graph gốc đã liên thông trước quotient đầu, '
        'khớp baseline chính xác và Ripser sau float32. Không cộng local H0 của các quotient.', '',
        'So giữa các iCC bằng multiset Decimal của endpoints đầy đủ đã in, giữ multiplicities và survivor cuối. '
        'Với Ripser 0.6.15 chiếu endpoints iCC sang float32, không tolerance. '
        'Chỉ bỏ survivor tạm ở các cửa sổ chưa cuối, theo actual schedule; parser verbose và compact khớp. '
        'Không tuyên bố đúng float64 với Ripser hoặc chứng minh PH tổng quát.', '',
        'Regression n100 và CTest 3/3 đạt; kiểm tra harness 3/3 đạt, gồm timeout, lỗi exit7, memory limit giả lập, '
        'resume và parser nhiều cửa sổ. Hai fixture cạnh đúng cutoff/điểm trùng vẫn là **KNOWN_DIFFERENCE**: '
        'iCC strict `<`, Ripser gồm equality; piecewise bỏ cạnh 0 và Ripser không xuất zero-length H0. Không đổi quy ước.', '',
        '## Thời gian và RAM so đối chứng chính', '',
        'Mỗi cấu hình đủ PH: một warmup + năm vòng đo, thứ tự xoay/đảo, từng process tuần tự, 4 luồng. '
        'Whole-run không verbose dùng monotonic clock và GNU Time; gồm mọi chi phí controller/telemetry/logging/quotient/rebuild. '
        'Bảng dưới là seed độc lập; mỗi ô là **trung vị [min–max]**. Tỷ số theo cặp = fixed/adaptive; >1 có lợi cho adaptive.', '',
        '| n / seed | Fixed (s) | Adaptive (s) | Tỷ số thời gian theo cặp | Fixed RSS (MiB) | Adaptive RSS (MiB) |',
        '|---|---:|---:|---:|---:|---:|']
    for n in (100,200):
        for seed in (20261010,20261011):
            f,a=row(n,seed,'n100'),row(n,seed,'adaptive')
            q=next(x for x in ratios if int(x['n'])==n and int(x['seed'])==seed and x['fixed']=='n100')
            ratio=f"{float(q['median_time_ratio']):.3f} [{float(q['min_time_ratio']):.3f}–{float(q['max_time_ratio']):.3f}]"
            lines.append(f'| {n} / {seed} | {fmt(f)} | {fmt(a)} | {ratio} | {mem(f)} | {mem(a)} |')
    adaptive_configs=[v for v in validations if v['method']=='adaptive']
    same=sum([w[1] for w in v['windows']]==['1.2','1.7','2'] for v in adaptive_configs)
    lines += ['', f"Trong {len(adaptive_configs)} ca adaptive, {same} ca chọn cùng lịch cố định đã khóa. "
        'Chênh lệch thời gian ở các ca trùng lịch không thể quy cho thay đổi mốc. '
        'Riêng n100 seed độc lập 20261011, adaptive chậm hơn khoảng 5,3% theo tỷ số hai trung vị '
        'nhưng giảm peak RSS khoảng 12,5%; đây là đánh đổi RAM/thời gian, chưa là cải thiện cả hai. '
        'Pilot này chưa đủ bằng chứng để khẳng định heuristic thắng lịch PwPH cố định cấu hình tốt.', '',
        'Đủ số đo từng seed phát triển/độc lập, min–max và tỷ số RAM theo cặp trong '
        '`summary.csv`, `paired-ratios.csv`; từng vòng và warmup trong `measurements.csv`. '
        'Không trộn với benchmark n100 cũ. Lịch đối chứng chuyển nguyên từ phần phát triển n100 sang n200, '
        'chưa chọn lại lịch tốt nhất riêng cho n200. Lịch tốt nhất theo từng ca nhìn lại chỉ là tham khảo hậu nghiệm; '
        'không thay lịch khóa bằng lịch thắng trên held-out.', '',
        '## Controller đã chọn thế nào', '',
        '`r = 1 − active_after/active_before`; `width = clamp(0.5*r/0.35, 0.2, 0.5)` epsilon. '
        'Lấy mốc xa nhất trong lưới tăng không vượt current+width, hoặc mốc kế tiếp để bảo đảm tiến lên, kết thúc 2.0. '
        'Nhãn giảm nhiều thì cho khoảng sau rộng hơn; ít giảm thì xử lý gần hơn để có dịp quotient tiếp. '
        'Chỉ dùng trạng thái hiện tại, không dùng PH/thời gian tương lai; không restart từ đầu.', '']
    examples=[]
    for p in (run/'cases').glob('validate-S2_R3_n100_seed20261011-adaptive-*/stdout.log'):
        examples=[x for x in p.read_text().splitlines() if x.startswith('[adaptive]')]
        if examples:break
    if examples:
        lines += ['Ví dụ seed độc lập n100/20261011 (log nonverbose, đơn vị ms và µs):', '', '```text',*examples,'```','']
    schedules={v['case']:v.get('windows') for v in validations if v['method']=='adaptive'}
    write(run/'adaptive-actual-schedules.json',schedules)
    lines += ['Đây là **heuristic dựa trên mức giảm nhãn**, chưa hiện thực hàm lợi ích G trong AC_PwPH. '
        'Cofacets và timing được ghi để giải thích/đánh giá, chưa tham gia quy tắc; RAM là ràng buộc giám sát bên ngoài. '
        'Không bật FFI/verbose để lấy tín hiệu điều khiển. Giữ nguyên kernel Morse, witness weights, điều kiện/thứ tự PV '
        'và quotient của từng cửa sổ; chỉ thay mốc tiếp theo. Log quyết định hoãn tới sau vòng xử lý, '
        'record cửa sổ sau là kết quả quan sát thực sau quyết định trước; không có dự báo.', '',
        '## Kết quả âm và lỗi', '',
        'Ba lịch có checkpoint 1.4/1.5 sai H2 ở n100 seed20261007/20261008 trên cả binary gốc và mới. '
        'Ví dụ `[1.2,1.5,2.0]`, seed07: birth H2 0.97696823381778153 thành 0.93668824329311717, '
        'death giữ 1.6602624057539093. Sai khác vẫn tồn tại sau float32. Loại khỏi bảng hiệu năng hợp lệ, không ép PASS.', '']
    for p in sorted((run/'shrunk').glob('*/result.json')):
        c=json.loads(p.read_text());lines.append(f"- Ca {c['seed']}/{c['method']}: tái lập, giảm 100→{c['retained_n']} điểm, {c['attempts']} lần thử, {c['elapsed_s']:.1f} s; giữ `{p.parent.relative_to(run)}/retained.csv` và toàn bộ input/log. Đây là ca đã thu nhỏ, không tuyên bố tối tiểu toàn cục.")
    lines += ['', 'Ripser mới trên cả hai ca thu nhỏ khớp baseline H1/H2; PwPH mất H2 (11 điểm), '
        'hoặc mất cả H1/H2 (10 điểm), không phải sai số float32. H0 các ca thu nhỏ không dùng để kết luận hiệu năng.', '',
        'Audit chỉ mục lần đầu exit1 vì kiểm hash cấu hình là duy nhất, quên các vòng lặp chủ ý. '
        'Đã sửa khóa thành tên ca/vòng + hash cấu hình, chạy lại audit đọc dữ liệu; giữ log lỗi và source v1. '
        'Không chạy lại phép đo PH/timing vì lỗi audit này.', '',
        'Không sửa lỗi PH nền trong lượt thử lịch này. Bộ đối chứng n100 cũ không bị reset/ghi đè. '
        'Lần tạo nhánh đầu exit128 do .git chỉ đọc; thực hiện lại qua quyền sandbox đã được phép và thành công. '
        'Probe môi trường phát hiện không có psutil: harness dùng /proc, không cài package. '
        'Matplotlib từng cảnh báo cache không ghi được; chuyển MPLCONFIGDIR sang /tmp cho các đợt sau. '
        'TIMEOUT/MEMORY_LIMIT giả lập chỉ thuộc unit test, không gọi là lỗi thực nghiệm.', '',
        '## Giới hạn và điều cần trao đổi với thầy', '',
        'Pilot chỉ thích nghi sau checkpoint đầu **1.2**, cutoff **2.0**, alpha **0.85**, simplex max-dim **3**, F2; '
        'S² đơn vị trong R³ tự sinh, không phải dữ liệu gốc tác giả. Generator '
        '`dSphere_gen.dSphere_uniformRandom(dimension=3, r1=r2=1, trials=0, selection=min)`, NumPy seed đã ghi, '
        'hash CSV/source và phiên bản package trong datasets/environment/context. Trials=0 khớp chính xác tọa độ n100 đã review; '
        'không thay phân phối hay sửa bộ sinh.', '',
        'Cần trao đổi về điều kiện bảo toàn H2 khi thay lịch, H0 nối qua quotient khi chưa liên thông, '
        'và thiết kế estimator G/ablation trước khi tuyên bố đóng góp. Chưa đánh giá nhiều phần cứng, dữ liệu khó, '
        'ngân sách RAM khác, hoặc ý nghĩa thống kê của chênh lệch nhỏ. Số test nhiều không là chứng minh tổng quát.', '',
        'Giới hạn RAM cây tiến trình 6,37 GiB (min(8 GiB,70% MemAvailable lúc bắt đầu)), '
        'giám sát /proc khoảng 50 ms, **gần đúng**, không hard cgroup limit; mỗi process thử tối đa 300 s. '
        'Không thay swap/WSL hoặc kill process ngoài nhiệm vụ. Peak RSS đo cho toàn process C++ bằng GNU Time, '
        'khác tổng RSS lấy mẫu để giám sát. Không cộng các timer lồng nhau.', '',
        'CLI auto-generated lịch đã được đọc: dùng `--eps-interval-count`/`--eps-interval-scale`, '
        'nhưng không ép được cả mốc đầu 1.2 và cutoff 2.0; không đưa vào bảng cùng điều kiện. '
        'Ripser chỉ kiểm PH, không benchmark tốc độ wrapper.', '',
        '## Bàn giao và tiếp tục', '',
        f'Run: `{run}`. `PROGRESS.md`, `index.jsonl`, `LOCK.json`, '
        '`correctness-overview.json` và raw `cases/` là điểm tiếp tục; không cần đọc từng log để lấy kết luận. '
        'Source/patch, dữ liệu, scripts, log, checksum được gói trong **một** archive '
        f'`{run.name}-handoff.tar.gz` cạnh báo cáo này; kiểm SHA256 cạnh archive và manifest bên trong. '
        'Không chứa build/venv/cache; build context và hash binary vẫn được lưu. Commit cục bộ được ghi trong `git-final.json`.', '',
        'Lệnh resume có tự lưu log, kiểm hash và bỏ qua các ca hoàn tất:', '',
        '```bash',
        '/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python experiments/pwph2026-stage2-pilot/resume.py',
        '```', '',
        'Resume chỉ trong ngân sách 6 giờ ban đầu và phải không còn runner khác. Khi run đã đủ ca, '
        'không tự benchmark lại. Sau khi hết ngân sách, nghiên cứu tiếp cần lượt mới, giữ nguyên evidence này.']
    text='\n'.join(lines)+'\n'
    (BASE/'BAO-CAO-TRIET.md').write_text(text);(run/'BAO-CAO-TRIET.md').write_text(text)
    print(json.dumps(dict(report=str(BASE/'BAO-CAO-TRIET.md'),counts=dict(counts)),ensure_ascii=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();report(a.run)
