#!/usr/bin/env python3
"""Summaries only: this file never runs PH or changes a chosen schedule."""
from collections import Counter
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import statistics as st

def read(p):return json.loads(Path(p).read_text())
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def table(p,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def gm(x):return math.exp(st.mean(math.log(v) for v in x)) if x else None

def summarize(run):
    rows=read(run/'final-timings.json');choices=read(run/'choices.json');lock=read(run/'LOCK.json');model=read(run/'model.json')
    plan=read(run/'PLAN.json');budget=read(run/'budget.json');survey=read(run/'survey-correctness.json')
    perf=[];correct=[];overhead=[];pairs=[]
    ids=[x['id'] for x in plan['datasets'] if x['split']=='final'];by={x['id']:x for x in plan['datasets']}
    for id in ids:
        for method in 'ABCDE':
            allrows=[r for r in rows if r['dataset']==id and r['method']==method];rs=[r for r in allrows if not r['warmup']]
            counts=Counter(r['status'] for r in allrows)
            status='PASS' if len(allrows)==6 and counts['PASS']==6 else ('PH_MISMATCH' if counts['PH_MISMATCH'] else 'UNVERIFIED')
            r=dict(dataset=id,geometry=by[id]['geometry'],n=by[id]['n'],method=method,schedule=allrows[0]['schedule'],PH_status=status,
                measured_repeats=len(rs),warmups=sum(x['warmup'] for x in allrows),
                median_s=st.median(r['elapsed_s'] for r in rs),min_s=min(r['elapsed_s'] for r in rs),max_s=max(r['elapsed_s'] for r in rs),
                median_rss_kib=st.median(r['rss_kib'] for r in rs),min_rss_kib=min(r['rss_kib'] for r in rs),max_rss_kib=max(r['rss_kib'] for r in rs))
            perf.append(r);correct.append(dict(dataset=id,method=method,schedule=r['schedule'],status=status,repeat_statuses=dict(counts)))
        p={r['method']:r for r in perf if r['dataset']==id}
        gr=[r for r in rows if r['dataset']==id and r['method']=='D' and not r['warmup']]
        er={r['round']:r for r in rows if r['dataset']==id and r['method']=='E' and not r['warmup']}
        diffs=[r['elapsed_s']-er[r['round']]['elapsed_s'] for r in gr]
        overhead.append(dict(dataset=id,schedule=p['D']['schedule'],G_median_s=p['D']['median_s'],replay_median_s=p['E']['median_s'],
                             difference_of_medians_s=p['D']['median_s']-p['E']['median_s'],paired_median_difference_s=st.median(diffs),paired_min_difference_s=min(diffs),paired_max_difference_s=max(diffs),
                             G_pre_exec_median_s=st.median(r['controller_pre_exec_s'] for r in gr),PH_status=p['D']['PH_status']))
        for method in 'ABCE':
            if p[method]['PH_status']!='PASS' or p['D']['PH_status']!='PASS':continue
            mr={r['round']:r for r in rows if r['dataset']==id and r['method']==method and not r['warmup']}
            ratios=[mr[r['round']]['elapsed_s']/r['elapsed_s'] for r in gr]
            pairs.append(dict(dataset=id,comparison=method+'/D',median_runtime_ratio=p[method]['median_s']/p['D']['median_s'],
                              paired_ratio_median=st.median(ratios),paired_ratio_min=min(ratios),paired_ratio_max=max(ratios)))
    comparisons={}
    for m in 'ABCE':
        x=[r for r in pairs if r['comparison']==m+'/D'];included=[r['dataset'] for r in x]
        comparisons[m+'/D']=dict(datasets=len(x),geometric_mean_median_runtime_ratio=gm([r['median_runtime_ratio'] for r in x]),
             geometric_mean_paired_median_ratio=gm([r['paired_ratio_median'] for r in x]),included=included,
             excluded=[dict(dataset=id,reasons={r['method']:r['PH_status'] for r in perf if r['dataset']==id and r['method'] in (m,'D') and r['PH_status']!='PASS'}) for id in ids if id not in included])
    saved_before_controller=[]
    for id in ids:
        p={r['method']:r for r in perf if r['dataset']==id}
        if p['C']['PH_status']==p['E']['PH_status']=='PASS':saved_before_controller.append(p['C']['median_s']/p['E']['median_s'])
    inventory=[]
    for p in sorted((run/'cases').glob('*/context.json')):
        c=read(p);inventory.append(dict(case=c['case'],status=c.get('status','INCOMPLETE'),exit_status=c.get('exit_status'),elapsed_s=c.get('elapsed_monotonic_s'),path=str(p.parent)))
    final_counts={m:dict(Counter(r['status'] for r in correct if r['method']==m)) for m in 'ABCDE'}
    summary=dict(created_at=datetime.now(timezone.utc).isoformat(),final_datasets=len(ids),final_configurations=len(correct),final_counts=final_counts,
                 final_repeat_counts=dict(Counter(r['status'] for r in rows)),survey_counts=dict(Counter(r['status'] for r in survey)),comparisons=comparisons,
                 fixed_over_replay_geomean=gm(saved_before_controller),fixed_schedule=lock['fixed_schedule'],choice_counts=dict(Counter(x['schedule'] for x in choices)),
                 mean_controller_difference_of_medians_s=st.mean(x['difference_of_medians_s'] for x in overhead),
                 inventory_statuses=dict(Counter(r['status'] for r in inventory)),raw_processes=len(inventory),
                 offline_phase_costs={s:read(run/(s+'-phase.json')) for s in ['prepare','build','survey','fit'] if (run/(s+'-phase.json')).exists()},
                 total_logged_process_elapsed_s=sum(x['elapsed_s'] or 0 for x in inventory),
                 tooling_errors=[p.name for p in run.glob('tooling-error-*.json')])
    summary['selected_schedule_prediction_MAE_s']=st.mean(abs(x['predictions_s'][x['schedule']]-next(r['median_s'] for r in perf if r['dataset']==x['dataset'] and r['method']=='E')) for x in choices)
    write(run/'SUMMARY.json',summary);table(run/'performance.csv',perf);table(run/'correctness.csv',correct);table(run/'controller-overhead.csv',overhead);table(run/'paired-ratios.csv',pairs);table(run/'inventory.csv',inventory)
    labels={'A':'Không quotient','B':'Lịch tác giả','C':'Fixed tốt nhất trên dev/chọn cấu hình','D':'G, đủ chi phí','E':'Chạy lại lịch G, bỏ dự đoán'}
    lines=['# Báo cáo G chọn lịch trước khi chạy — 16/10/2026','',
      'Đã triển khai, biên dịch, khảo sát hữu hạn, khóa mô hình trước khi sinh dữ liệu cuối và đo riêng toàn bộ chi phí G. '+
      f"Tập cuối có {len(ids)} bộ độc lập, H0–H1 (đồng điều bền vững chiều 0 và 1), trường F2 gồm hai phần tử; G đạt kiểm PH {final_counts['D'].get('PASS',0)}/{len(ids)} bộ. G chọn [1.0,2.0] cho cả sáu bộ, nhưng chậm hơn lịch cố định tốt trung bình nhân {1/comparisons['C/D']['geometric_mean_median_runtime_ratio']:.2f} lần khi tính đủ chi phí. Đây là kết quả bất lợi, chưa có bảo đảm PH cho mọi dữ liệu.",'',
      f"Ngân sách cố định bắt đầu **{budget['started_at']}**, dừng mở thí nghiệm **{budget['experiments_stop']}**, hạn **{budget['deadline']}**; không đặt lại khi lỗi. Thời điểm lập bảng {summary['created_at']}. Repo/HEAD đầu lượt trong initial-context.json; commit bàn giao trong git-final.json. Không push, không reset/clean.",'',
      '1. **G chọn lịch dựa trên gì?** Quotient là phép lấy phức thương bằng các đỉnh giả PV. Bộ chọn đọc tọa độ trước PH, lấy n và mật độ cạnh ở các ngưỡng k/10 (k nguyên từ 8 đến 18), tính vector 10 đặc trưng: hệ số chặn, log(n/100), mật độ cạnh tại 0.8, số mốc quotient, có quotient hay không, mốc quotient đầu/2, tổng mật độ cạnh ở các mốc quotient, hai tương tác với log(n) và mật độ 0.8, cùng proxy biên Σdε(1−dε)(2−ε). Proxy là đại lượng thay thế từ khoảng cách gốc; không phải số witness, PV hay cofacet tương lai đã biết. Cofacet là đơn hình kề trên một mặt; witness là tổ hợp đỉnh gốc đại diện dùng xác định đơn hình và trọng số.',
      '',
      'Hồi quy ridge (bình phương tối thiểu có phạt hệ số) dùng 10 hệ số cho log thời gian toàn tiến trình của engine (bộ tính toán Morse/iCC). Lấy hàm mũ để dự đoán giây cho từng lịch, G(x,L)=T̂(noQ)−T̂(L). Các dự đoán này cùng phạm vi đọc dữ liệu, Morse/iCC tới 2, xuất PH và giải phóng của engine; phần Python chọn lịch là chi phí chung cộng thêm của phương án D, đo trọn bên ngoài và báo riêng, chưa được hồi quy thành dự đoán lợi ích ròng. Đây là giới hạn của G hiện tại: G dự đoán dương chưa đủ kết luận tổng thời gian D thắng A. Chọn lịch có thời gian dự đoán thấp nhất trong tập đủ kiểm tra phát triển/chọn cấu hình. Không gọi đây là controller (bộ điều khiển) thích nghi sau mỗi quotient. Không dựng các checkpoint (mốc kiểm tra) ngoài lịch để quyết định bỏ qua.',
      '',
      f"Có 12 bộ development (phát triển), 6 bộ selection (chọn cấu hình), 6 bộ final (đánh giá cuối); seed (giá trị khởi tạo ngẫu nhiên) khác nhau và sinh final sau LOCK.json. Chỉ thử ridge 0.01 và 1, chọn theo MAE (sai số tuyệt đối trung bình) trên selection. Lambda được chọn {model['ridge']}; MAE selection {model['selection_MAE_s']:.6f} giây. Không huấn luyện bằng final. Bộ lọc lịch chỉ là bằng chứng trên mẫu đã thấy; không phải chứng chỉ topo. Lịch đủ điều kiện: {model['eligible_schedules']}. Ngoài miền n/mật độ đã khóa, điểm trùng hoặc sai hình học đơn vị: dự phòng noQ. Giới hạn miền trong model.json.",
      '',
      '2. **Giữ DMT và luật PV?** Có. DMT là lý thuyết Morse rời rạc; PV là đỉnh giả. Engine được xuất nguyên mã commit (phiên bản ghi nhận trong Git) 259d4bf, biên dịch Release (cấu hình tối ưu) với GCC 13.3.0, C++23, -O3 -DNDEBUG -fopenmp, 4 luồng. Mã băm SHA-256 của binary (tệp thực thi) `'+lock['engine']['sha256']+'` trùng bản đóng băng trước lượt này. Không sửa matching (ghép cặp Morse), quotient, chọn PV, hỗ trợ witness hoặc cách dựng phức. Bộ chọn chỉ chuyển lịch sang CLI (giao diện dòng lệnh). Không dùng f2-stream/f2-rebuild; CONTINUE giữ trạng thái tạm dừng. Hồ sơ cũ được bảo toàn.',
      '',
      'Đã đọc đủ 7 trang hi.pdf (hash/text trong hi-proposal.*). Các tỷ lệ tăng tốc, chứng chỉ an toàn, tham số và mô tả Gcurr=|A|−1 trong đó chưa được xác minh trên code này. G v2 cũ đã dự đoán chênh lệch thời gian. Code PV hiện tại đã thu nhiều vùng rồi cập nhật chung; không nhận gom lô là cải tiến mới. Không áp dụng chính sách đổi PV trong tài liệu. Chi phí đặc trưng hiện tại O(n²×11), không phải O(1), và được tính trong D.',
      '',
      '3. **G chọn lịch nào?** Trong mọi danh sách, 2.0 là cutoff (ngưỡng kết thúc), không quotient ở đó. Mốc được tạo từ k nguyên chia 10. Có đúng 10 ứng viên, không quét tổ hợp khác.','',
      '| Bộ cuối | Lịch G | Trong miền | T̂ lịch (s) | G dự đoán (s) |','|---|---|---|---:|---:|']
    for x in choices:lines.append(f"| {x['dataset']} | {[k/10 for k in x['ticks']] if x['ticks'] else 'không quotient'} | {'Có' if x['in_domain'] else 'Không'} | {x['predictions_s'][x['schedule']]:.5f} | {x['G_seconds'][x['schedule']]:.5f} |")
    lines += ['',f"4. **So sánh tốc độ và RAM.** Lịch cố định (fixed) tốt nhất chọn trước final là `{lock['fixed_schedule']}` = {[k/10 for k in plan['schedules'][lock['fixed_schedule']]] if plan['schedules'][lock['fixed_schedule']] else 'không quotient'}. A: không quotient; B: tác giả [1.2,1.7,2]; C: lịch cố định đó; D: G gồm dự đoán; E: chạy lại đúng lịch G bỏ dự đoán. F: Ripser chỉ kiểm PH, không chọn lịch và không công bố benchmark (đo đối chứng hiệu năng) tốc độ/RAM với Ripser.",'',
      'Mỗi phương án có 1 warm-up (làm nóng) và 5 lượt đo, tuần tự, xoay và đảo thứ tự giữa lượt. Bảng dưới là median (trung vị) toàn tiến trình, giây; khoảng [min,max] là nhỏ nhất–lớn nhất. RSS (bộ nhớ thường trú) là peak (mức cao nhất) GNU time, KiB; không dùng peak theo cửa sổ verbose (log chi tiết) làm số đo toàn tiến trình. G dùng execv: Python nhường cùng tiến trình cho engine, không giữ thêm tiến trình Python chờ PH. Phép đo bắt đầu trước khởi động Python, gồm đọc, đặc trưng, dự đoán, ghi lựa chọn, PH, xuất và giải phóng. Diagnostic (lượt chẩn đoán) bật verbose riêng, không nằm trong bảng thời gian.',
      '', '| Bộ | Phương án | PH | Trung vị [min,max] s | RSS trung vị [min,max] KiB |','|---|---|---|---:|---:|']
    for r in perf:lines.append(f"| {r['dataset']} | {r['method']} | {r['PH_status']} | {r['median_s']:.5f} [{r['min_s']:.5f},{r['max_s']:.5f}] | {r['median_rss_kib']:.0f} [{r['min_rss_kib']},{r['max_rss_kib']}] |")
    lines += ['', 'Tỷ số thời gian X/D >1 nghĩa là G nhanh hơn X. Chỉ dùng giao các bộ mà cả X và D đều PASS; tổng hợp bằng trung bình nhân của tỷ số trung vị từng bộ. Không lấy ca sai PH làm thành tích phương pháp PH chính xác.']
    for m,x in comparisons.items():
        ratio=x['geometric_mean_median_runtime_ratio'];lines += ['',f"- {m}: {ratio:.4f} trên {x['datasets']}/{len(ids)} bộ; loại {len(x['excluded'])} bộ: {x['excluded']}." if ratio is not None else f"- {m}: chưa có cặp đủ PASS; loại: {x['excluded']}."]
    lines += ['', '5. **Đúng, sai, chưa xác minh.** PASS = khớp kiểm chứng; PH_MISMATCH = sai PH; UNVERIFIED = chưa xác minh; TIMEOUT/MEMORY_LIMIT = vượt thời gian/RAM. So multiset (tập giữ số lần lặp) của mọi khoảng; C++ xuất endpoint (đầu mút) 17 chữ số, so Decimal chính xác giữa các lịch. Ripser 0.6.15, F2, khoảng cách Euclid, maxdim=1, thresh=2, không chọn tập điểm đại diện, không xuất đối chu trình; đối chiếu sau chuyển endpoint iCC sang float32 (số thực 32 bit), không thêm tolerance (sai số chấp nhận). H0 lấy từ MSF (rừng khung nhỏ nhất) thật do clearing (loại cạnh đã giải H0) của iCC tạo ở từng cửa sổ; nối endpoint bằng ledger (sổ theo dõi), kiểm số thành phần liên thông và so chính xác với phương án tham chiếu không quotient và Ripser. H0 được kiểm trong lượt chẩn đoán riêng của đúng lịch; mọi lượt đo H1 đều so với lượt chẩn đoán đó. Đây là kiểm H0 có điều kiện khoảng cách dương, alpha<1; không suy ra bảo toàn H1 tổng quát.',
      '',
      'iCC lọc khoảng cách<epsilon, Ripser gồm equality. Survivor (khoảng còn sống) cuối có death=-1 ở iCC, +∞ ở Ripser; chỉ nói sống tới cutoff, không khẳng định sống mãi. Bỏ survivor tạm ở các cửa sổ chưa cuối. Không tự bỏ khoảng dài 0 ở parser: giữ mọi khoảng engine thực xuất; Ripser thường không xuất khoảng 0. Kernel (nhân tính toán) Morse bỏ cặp H1 có birth=death, parser (bộ đọc kết quả) không lọc thêm khoảng nào; số cặp dài 0 được kiểm riêng trong AUDIT.json. Đặc trưng kiểm điểm phân biệt, không có cạnh đúng cutoff 2 trong dữ liệu ngẫu nhiên này; các khác biệt cutoff/điểm trùng trong hồ sơ cũ vẫn tồn tại. Không nới sai số để đổi FAIL thành PASS.',
      '',f"Trên final, {len(correct)} cấu hình phương án/dataset (các alias vẫn được đo riêng); trạng thái từng phương án: {final_counts}. Tổng {len(rows)} lượt gồm làm nóng: {dict(Counter(x['status'] for x in rows))}. Khảo sát 180 cấu hình lịch/dataset: {summary['survey_counts']}. Bảng correctness.csv và survey-correctness.* giữ toàn bộ ca sai/chưa kiểm; final-correctness.json giữ chẩn đoán và H0.",
      '',f"Ca theo dõi 8 điểm tách khỏi tập huấn luyện/đánh giá cuối: {[(x['schedule'],x['status']) for x in read(run/'known-n8-sentinel.json')]}. Không phân tích lại nguyên nhân/thu nhỏ. Ca [0.5,2] nằm ngoài 10 ứng viên, dùng tái hiện lỗi đã biết. Các ca sai không làm dừng khảo sát các lịch khác, và không đổi lựa chọn final sau khi xem Ripser.",
      '', '6. **Chi phí controller có bù lợi ích?** Chi phí D−E dưới đây là chênh lệch hai trung vị toàn tiến trình, có nhiễu lịch hệ điều hành; paired (ghép theo vòng) min/max nằm trong controller-overhead.csv. pre_exec_s chỉ đo phần bên trong Python trước exec, thiếu khởi động interpreter và bước ghi cuối, nên không thay thế D−E.','',
      '| Bộ | D−E (s) | Phần trước exec (s) |','|---|---:|---:|']
    for x in overhead:lines.append(f"| {x['dataset']} | {x['difference_of_medians_s']:.6f} | {x['G_pre_exec_median_s']:.6f} |")
    c_ratio=comparisons['C/D']['geometric_mean_median_runtime_ratio']
    lines += ['',f"Trung bình chênh lệch trung vị D−E: {summary['mean_controller_difference_of_medians_s']:.6f}s. C/E (giá trị lịch trước phí G): {summary['fixed_over_replay_geomean']}; C/D (sau đủ phí G): {c_ratio}. "+('G thắng C trên giao PASS theo chỉ số tổng hợp này, nhưng mẫu nhỏ và một máy.' if c_ratio is not None and c_ratio>1 else 'G chưa thắng lịch cố định tốt trên giao PASS; không có bằng chứng chi phí chọn lịch được bù so với C.'),
      '',f"MAE dự đoán thời gian lịch đã chọn trên 6 bộ cuối, so với E: {summary['selected_schedule_prediction_MAE_s']:.6f}s. Đây là phân tích sau khóa, không dùng đổi lịch.",
      '',f"Học và khảo sát không miễn phí: từ thời điểm bắt đầu bảo thủ tới khóa mô hình mất {(datetime.fromisoformat(lock['locked_at'])-datetime.fromisoformat(budget['started_at'])).total_seconds():.1f}s, gồm đọc ngữ cảnh, triển khai, lỗi thiết lập, biên dịch, khảo sát và học. Các số phase (giai đoạn) dưới đây đo trong harness (chương trình điều phối) sau khi nạp thư viện Python, chưa gồm khởi động/nạp thư viện của chính harness. Mọi ca PH/Ripser bên trong vẫn có số đo toàn tiến trình riêng. Tổng thời gian các tiến trình đã ghi: {summary['total_logged_process_elapsed_s']:.3f}s (gồm cả build/kiểm chứng/đo, không cộng nhầm vào thời gian G).",
      '', '| Giai đoạn chuẩn bị | Thời gian harness (s) |','|---|---:|',
      *[f"| {s} | {x['elapsed_s']:.6f} |" for s,x in summary['offline_phase_costs'].items()],
      '', '7. **Điều đáng báo thầy và giới hạn.** Có một bộ chọn toàn lịch chạy được, mã nhỏ, dữ liệu độc lập trước/sau khóa, giữ nguyên thuật toán Morse/iCC và PV, đối chiếu H0–H1 riêng, đo đầy đủ phí điều khiển và giữ ca sai. Kết quả tốc độ và coverage (tỷ lệ bộ vượt kiểm PH) ở trên là kết quả thực, gồm cả bất lợi. Mô hình có thể chọn cùng một lịch cho nhiều hình học; khi đó chưa chứng minh lợi ích chọn lịch theo dữ liệu. Chưa có chứng chỉ an toàn quotient, chưa khắc phục upstream, chưa đánh giá H2 (đồng điều chiều 2), chưa chứng minh khái quát ngoài hai họ hình học/n64–200 và máy này. Không dùng thành tích phát triển thay cho final, không gọi đây là tái lập toàn bộ bài báo.',
      '',f"Tài nguyên: giới hạn process-group (nhóm tiến trình riêng của ca) {read(run/'initial-context.json')['memory_limit_bytes']/1024**3:.2f} GiB, timeout (hạn ca) tối đa 300s hoặc thời gian còn lại trước mốc đóng thí nghiệm. Monitor (bộ theo dõi) khoảng 50ms, không phải giới hạn cgroup cứng. Không đổi swap/WSL, không chạy đồng thời các ca nặng. Không dừng tiến trình người dùng.",
      '',f"Lỗi công cụ có {len(summary['tooling_errors'])} traceback lưu riêng, cộng initial-tooling-errors.json. Mọi raw command (lệnh gốc), stdout/stderr, exit code (mã kết thúc), mã băm dữ liệu đầu vào trước/sau và RSS ở cases/. Không sửa log cũ. INDEX.md, REPRODUCE.md và gói kèm manifest SHA-256 phục vụ bàn giao."]
    (run/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:summary[k] for k in ['final_counts','survey_counts','comparisons','fixed_schedule','choice_counts','mean_controller_difference_of_medians_s']},indent=2,ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();summarize(a.run.resolve())
