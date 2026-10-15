#!/usr/bin/env python3
"""A short Vietnamese handoff, derived from saved evidence."""
from summarize import *

def report(run):
    inv=json.loads((run/'inventory.json').read_text());lock=json.loads((run/'LOCK-v2.json').read_text())
    checks=json.loads((run/'correctness.json').read_text());final=[x for x in checks if x['split']=='final']
    perf=list(csv.DictReader((run/'performance.csv').open()));ratios=list(csv.DictReader((run/'paired-ratios.csv').open()))
    primary=[x for x in final if x['method']=='G' and x['maxdim']==2];ext=[x for x in final if x['method']=='G' and x['maxdim']==3]
    def outcomes(rs):return ', '.join(f'{k} {v}/{len(rs)}' for k,v in sorted(Counter(x['status'] for x in rs).items()))
    def get(name,d,m):return next((x for x in perf if x['dataset']==name and int(x['dimension'])==d and x['method']==m and x['status']=='PASS'),None)
    def val(x,k='median_s',scale=1):return '—' if x is None else f"{float(x[k])/scale:.3f}"
    good=[x for x in ratios if int(x['dimension'])==2 and x['reference']=='best_fixed']
    faster=sum(float(x['median_time_ratio'])>1 for x in good)
    aggregate=list(csv.DictReader((run/'aggregate-ratios.csv').open()))
    def gm(d,reference):return float(next(x['geometric_mean_time_ratio'] for x in aggregate if int(x['dimension'])==d and x['reference']==reference))
    rp=run/'source-reproduction.json';reproduction=json.loads(rp.read_text()) if rp.exists() else {}
    reproduced=reproduction.get('status')=='PASS' and reproduction.get('binary_sha256')==reproduction.get('benchmark_binary_unchanged')
    reproline='Source+CLI11 đóng gói build lại offline đạt 3/3 CTest, binary SHA256 trùng binary benchmark. ' if reproduced else 'Chưa xác nhận build offline từ gói. '
    text=[
      '# Báo cáo cho Triết — PwPH v2, 16/10/2026','',
      f"**Chưa sửa được mất thông tin PH nền. Đã xây G điều khiển quotient/defer thật. Final G: max-dim 2 (H0–H1) {outcomes(primary)}; max-dim 3 (H0–H2) {outcomes(ext)}.** "
      f"Chưa có lợi ích tốc độ so fixed tốt: fixed/G={gm(2,'best_fixed'):.3f} ở H0–H1 và {gm(3,'best_fixed'):.3f} ở H0–H2 (<1 bất lợi cho G). Không suy ra bảo toàn PH tổng quát.",'',
      f"Có {inv['unique_new_datasets']} dataset mới ({inv['datasets_by_split']}), n32/64/100/200, S2 đều, S2 tập trung cực và cụm đa thang đo. "
      'Seed07–11 chỉ là regression; seed final30–37 chỉ sinh sau lock. Không chỉnh thiết kế sau final. '
      f"Đến tổng hợp đã dùng {inv['budget_elapsed_s']/60:.1f} phút/240 phút (monotonic), gồm truy nguyên, build, calibration và đo; deadline cũ giữ nguyên khi resume.",'',
      '## Độ đúng đã kiểm, lỗi còn lại','',
      'Hai phản ví dụ cũ tái lập ở 1/4 luồng. Oracle độc lập khớp từng simplex/weight witness và reduction F2 của C++; '
      'nhưng phức sau quotient không còn đủ thông tin VR: ca 11 điểm còn ba nhãn trước H2 mới; ca 10 điểm còn hai nhãn trước hai H1 và một H2 mới. '
      'PV vẫn qua các bộ lọc đang cài đặt. Không có căn cứ vá parser/matching hoặc đổi alpha để ép PASS. Sáu lỗi n100 đầy đủ vẫn giữ nguyên; chưa chứng minh cùng nguyên nhân với ca nhỏ. '
      'Max-dim 2 là cấu hình kiểm riêng, không che lỗi H1 đã biết ở max-dim 3.','',
      'Bổ sung riêng bộ kiểm H0 từ MSF thật từng cửa sổ, kiểm conservation và toàn multiset với baseline/Ripser. '
      'Giả thiết: khoảng cách dương, alpha∈(0,1), luật disjoint/frozen PV hiện tại. Sáu unit tests và audit 576 cấu hình đạt; 18 H0 chưa xác minh cũ được bổ sung, evidence gốc giữ nguyên. '
      'Đây là kiểm chứng H0 có điều kiện, không sửa H1/H2. CTest, n100 chuẩn, regression cũ và nhánh toàn DEFER được kiểm. '+reproline+
      'Endpoints/multiplicities giữa iCC so đầy đủ; Ripser theo float32 đã chốt, không nới tolerance. Fixture cutoff/điểm trùng vẫn KNOWN_DIFFERENCE.','',
      '## Đo final: thời gian toàn tiến trình và RSS','',
      'Mỗi nhóm đủ PH có 1 warmup +5 vòng xoay/đảo thứ tự, tuần tự 4 luồng, cùng binary, không verbose. '
      'Bảng max-dim 2: giây trung vị; tỷ số paired fixed/G >1 lợi cho G; RSS MiB. — là không đủ điều kiện, không phải 0.','',
      '| Final (hình học/n/seed) | Fixed / G (s) | Fixed/G [min,max] | RSS fixed / G | G PH |',
      '|---|---:|---:|---:|---|']
    for item in lock['final_plan']:
        name=item['id'];f=get(name,2,'best_fixed');g=get(name,2,'G');p=next((x for x in good if x['dataset']==name),None)
        rr='—' if p is None else f"{float(p['median_time_ratio']):.2f} [{float(p['min_time_ratio']):.2f},{float(p['max_time_ratio']):.2f}]"
        status=next(x['status'] for x in primary if x['dataset_id']==name)
        text.append(f"| {item['geometry']}/{item['n']}/{str(item['seed'])[-2:]} | {val(f)} / {val(g)} | {rr} | {val(f,'median_rss_kib',1024)} / {val(g,'median_rss_kib',1024)} | {status} |")
    text+=['',f"Đã khóa theo chiều và nhóm n≤64/n≥100: `{lock['selected']}`. "
      'Chọn bằng coverage PH, rồi H0, rồi thời gian development/validation. Bảng đầy đủ gồm fixed gốc, fixed tốt, v1, G, H2 và min–max ở performance.csv/paired-ratios.csv. '
      f"Các tỷ số tổng hợp là geometric mean của median paired ratio từng dataset: v1/G={gm(2,'v1'):.3f} (D2), {gm(3,'v1'):.3f} (D3); fixed gốc/G={gm(2,'original'):.3f}/{gm(3,'original'):.3f}. "
      'G chọn cùng hành động với fixed tốt ở 7/8 ca mỗi chiều. Final chỉ có một checkpoint nên không đo lợi ích EMA cho quyết định sau. '
      'Cùng lịch/hành động mà lệch nhỏ chưa chứng minh giá trị scheduling; tỷ số có điều kiện trên ca đủ đúng, không bao gồm ca sai.','',
      '## G và chất lượng dự báo','',
      '`G = Ccontinue − (Cselect+quotient + Crebuild + Cfuture_afterQ)` từ cùng trạng thái tới cùng h, đơn vị giây. '
      'h là checkpoint kế tiếp; G chọn có quotient hay trì hoãn, không chọn PV theo chi phí. Ridge đơn giản theo nhãn/cofacets/PV/support/độ rộng, cập nhật scale từ quá khứ. '
      'Quotient khi G>max(0.5ms; margin×Ccontinue); RAM KiB là ràng buộc riêng với 30% headroom. Ngoài miền dự báo fallback DEFER; không gọi là an toàn PH tổng quát. '
      'Toàn telemetry/controller/checkpoint/logging tính vào runtime.','']
    example=None
    for x in primary+ext:
        log=(Path(x['path'])/'stdout.log').read_text();gs=[dict(t.split('=') for t in l.split()[1:]) for l in log.splitlines() if l.startswith('[g]')];ss=parse_states(log)
        if gs and len(ss)>1:example=x,gs[0],ss[0],ss[1];break
    if example:
        x,p,s,n=example
        text += [f"Ví dụ `{x['dataset_id']}` D{x['maxdim']}, ε={p['eps']}→h={p['h']}, A={s['A']}, C={s['C']}: "
          f"continue dự báo {float(p['continue_s'])*1000:.2f}ms; Q/rebuild/future={float(p['select_quotient_s'])*1000:.2f}/{float(p['rebuild_s'])*1000:.2f}/{float(p['future_s'])*1000:.2f}ms; "
          f"G={float(p['G_s'])*1000:.2f}ms, ngưỡng {float(p['threshold_s'])*1000:.2f}ms → **{p['action']}** ({p['reason']}). "
          f"Work cửa sổ sau quan sát {float(n['work_s'])*1000:.2f}ms; final chỉ đo nhánh đã chọn.",'']
    models=json.loads((run/'G-prediction-summary.json').read_text())
    for source,label in [('estimator-validation.json','Ridge .01 ban đầu'),('estimator-global-validation.json','Bổ sung chân trời cuối')]:
        rows=[x for x in models if x['source']==source and x['ridge']==.01]
        text.append(label+': '+ '; '.join(f"D{x['dimension']} {x['samples']} cặp, MAE {x['G_MAE_s']*1000:.1f}ms, đúng dấu {x['sign_accuracy']:.1%} (luôn dương {x['always_positive_accuracy']:.1%})" for x in rows)+'.')
    text+=['Bỏ nhóm PV/support làm MAE validation xấu hơn nên giữ feature đầy đủ; đúng dấu riêng có thể tốt hơn, không giấu kết quả đó. '
      'Nhãn G_actual là hiệu chi phí các phase không chồng lấn, chưa phải chi phí phản thực whole-process chính xác: phí hủy support-container Q và phần controller/log chưa tách. '
      'Biên 0,5ms chưa là cận nhiễu được chứng minh. Chi phí replay/fitting offline ghi riêng; online không chạy thử hai nhánh.','',
      '## Kiểm kê và bàn giao','']
    shr=json.loads((run/'new-shrink-summary.json').read_text())
    text.append('Thu nhỏ lỗi mới: '+ ('; '.join(f"{x['case']} D{x['dimension']}: {x['original_n']}→{x['retained_n']} điểm" for x in shr) if shr else 'không có ca phù hợp')+'. Input/trace giữ nguyên, không tuyên bố tối tiểu.')
    ram=json.loads((run/'ram-study.json').read_text())
    text.append(f"RAM study: {ram['status']}"+(f", đổi hành động={ram['action_changed']}, trạng thái tương đương={ram['equivalent_initial_state']}; đây là phản ứng dự báo, chưa chứng minh tiết kiệm RAM thực." if 'action_changed' in ram else '; không chạy hai nhãn ngân sách vô nghĩa.'))
    stability=json.loads((run/'stability.json').read_text())
    text.append('Hai cặp validation khó được lặp để đo nhiễu: '+ '; '.join(f"D{x['dimension']} G thực min–max {x['min_G_actual_s']*1000:.1f}–{x['max_G_actual_s']*1000:.1f}ms" for x in stability)+'.')
    text += [f"{inv['validation_configurations']} cấu hình logic ({inv['unique_executed_validation_configurations']} thực thi riêng), trạng thái {inv['validation_statuses']}; "
      f"{inv['warmups']} warmup +{inv['measured']} đo final, trạng thái {inv['measurement_statuses']}. "
      'Alias không tính thêm process/dataset. TIMEOUT/MEMORY_LIMIT/ERROR giữ riêng trong inventory/raw logs. '
      'Lỗi tooling snapshot/replay/orchestration và lần auto-review commit timeout được lưu; retry không xóa evidence.','',
      'Chưa có điều kiện đủ bảo toàn PH mọi quotient; estimator còn sai số, mẫu final nhỏ và chỉ một máy. Chưa làm giai đoạn 3 hoặc hoàn tất giai đoạn 2. '
      'Đọc PH-ROOT-CAUSE.md, H0-LEDGER.md, CONTROLLER-v2.md; inventory/correctness/performance/G-predictions/AUDIT trong run. '
      'README có resume giữ deadline. Một archive kèm manifest/SHA256 chứa source, patch, tests, input, model và raw logs, không build/venv/cache. Commit cục bộ, không push.']
    content='\n'.join(text)+'\n';(BASE/'BAO-CAO-TRIET.md').write_text(content);(run/'BAO-CAO-TRIET.md').write_text(content)
    print(json.dumps(dict(report=str(BASE/'BAO-CAO-TRIET.md'),words=len(content.split()),primary=outcomes(primary),extended=outcomes(ext)),ensure_ascii=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();report(a.run)
