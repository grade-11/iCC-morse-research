#!/usr/bin/env python3
"""Summarize immutable evidence; never replace raw outcomes or endpoints."""
from workflow import *
from collections import Counter
import statistics,re

def main():
    def read(name):return json.loads((RUN/name).read_text())
    def probe_rss(val):
        context=json.loads((Path(val['path'])/'context.json').read_text())
        sampled=context.get('tree_rss_sampled_bytes',0)/1024
        gnu=int(context.get('resources',{}).get('max_rss_kib',-1))
        text=(Path(val['path'])/'stdout.log').read_text()
        internal=max((int(v) for v in re.findall(r'rss_peak_kernel_kib=(\d+)',text)),default=-1)
        return dict(probe_gnu_process_rss_kib=gnu,probe_sampled_group_rss_kib=sampled,probe_internal_window_max_kib=internal,probe_rss_kib=max(gnu,sampled,internal),probe_rss_source='max GNUtime/process-group sampled~50ms/internal window peaks; verbose can reset VmHWM')
    perf=list(csv.DictReader((RUN/'performance.csv').open()));vals=read('survey-validation.json');runs=read('survey-runs.json')
    for row in perf:
        val=next(v for v in vals if v['dataset']==row['dataset'] and str(v['dimension'])==row['dimension'] and v['method']==row['method'])
        c=read(str(Path(val['path']).relative_to(RUN)/'context.json'))
        rss=probe_rss(val)
        if val['status'].startswith('NOT_RUN'):rss={k:('Not run; blocking baseline context linked' if k=='probe_rss_source' else '') for k in rss}
        row.update(study='Native-exact baseline survey; follow-up n300 H2 float32 study is separate',warmups=sum(x['dataset']==row['dataset'] and str(x['dimension'])==row['dimension'] and x['method']==row['method'] and x['phase']=='warmup' for x in runs),probe_elapsed_s=val.get('elapsed_s',''),probe_exit_status=c.get('exit_status') if not val['status'].startswith('NOT_RUN') else '',**rss)
    table(RUN/'performance.csv',perf)
    h2perf=list(csv.DictReader((RUN/'n300-h2-performance.csv').open()))
    for row in h2perf:row.update(**probe_rss(next(x for x in read('n300-h2-validation.json') if x['method']==row['method'])))
    table(RUN/'n300-h2-performance.csv',h2perf)
    ratio=[]
    for dataset in sorted({x['dataset'] for x in runs}):
      for numerator,denominator,scope in [('native-defer-grid01','native-one','native checkpoint rebuild overhead'),('f2-stream-grid01','f2-rebuild-grid01','same explicit F2 kernel, retained state vs reset'),('f2-stream-grid01','native-one','different kernels: serial explicit F2 vs native4'),('native-author','native-one','fixed quotient, individually validated'),('native-prior-fixed','native-one','prior selected fixed quotient, individually validated')]:
        cells={m:{x['round']:x['elapsed_s'] for x in runs if x['dataset']==dataset and x['dimension']==2 and x['method']==m and x['phase']=='measured'} for m in (numerator,denominator)}
        shared=sorted(cells[numerator].keys()&cells[denominator].keys())
        if shared:
            ratios=[cells[numerator][k]/cells[denominator][k] for k in shared]
            ratio.append(dict(dataset=dataset,dimension=2,numerator=numerator,denominator=denominator,rounds=len(shared),median_time_ratio=statistics.median(ratios),min_time_ratio=min(ratios),max_time_ratio=max(ratios),scope=scope))
    table(RUN/'paired-ratios.csv',ratio)
    contexts=[]
    for p in (RUN/'cases').glob('*/context.json'):
        c=json.loads(p.read_text());assert 'exit_status' in c,'Incomplete child: '+str(p)
        contexts.append(c)
    stages=[]
    for p in RUN.glob('stage-*/context.json'):
        c=json.loads(p.read_text());assert 'exit_status' in c,'Incomplete outer stage: '+str(p)
        if c['exit_status']:stages.append(dict(path=str(p.parent),exit_status=c['exit_status'],error=c.get('error')))
    problems=[dict(case=c['case'],status=c['status'],exit_status=c['exit_status']) for c in contexts if c['status']!='OK']
    write(RUN/'TOOLING-ERRORS.json',dict(nonzero_or_resource_cases=problems,failed_stage_captures=stages,explanations={
      'remote_DNS':'Sandbox DNS failed exit128; direct read-only escalated ls-remote succeeded, raw both preserved.',
      'prototype_format':'Initial Python shortest-repr vs C++17digit decimal assertion was a harness comparison error; exact float64 comparison and .17g serialization corrected, old source/failure retained.',
      'CSV_export':'288 locked-action checks and52 G/replay checks completed before a duplicate-keyword CSV export TypeError. CSV recovered from saved JSON, no experiment rerun or result editing.',
      'duplicate_hypothesis':'Initial H1(max-simplex2) mismatch hypothesis was rejected. Actual zero-edge exclusion produces spurious H2 only in max-simplex3 no-quotient piecewise on duplicates4; scopes corrected explicitly.',
      'sanitizer':'LeakSanitizer explicitly unsupported under ptrace, initial CTest exit8. ASAN/UBSAN with detect_leaks=0 passed4/4; leak check remains NOT_VERIFIED.',
      'expected_rejections':'Three new-mode CLI unsupported-feature cases intentionally return nonzero; not successful PH runs.',
      'G_forecast_refusals':'12 predeclared old-grid low-budget stress cases return G_MEMORY_LIMIT exit1 (controller forecast refusal, not physical RSS kill). Initial harness asserted completion; corrected to NOT_VERIFIED incomplete PH, old source and failed stage retained; original cached failed native process not rerun.',
      'runner_guard':'Concurrent-stage smoke intentionally refused with exit1; no second experiment launched.',
      'threshold_metadata':'Only lossy window metadata reconstructed from exact argv binary64 (.17g); pair endpoints/raw logs unchanged. Audits preserve original340 multisets and check H0 at the decimal boundary.',
      'auxiliary_reads':'A read-only old-RSS helper compared a string to int (exit1), corrected via int cast; cat piped to head emitted SIGPIPE under an outer exit0. Neither changed experimental evidence.'}))
    inventory=dict(at=now(),original_matrix_configurations=len(read('correctness.json')),original_matrix_statuses=dict(Counter(x['status'] for x in read('correctness.json'))),cross_version=read('cross-version-summary.json'),counterexample8_quotient_configurations=len(read('counterexample8.json')),counterexample8_mismatches=sum(not x['quotient_equal_baseline'] for x in read('counterexample8.json')),fixed_instrumentation_equal=sum(x['PH_equal_with_locked_actions'] for x in read('instrumentation.json')),G_replay_equal=len(read('g-replay.json')),observed_G_trajectory_differences=len(read('g-trajectory-differences.json')),actual_H0_checks=len(read('h0-verification.json')),cpp_stream_configurations=len(read('stream-correctness.json')),cpp_stream_statuses=dict(Counter(x['status'] for x in read('stream-correctness.json'))),native_regression_unchanged=len(read('stream-native-regression.json')),survey_validations=dict(Counter(x['status'] for x in vals)),survey_warmups=sum(x['phase']=='warmup' for x in runs),survey_measured=sum(x['phase']=='measured' for x in runs),new_child_processes=len(contexts),child_process_states=dict(Counter(x['status'] for x in contexts)),new_process_elapsed_sum_s=sum(c.get('elapsed_monotonic_s',0) for c in contexts),reused_n200=read('n200-dossier.json'),n300_H2=read('n300-h2-summary.json'),G_grid01=read('G-GRID01-GATE.json'),budget_elapsed_s=time.monotonic()-budget()['start_monotonic'])
    write(RUN/'inventory.json',inventory)
    gb=read('G-boundary-summary.json');inventory['G_old_grid_boundary_diagnostic']=gb;write(RUN/'inventory.json',inventory)
    inventory.update(new_G_offline_training_s=0,new_G_phase='NOT_RUN_FOUNDATION_PH_MISMATCH');write(RUN/'inventory.json',inventory)
    expanded=read('grid01-expanded-upstream.json');verdicts=list(csv.DictReader((RUN/'case-verdicts.csv').open()))
    for x in expanded:
        if x['version']!='U17':continue
        case=x['dataset']+'-allQ-grid01'
        if any(r['case']==case for r in verdicts):continue
        verdicts.append(dict(case=case,dimension=2,threads=4,status=x['status'],verdict='Reproduced on upstream734, fixed allQ grid0.1; not G',H0='UPSTREAM_NOT_EXPORTED',upstream_U_path=next(y['quotient']['path'] for y in expanded if y['dataset']==x['dataset'] and y['version']=='U'),upstream_U17_path=x['quotient']['path']))
    upstreamG=read('G-selected-upstream.json')
    grows=read('G-boundary.json')+[{**x,'PH_reference_status':x['status']} for x in read('G-existing-reference.json')]
    for row in grows:
        if row['PH_reference_status']!='PH_MISMATCH':continue
        case='G-oldgrid-'+row['case']+'-mask'+''.join(map(str,row['mask']))
        if any(x['case']==case for x in verdicts):continue
        allq=all(row['mask']) and bool(upstreamG) and row['case']=='standard100'
        verdict='Reproduced on upstream734 actual allQ schedule' if allq else 'Not equivalent: upstream lacks mixed-action mask; group-vs-upstream cause not established for this mask' if not all(row['mask']) else 'Insufficient evidence: this actual G allQ schedule not yet tested on upstream'
        verdicts.append(dict(case=case,dimension=3,threads=row.get('threads',4),status='PH_MISMATCH',verdict=verdict,H0='UPSTREAM_NOT_EXPORTED; D replay actual ledger PASS',upstream_U_path=next((x['result']['path'] for x in upstreamG if x['version']=='U'),'') if allq else '',upstream_U17_path=next((x['result']['path'] for x in upstreamG if x['version']=='U17'),'') if allq else ''))
    table(RUN/'case-verdicts.csv',verdicts)
    def value(dataset,method,key='median_s'):
        r=next(x for x in perf if x['dataset']==dataset and x['method']==method and x['dimension']=='2')
        return f'{float(r[key]):.3f}' if r[key] not in ('',None) else r['status']
    rows='\n'.join('   | '+str(n)+' | '+' | '.join(value('standard'+str(n),m) for m in ('native-one','native-defer-grid01','f2-stream-grid01','f2-rebuild-grid01'))+' |' for n in (100,200,300))
    h2=read('n300-h2-summary.json');h2lines='; '.join(x['method']+': '+x['status']+(f", median {x['median_s']:.3f}s, {x['rounds']} lượt" if x['rounds'] else ', không đo lặp') for x in h2['performance'])
    gtext=f"Đối chiếu độc lập52 G cũ: 48 PASS,4 PH_MISMATCH; bốn allQ n100 sai đã tái hiện U/U17. Replay bằng nhau không đồng nghĩa PH đúng. Stress16 G với forecast budget12/16/32/48MiB:12 forecast-refusal (PH chưa xác minh),4 hoàn tất nhưng PH_MISMATCH; không quan sát khác trajectory do verbose. Đây là grid/model cũ, không phải G0.1; xem [chi tiết](G-boundary-summary.json), [kiểm PH](G-existing-reference.csv)."
    boot=budget();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO).decode().strip()
    report=f'''# Báo cáo truy nguyên PH và CONTINUE — 16/10/2026

**Lỗi H1 tái hiện trên upstream734; đã có engine giữ state H0–H1 riêng;
G trên lưới0.1 dừng theo gate tính đúng.**

Repo `/home/alex_1709/projects/iCC-team`; HEAD đầu lượt `259d4bf`, feature commit
`50636a5`; HEAD hồ sơ `{head[:12]}`. Bắt đầu {boot['started_at']}, hạn
{boot['deadline_at']}; ngân sách không được đặt lại. Không push, không reset/clean.
Nguồn, raw logs, lệnh, exit code và hash nằm trong gói; đọc [INDEX](INDEX.md).
Các so sánh chính dùng F2, Euclidean, strict cutoff2, alpha0.85;
cột dimension trong CSV là chiều simplex tối đa.

1. **Code gốc nào?** Đã xác minh trực tiếp remote [wilseypa/iCC commit734be3d](https://github.com/wilseypa/iCC/commit/734be3dd396f71217b7693b6acb382b81d7db6a9)
   là HEAD/main/parallel tại lúc kiểm tra. U là source nguyên trạng; U17 chỉ vá
   precision output6→17. Nhóm bắt đầu từ nền này; tới `3c59f45`, C++ không đổi.
   B17=`36b9db8` thêm precision/steady clock; C=`c8cd0df` có instrumentation;
   D=`259d4bf` là Gv2. Hai PDF chỉ dẫn URL, **chưa chốt được commit artifact
   của bài2025 hoặc2026**. Snapshot2025 `044f443` có entrypoint hardcode và Q&E
   khác, ghi NOT_EQUIVALENT, không tự sửa nó thành upstream nhiều cửa sổ.
   Xem [nguồn gốc](provenance.csv), [diff/history](provenance-git-commands.json).

2. **Lỗi có ở upstream?** Có trên các cấu hình phản ví dụ đã kiểm tra, gồm
   10/11 điểm cũ và ca mới v2. Ma trận gốc340 cấu hình, đối chiếu chéo392 lần
   đều khớp giữa các phiên bản ở precision tương ứng. U17/B17/C/D khớp multiset
   chính xác; U nguyên trạng khớp projection6 chữ số, không gọi là full64-bit PASS.
   Không tìm thấy cặp group PASS→FAIL nên không có cơ sở bisect lỗi của nhóm.
   Phản ví dụ mới **8 điểm phân biệt trên S² đơn vị**, alpha0.85, max-simplex2,
   lịch `[0.5,2]`: hai PV được implementation chấp nhận có đường kính0.40386759<0.425; còn hai nhãn nên
   mất H1 tương lai `(0.66666666666666652,0.71802197428460035)`. U và mọi bản nhóm
   tái hiện ở1/4 luồng; oracle kiểm đúng từng simplex/weight của traceS. AllQ grid0.1 cũng
   sai trên n8 và n200/n300 chuẩn seed20261007, đã chạy U/U17/D đối chứng trực tiếp.
   Lịch tác giả vẫn PASS trên những ô ghi PASS. Xem [verdict từng ca](case-verdicts.csv),
   [ma trận](correctness.csv), [n8](counterexample8-summary.json), [grid mở rộng](grid01-expanded-upstream.json).
   [Hình n8](counterexample8.png) mô tả hai bridge bị gộp thành một cạnh.
   Đây không phải chứng minh mọi ca lỗi cùng nguyên nhân hoặc mọi mệnh đề bài báo sai.

3. **Instrumentation ảnh hưởng gì?** 288 phép thử khóa hành động (timer,
   internal RAM, verbose, telemetry/trace;1/4 luồng, lặp2) giữ cùng PH. 52 G
   trajectory/replay trên grid cũ khớp PH và hash nhãn/PV/support; chưa quan sát trajectory
   đổi khi bật verbose/lặp trong mẫu này. Chưa loại trừ tác động timing ở mẫu
   khác: hệ số common-scale triệt trong margin tương đối, nhưng ngưỡng tuyệt đối
   0.5ms vẫn có thể đổi quyết định. Verbose còn reset VmHWM từng window, có thể
   đổi quan sát RAM. Clock-off/RAM-off chỉ dùng đối chứng fixed,
   không benchmark/huấn luyện G. H0 actual-MSF ledger có210 kiểm tra PASS trong
   phạm vi khoảng cách dương, alpha<1; upstream không xuất H0 nên ghi NOT_EXPORTED.
   Xem [instrumentation](instrumentation.csv), [G/replay](g-replay.csv), [H0](h0-verification.csv).
   {gtext}

4. **CONTINUE tới đâu?** Đã hiện thực mode C++ `f2-stream` giữ reduced pivots,
   live births và chỉ xử lý edges/triangles mới; `f2-rebuild` là control cùng
   kernel reset mỗi checkpoint. Đây là **explicit F2 kernel tuần tự, H0–H1**,
   chưa phải continuation native Morse. Quotient/G/v1/H2 bị từ chối rõ; chưa
   có migration sau quotient. Default native không đổi trên24 hồi quy.
   {inventory['cpp_stream_configurations']} cấu hình C++ kiểm ties, duplicates, strict cutoff,
   nhiều partition/k/10, lớp sống qua checkpoint rồi chết; state hash từng prefix
   khớp fresh engine, H1 khớp native, H0 khớp actual MSF. Các khác biệt quy ước
   zero/cutoff được ghi riêng trong [bảng kiểm](stream-correctness.csv).
   Release CTest4/4; ASAN/UBSAN4/4 khi tắt LSan không chạy được dưới ptrace;
   leak check chưa xác minh. Đọc [thiết kế](source/experiments/pwph2026-provenance-continue/CONTINUE-DESIGN.md),
   [phạm vi sanitizer/lỗi harness](TOOLING-ERRORS.json).

5. **Lưới0.1 có lợi?** Không chạy G mới vì nền quotient đã sai H1; không reuse
   model cũ cho miền mới, không fit hay mở final mới. Chỉ khảo sát checkpoint
   và fixed schedule có kiểm PH từng ô. Bảng sau là median giây H0–H1, một warm-up
   +5 lượt; toàn tiến trình gồm input, sort, checkpoint, output và teardown.
   RAM các lượt đo lấy GNUtime với verbose tắt, không reset VmHWM; RAM probe
   verbose ghi riêng các nguồn GNUtime, sampled process-group và kernel-window.

   | n | Native một lần | Native20 DEFER | F2 giữ state20 | F2 dựng lại20 |
   |---:|---:|---:|---:|---:|
{rows}

   Grid dày làm tăng phí rebuild native; giữ state giảm phí so cùng F2 reset,
   nhưng F2 hiện chậm hơn native về thời gian tuyệt đối. Không gán khác biệt
   giữa hai kernel chỉ cho state reuse. Native cấu hình4 luồng, F2 reduction serial.
   [Thời gian/RAM/min–max và ô thất bại](performance.csv), [tỷ số từng lượt](paired-ratios.csv),
   [validation](survey-validation.json), [gate G](G-GRID01-GATE.json). Ripser dùng
   kiểm PH float32/F2, không báo so tốc độ với Ripser; C++17digit so exact multiset,
   không nới tolerance. Epsilon tạo bằng k/10, không cộng dồn0.1.

6. **n200/n300?** N200 cũ seed20261033 tái sử dụng8 cấu hình PASS H0–H1/H0–H2,
   8 warm-up+40 lượt,508 file đã kiểm manifest cũ. Binary D dựng lại khớp hash
   `a950a64e…`; không gán timing này cho S. H0–H2 median cũ: lịch tác giả1.870s,
   fixed2.938s, G2.853s — G không thắng lịch tác giả trên ô này.
   [Hồ sơ n200](n200-dossier.json), [bảng cũ chọn lọc](n200-reused/performance-selected.csv).
   N200/n300 chuẩn seed20261007 có khảo sát H0–H1 mới trong bảng trên.
   N300 baseline H0–H2: **{h2['baseline_d3_status']}**; không gọi partial output là PH_PASS.
   Probe H2 quotient: {h2lines}. Khi thiếu baseline native H2 hoàn chỉnh,
   H2 chỉ đối chiếu Ripser float32, H0/H1 vẫn có exact native; ghi rõ phạm vi,
   không gọi là exact-H2 CPP PASS. [Chi tiết n300 H2](n300-h2-validation.json),
   [thời gian/RAM](n300-h2-performance.csv). Giới hạn mỗi process300s/~6.05GiB,
   monitor process-group~50ms, không phải hard cgroup; không đổi swap/WSL.
   N300 là bước trung gian, không phải cấu hình S² TableI n100/n200/n400.

7. **Có thể báo thầy:** lỗi PH cụ thể tái hiện trên upstream734; cap/PV được implementation chấp nhận
   không đủ bảo toàn H1 tương lai; engine giữ state H0–H1 riêng đã có đối chứng;
   grid dày/rebuild có phí đo được và có ô allQ sai PH. **Chưa được kết luận:**
   artifact bài báo đã được chốt, quotient/G đã sửa đúng tổng quát, native Morse
   CONTINUE/H2 hoặc migration đã hoàn thành, G0.1 có lợi, hay đã tái lập đầy đủ
   bài báo. Các lỗi, giới hạn RAM, giả thuyết bị bác và kết quả bất lợi đều giữ lại.
'''
    (RUN/'REPORT.md').write_text(report)
    (RUN/'INDEX.md').write_text('''# Thứ tự đọc gói bàn giao

1. [REPORT.md](REPORT.md): bảy câu trả lời và phạm vi kết luận.
2. [provenance.csv](provenance.csv), [case-verdicts.csv](case-verdicts.csv), [correctness.csv](correctness.csv), [cross-version.csv](cross-version.csv): nguồn gốc và đối chứng.
3. [counterexample8-summary.json](counterexample8-summary.json), [counterexample8-plan.json](counterexample8-plan.json), [grid01-expanded-upstream.json](grid01-expanded-upstream.json): lỗi H1 mới, construction/config/hash và kết quả upstream.
4. [CONTINUE-DESIGN](source/experiments/pwph2026-provenance-continue/CONTINUE-DESIGN.md), [stream-correctness.csv](stream-correctness.csv), [README/lệnh tái hiện](source/experiments/pwph2026-provenance-continue/README.md): engine H0/H1, giới hạn và build offline.
5. [performance.csv](performance.csv), [paired-ratios.csv](paired-ratios.csv), [survey-validation.json](survey-validation.json): chi phí checkpoint, từng ô PASS/FAIL/resource; không có G mới.
6. [n200-dossier.json](n200-dossier.json), [n200-reused/copy-map.json](n200-reused/copy-map.json), [n300-h2-summary.json](n300-h2-summary.json): bằng chứng n200 tái dùng và giới hạn n300.
7. [inventory.json](inventory.json), [TOOLING-ERRORS.json](TOOLING-ERRORS.json), [PROGRESS.md](PROGRESS.md), [budget.json](budget.json): số lượng, lỗi/phục hồi, một tiến độ và ngân sách ban đầu.
8. cases/: command.txt, stdout/stderr, context/exit-status/resources, PH checks. inputs/: nguyên bytes/hash; dataset-inventory.csv ghi hình học. references/: config+NPZ+SHA; reused-reference-map.json chỉ tới raw tham chiếu cũ đã copy chọn lọc. versions/: U/U17/B17/C/D, hai diagnostic variants và historical2025; source/: snapshot S hiện tại; third-party/CLI11/: source dependency. build-context/: flags/cache/test logs, không executable. offline-smoke.json và evidence-audit.json ghi kiểm build/integrity.
9. changes.patch, git-final.json, MANIFEST.sha256.json: thay đổi và kiểm integrity. Archive nằm cạnh thư mục runs/, có file .sha256 bên cạnh; packaging-capture/ lưu kết quả kiểm gói bên ngoài archive để tránh tự thay đổi manifest.
''')
    progress('research-summary','COMPLETE_WITH_SCOPE','Seven questions answered; new G gated, no native quotient fix/migration or full paper artifact claim')
    print(json.dumps(inventory,indent=2))

if __name__=='__main__':main()
