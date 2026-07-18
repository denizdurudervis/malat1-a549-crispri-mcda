#!/usr/bin/env python3
from pathlib import Path
import csv, json, math, platform, sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "generated"
REF = ROOT / "reference"
REPORT_JSON = GEN / "reproducibility_validation_report.json"
REPORT_MD = GEN / "reproducibility_validation_report.md"
TOL = 5e-8

def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def numeric_value(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

def compare_tsv(generated, reference, keys):
    g=read_tsv(generated); r=read_tsv(reference)
    result={"generated":str(generated.relative_to(ROOT)),"reference":str(reference.relative_to(ROOT)),"rows_generated":len(g),"rows_reference":len(r),"passed":True,"issues":[]}
    if len(g)!=len(r): result["passed"]=False; result["issues"].append("row_count_mismatch")
    if (list(g[0]) if g else []) != (list(r[0]) if r else []):
        result["passed"]=False; result["issues"].append("column_mismatch")
    gm={tuple(x[k] for k in keys):x for x in g}; rm={tuple(x[k] for k in keys):x for x in r}
    if len(gm)!=len(g): result["passed"]=False; result["issues"].append("duplicate_generated_keys")
    if len(rm)!=len(r): result["passed"]=False; result["issues"].append("duplicate_reference_keys")
    if set(gm)!=set(rm): result["passed"]=False; result["issues"].append("key_set_mismatch")
    maxdiff=0.0
    for key in set(gm)&set(rm):
        for col in rm[key]:
            if col in keys:
                continue
            gv=gm[key].get(col); rv=rm[key].get(col)
            gn=numeric_value(gv); rn=numeric_value(rv)
            if gn is not None and rn is not None:
                if math.isnan(gn) and math.isnan(rn):
                    diff=0.0
                else:
                    diff=abs(gn-rn)
                maxdiff=max(maxdiff,diff)
                if diff>TOL:
                    result["passed"]=False; result["issues"].append(f"numeric_difference:{key}:{col}:{diff}")
            elif gv != rv:
                result["passed"]=False; result["issues"].append(f"text_difference:{key}:{col}")
            if len(result["issues"])>20: break
        if len(result["issues"])>20: break
    result["maximum_numeric_difference"]=maxdiff
    return result

def main():
    locked_outputs=[
        ('normalized/malat1_a549_normalized_decision_matrix_v0_3.tsv',["candidate_id"]),
        ('profile_smaa/malat1_a549_smaa_central_weights_v0_1.tsv',["profile","candidate_id"]),
        ('profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv',["candidate_id"]),
        ('profile_smaa/malat1_a549_smaa_profile_summary_v0_1.tsv',["profile","candidate_id"]),
        ('profile_smaa/malat1_a549_smaa_rank_acceptability_v0_1.tsv',["profile","candidate_id","rank"]),
        ('extended/global_weight_smaa_summary.tsv',["candidate_id"]),
        ('extended/integrated_candidate_technical_status.tsv',["candidate_id"]),
        ('extended/monte_carlo_convergence.tsv',["iterations"]),
        ('extended/offtarget_model_summary.tsv',["offtarget_model","profile"]),
        ('extended/pareto_front_and_layers.tsv',["candidate_id"]),
        ('extended/position_model_summary.tsv',["position_model","profile"]),
        ('extended/random_seed_summary.tsv',["seed_index"]),
        ('extended/sequence_qc_summary.tsv',["sequence_scenario","profile"]),
    ]
    checks=[]
    for rel,keys in locked_outputs:
        if rel.startswith('extended/'):
            generated=GEN/'extended/tables'/rel.removeprefix('extended/')
        else:
            generated=GEN/rel
        checks.append(compare_tsv(generated,REF/rel,keys))
    norm_g=GEN/'normalized/malat1_a549_normalized_decision_matrix_v0_3.tsv'
    raw=read_tsv(ROOT/'inputs/malat1_a549_raw_decision_matrix_v0_2.tsv')
    normalized=read_tsv(norm_g)
    overall=read_tsv(GEN/'profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv')
    integrated=read_tsv(GEN/'extended/tables/integrated_candidate_technical_status.tsv')
    expected_extended_figures=[
        'Figure_7_position_model_rank_ranges.png',
        'Figure_8_sequence_qc_rank_ranges.png',
        'Figure_9_global_weight_smaa_top10.png',
        'Figure_10_pareto_efficacy_local_safety.png',
        'Figure_11_mds_criterion_space.png',
        'Figure_12_uncertainty_source_rank_impact.png',
        'Figure_13_monte_carlo_convergence.png',
        'Figure_14_historical_guide_tss_distance.png',
    ]
    expected_main_figures=[
        'Figure_1_A549_TSS_evidence_and_candidate_window',
        'Figure_2_role_structured_panel_criteria',
        'Figure_3_panel_profile_and_global_robustness',
        'Figure_4_tradeoff_and_four_domain_pareto_projection',
        'Figure_5_uncertainty_source_rank_impact',
        'Figure_6_historical_guide_distance_full_and_zoom',
    ]
    expected_supplementary_figures=[
        'Figure_S1_balanced_rank_acceptability',
        'Figure_S2_global_weight_top10',
        'Figure_S3_position_model_rank_ranges',
        'Figure_S4_sequence_QC_rank_ranges',
        'Figure_S5_offtarget_model_rank_ranges',
        'Figure_S6_seed_stability_and_convergence',
        'Figure_S7_MDS_criterion_space',
        'Figure_S8_leave_one_domain_out',
    ]
    extended_dir=GEN/'extended/figures'
    main_dir=GEN/'manuscript_figures/main_figures'
    supplementary_dir=GEN/'manuscript_figures/supplementary_figures'
    pngs=(
        [extended_dir/name for name in expected_extended_figures]
        + [main_dir/f'{name}.png' for name in expected_main_figures]
        + [supplementary_dir/f'{name}.png' for name in expected_supplementary_figures]
    )
    pdfs=(
        [main_dir/f'{name}.pdf' for name in expected_main_figures]
        + [supplementary_dir/f'{name}.pdf' for name in expected_supplementary_figures]
    )
    pngs_valid=all(
        path.is_file()
        and path.stat().st_size>100
        and path.read_bytes()[:8]==b'\x89PNG\r\n\x1a\n'
        for path in pngs
    )
    pdfs_valid=all(
        path.is_file()
        and path.stat().st_size>100
        and path.read_bytes()[:4]==b'%PDF'
        for path in pdfs
    )
    invariants={
        "raw_candidate_count_86":len(raw)==86,
        "normalized_candidate_count_86":len(normalized)==86,
        "unique_candidate_ids":len({r['candidate_id'] for r in normalized})==86,
        "all_jost_present":all(r['jost_crispri_specificity'] not in ('','NA','nan') for r in raw),
        "one_mismatch_candidate_count_6":sum(int(r['one_mismatch_hits'])>0 for r in raw)==6,
        "no_extra_exact_offtarget":all(int(r['exact_offtargets_excluding_intended'])==0 for r in raw),
        "core_integrated_top3":[r['candidate_id'] for r in sorted(integrated,key=lambda x:int(x['integrated_order']))[:3]]==['MALAT1_A549_NGG_0046','MALAT1_A549_NGG_0047','MALAT1_A549_NGG_0048'],
        "profile_robustness_top1":sorted(overall,key=lambda x:int(x['robustness_order']))[0]['candidate_id']=='MALAT1_A549_NGG_0046',
        "all_22_expected_png_figures_valid":pngs_valid,
        "all_14_expected_pdf_figures_valid":pdfs_valid,
    }
    passed=all(c['passed'] for c in checks) and all(invariants.values())
    report={
        "passed":passed,
        "numeric_tolerance":TOL,
        "python":sys.version,
        "platform":platform.platform(),
        "numpy":np.__version__,
        "locked_output_count":len(checks),
        "comparisons":checks,
        "invariants":invariants,
        "generated_confidence_factor_file_exists":(GEN/'profile_smaa/malat1_a549_smaa_confidence_factors_v0_1.tsv').exists(),
        "interpretation":"PASS means every column in all locked post-FlashFry scientific TSV outputs reproduces the reference within tolerance and all expected PNG/PDF figure artifacts are valid; figure byte equality is intentionally not required.",
    }
    REPORT_JSON.write_text(json.dumps(report,indent=2),encoding='utf-8')
    lines=["# Reproducibility validation report","",f"**Overall: {'PASS' if passed else 'FAIL'}**","",f"Tolerance: `{TOL}`","","## Invariants"]
    for k,v in invariants.items(): lines.append(f"- {'PASS' if v else 'FAIL'} — {k}")
    lines += ["","## Output comparisons"]
    for c in checks: lines.append(f"- {'PASS' if c['passed'] else 'FAIL'} — `{c['generated']}`; max difference `{c['maximum_numeric_difference']}`")
    lines += ["",f"Validated locked scientific TSV outputs: `{len(checks)}`.","","All expected figure files were checked for valid nonempty PNG or PDF artifacts. Figure byte equality is not used as a scientific equality test because rendering and font metadata can differ across systems."]
    REPORT_MD.write_text('\n'.join(lines),encoding='utf-8')
    print(f"Validation {'PASS' if passed else 'FAIL'} -> {REPORT_JSON}")
    raise SystemExit(0 if passed else 1)

if __name__=='__main__': main()
