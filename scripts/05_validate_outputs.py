#!/usr/bin/env python3
from pathlib import Path
import csv, json, math, hashlib, platform, sys
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

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def compare_tsv(generated, reference, keys, numeric_cols):
    g=read_tsv(generated); r=read_tsv(reference)
    result={"generated":str(generated.relative_to(ROOT)),"reference":str(reference.relative_to(ROOT)),"rows_generated":len(g),"rows_reference":len(r),"passed":True,"issues":[]}
    if len(g)!=len(r): result["passed"]=False; result["issues"].append("row_count_mismatch")
    gm={tuple(x[k] for k in keys):x for x in g}; rm={tuple(x[k] for k in keys):x for x in r}
    if set(gm)!=set(rm): result["passed"]=False; result["issues"].append("key_set_mismatch")
    maxdiff=0.0
    for key in set(gm)&set(rm):
        for col in numeric_cols:
            try: diff=abs(float(gm[key][col])-float(rm[key][col])); maxdiff=max(maxdiff,diff)
            except (ValueError,TypeError,KeyError):
                result["passed"]=False; result["issues"].append(f"non_numeric_or_missing:{key}:{col}"); continue
            if diff>TOL:
                result["passed"]=False; result["issues"].append(f"numeric_difference:{key}:{col}:{diff}")
                if len(result["issues"])>20: break
        if len(result["issues"])>20: break
    result["maximum_numeric_difference"]=maxdiff
    return result

def main():
    checks=[]
    norm_g=GEN/'normalized/malat1_a549_normalized_decision_matrix_v0_3.tsv'
    norm_r=REF/'normalized/malat1_a549_normalized_decision_matrix_v0_3.tsv'
    checks.append(compare_tsv(norm_g,norm_r,["candidate_id"],["C1_position_soft_support","C6_Jost_CRISPRi_specificity","balanced_score","efficacy_focused_score","safety_focused_score"]))
    checks.append(compare_tsv(GEN/'profile_smaa/malat1_a549_smaa_profile_summary_v0_1.tsv',REF/'profile_smaa/malat1_a549_smaa_profile_summary_v0_1.tsv',["profile","candidate_id"],["p_rank_1","p_top_5","mean_rank","mean_utility","strict_no_1mm_p_top_5"]))
    checks.append(compare_tsv(GEN/'profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv',REF/'profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv',["candidate_id"],["mean_p_top_5_across_profiles","minimum_p_top_5_across_profiles","mean_rank_across_profiles"]))
    checks.append(compare_tsv(GEN/'extended/tables/global_weight_smaa_summary.tsv',REF/'extended/global_weight_smaa_summary.tsv',["candidate_id"],["p_rank_1","p_top_5","mean_rank","strict_no_1mm_p_top_5"]))
    checks.append(compare_tsv(GEN/'extended/tables/integrated_candidate_technical_status.tsv',REF/'extended/integrated_candidate_technical_status.tsv',["candidate_id"],["integrated_order","evidence_flag_count_0_to_7","three_profile_mean_p_top5","global_weight_p_top5"]))
    raw=read_tsv(ROOT/'inputs/malat1_a549_raw_decision_matrix_v0_2.tsv')
    normalized=read_tsv(norm_g)
    overall=read_tsv(GEN/'profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv')
    integrated=read_tsv(GEN/'extended/tables/integrated_candidate_technical_status.tsv')
    invariants={
        "raw_candidate_count_86":len(raw)==86,
        "normalized_candidate_count_86":len(normalized)==86,
        "unique_candidate_ids":len({r['candidate_id'] for r in normalized})==86,
        "all_jost_present":all(r['jost_crispri_specificity'] not in ('','NA','nan') for r in raw),
        "one_mismatch_candidate_count_6":sum(int(r['one_mismatch_hits'])>0 for r in raw)==6,
        "no_extra_exact_offtarget":all(int(r['exact_offtargets_excluding_intended'])==0 for r in raw),
        "core_integrated_top3":[r['candidate_id'] for r in sorted(integrated,key=lambda x:int(x['integrated_order']))[:3]]==['MALAT1_A549_NGG_0046','MALAT1_A549_NGG_0047','MALAT1_A549_NGG_0048'],
        "profile_robustness_top1":sorted(overall,key=lambda x:int(x['robustness_order']))[0]['candidate_id']=='MALAT1_A549_NGG_0046',
    }
    passed=all(c['passed'] for c in checks) and all(invariants.values())
    report={
        "passed":passed,
        "numeric_tolerance":TOL,
        "python":sys.version,
        "platform":platform.platform(),
        "numpy":np.__version__,
        "comparisons":checks,
        "invariants":invariants,
        "generated_confidence_factor_file_exists":(GEN/'profile_smaa/malat1_a549_smaa_confidence_factors_v0_1.tsv').exists(),
        "interpretation":"PASS means the post-FlashFry scientific tables reproduce the locked reference within tolerance; PNG/XLSX binary hashes are intentionally not required.",
    }
    REPORT_JSON.write_text(json.dumps(report,indent=2),encoding='utf-8')
    lines=["# Reproducibility validation report","",f"**Overall: {'PASS' if passed else 'FAIL'}**","",f"Tolerance: `{TOL}`","","## Invariants"]
    for k,v in invariants.items(): lines.append(f"- {'PASS' if v else 'FAIL'} — {k}")
    lines += ["","## Output comparisons"]
    for c in checks: lines.append(f"- {'PASS' if c['passed'] else 'FAIL'} — `{c['generated']}`; max difference `{c['maximum_numeric_difference']}`")
    lines += ["","PNG and XLSX binary hashes are not used as scientific equality tests because rendering and archive metadata can differ across systems."]
    REPORT_MD.write_text('\n'.join(lines),encoding='utf-8')
    print(f"Validation {'PASS' if passed else 'FAIL'} -> {REPORT_JSON}")
    raise SystemExit(0 if passed else 1)

if __name__=='__main__': main()
