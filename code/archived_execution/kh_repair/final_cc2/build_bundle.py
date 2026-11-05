from pathlib import Path
import json,hashlib,difflib
D=Path(__file__).resolve().parent;OLD=D.parent;W=OLD.parent;ROOT=W/'refit_final_20261104';B=ROOT/'04_분석/BRCA'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
ready=json.loads((ROOT/'DATA_READY.json').read_text());panel=B/'panel_BRCA_std.csv';graph=ROOT/'input_graph.graph'
for rel in ['input_graph.graph','04_분석/BRCA/panel_BRCA_std.csv','04_분석/BRCA/config_BRCA_std.json','04_분석/BRCA/panel_META_BRCA_std.csv','04_분석/BRCA/panel_INFO_BRCA_std.json']:
 assert sha(ROOT/rel)==ready['sha256'][rel]
assert sha(graph)==sha(OLD/'proposed_primary_nb_BRCA.graph')
for name in ['kh_constraints.R','spatial_plus_components.py','region_order.json']:(D/name).write_bytes((OLD/name).read_bytes())
(D/'chosen_nb_BRCA.graph').write_bytes(graph.read_bytes())
m=dict(graph_sha256=sha(graph),chosen_graph_source=str(graph.resolve()),graph_role='primary_contact593_cc2',graph_components=[227,2],panel_authorized=True,authorized_panel_sha256=sha(panel),panel_source=str(panel.resolve()),panel_gate_reason='Final DATA_READY hashes checked; user/parent confirms direct KOSIS verification and independent audit passed.',original_c3_audit=str(OLD.resolve()),authorized_input_hashes={p.name:sha(p) for p in [B/'config_BRCA_std.json',B/'panel_META_BRCA_std.csv',B/'panel_INFO_BRCA_std.json',B/'sets_BRCA_std.json']})
(D/'dependency_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n')
s=(D/'kh_constraints.R').read_text();s=s.replace('  invisible(TRUE)\n}', '  invisible(TRUE)\n}')
needle='  stopifnot(digest::digest(file=panel_file,algo="sha256")==m$authorized_panel_sha256)'
s=s.replace(needle,needle+'''\n  for (fname in names(m$authorized_input_hashes))
    stopifnot(digest::digest(file=file.path(dirname(panel_file),fname),algo="sha256")==m$authorized_input_hashes[[fname]])
  stopifnot(identical(as.integer(lengths(kh_components(g)$nodes)),c(227L,2L)))''')
(D/'kh_constraints.R').write_text(s)
s=(OLD/'88k_struct_kh.R').read_text()
start=s.index('## 그래프(');end=s.index('DD<-DD[order',start)
s=s[:start]+'''## Frozen, independently reviewed primary graph override (no re-derived queen graph).
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
KH_GRAPH_FILE<-Sys.getenv("KH_GRAPH_FILE",unset=file.path(ROOT,"input_graph.graph"))
g<-inla.read.graph(KH_GRAPH_FILE)
kh_assert_dependencies(g,KH_GRAPH_FILE,G$region,
  file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),KH_HELPER_DIR)
'''+s[end:]
s=s.replace('legacy graph; component constraints','contact593 cc2 graph; component constraints')
(D/'88k_struct_kh.R').write_text(s)
a=(ROOT/'01_스크립트/88k_struct_kh.R').read_text()
(D/'88k_struct_kh.patch').write_text(''.join(difflib.unified_diff(a.splitlines(True),s.splitlines(True),fromfile='a/01_스크립트/88k_struct_kh.R',tofile='b/01_스크립트/88k_struct_kh.R')))
p=ROOT/'01_스크립트/88g_sens_suite.R';a=p.read_text();b=a.replace('c("splus","lag1","lag2","extreme","ma3")','c("splus_cc2_K23","lag1","lag2","extreme","ma3")').replace('c(splus="D1_spatial_plus",','c(splus_cc2_K23="D1_spatial_plus_cc2_K23",')
assert a!=b
needle='    if(dn(an))next;X<-rdp(sprintf("S41_panel_%s.csv",nm));'
guard='''    if(nm=="splus_cc2_K23"){
      md<-jsonlite::fromJSON(file.path(RES,"S41_splus_cc2_K23_definition.json"))
      stopifnot(BEST==md$main_set,setequal(MAINV,md$covariates),
        digest::digest(file=file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),algo="sha256")==md$panel_sha256,
        digest::digest(file=file.path(RES,sprintf("nb_%s.graph",TAG)),algo="sha256")==md$graph_sha256,
        digest::digest(file=file.path(RES,"S41_panel_splus_cc2_K23.csv"),algo="sha256")==md$output_sha256)
    }
'''
assert b.count(needle)==1;b=b.replace(needle,guard+needle)
(D/'88g_splus_suite_cc2.patch').write_text(''.join(difflib.unified_diff(a.splitlines(True),b.splitlines(True),fromfile='a/01_스크립트/88g_sens_suite.R',tofile='b/01_스크립트/88g_sens_suite.R')))
(D/'source_snapshot_hashes.json').write_text(json.dumps({str(p):sha(p) for p in [ROOT/'01_스크립트/88k_struct_kh.R',ROOT/'01_스크립트/88g_sens_suite.R',ROOT/'DATA_READY.json',graph,panel]},indent=2)+'\n')
print('Final cc2 bundle staged in owned folder only; final panel',sha(panel),'graph',sha(graph))
