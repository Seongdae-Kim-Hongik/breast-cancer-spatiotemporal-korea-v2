# ONLY audited original88k nb/struct specifications; no KH/conf/period/combine.
source('BRCA_재투고준비_20261104/kh_repair/final_cc2/numerical_review/common.R')
OUT<-file.path(N,'nb_struct');dir.create(file.path(OUT,'graphs'),showWarnings=FALSE)
# canonical contact593 g was hash/order checked by common.R. Never rebuild queen.
G<-st_read('${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp',quiet=TRUE)
stopifnot(identical(enc2utf8(as.character(G$region)),enc2utf8(regions)))
cen<-st_coordinates(st_centroid(st_geometry(G)))
# Do not assign an undocumented CRS; KNN reproduces original coordinate metric.
runone<-function(id,sp=SPf(),note='canonical_contact593',gf=file.path(D,'chosen_nb_BRCA.graph')){
 m<-fit_once(id,OUT,sp=sp,graphnote=note,graph_file=gf)
 if(is.null(m))m<-fit_once(paste0(id,'_EB_retry'),OUT,sp=sp,ci=list(control.vb=list(emergency=50),strategy='gaussian',int.strategy='eb'),graphnote=note,graph_file=gf)
 invisible(m)
}
runone('NB_main_contact593')
for(k in 3:5){
 nbk<-make.sym.nb(knn2nb(knearneigh(cen,k=k)));stopifnot(all(card(nbk)>0))
 gf<-file.path(OUT,'graphs',paste0('knn',k,'.graph'));nb2INLA(gf,nbk);gk<-inla.read.graph(gf)
 stopifnot(gk$n==n,all(vapply(1:n,function(i)identical(sort(as.integer(gk$nbs[[i]])),sort(as.integer(nbk[[i]]))),logical(1))))
 write_json(list(policy='Original88k symmetric centroid KNN; deliberate alternative topology; coordinate metric unchanged and CRS unverified',k=k,edges=sum(gk$nnbs)/2,components=lengths(gk$cc$nodes),graph_sha256=digest::digest(file=gf,algo='sha256')),paste0(gf,'.json'),auto_unbox=TRUE,pretty=TRUE)
 runone(paste0('NB_knn',k),sp=SPf('gk'),note=paste0('symmetric_centroid_knn',k),gf=gf)
}
runone('STRUCT_icar_besag',sp="f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)")
runone('STRUCT_iid_only',sp="f(area,model='iid',hyper=pcp)")
runone('STRUCT_no_spatial',sp='')
pcC<-pcU(.5,.01,.5,.5);pcA<-pcU(1,.01,.5,.5);pcB<-pcU(2,.05,.5,.5)
runone('STRUCT_pc_phi_flat',sp=SPf(hyper='pcC'))
runone('STRUCT_pc_relaxA',sp=SPf(hyper='pcA'))
runone('STRUCT_pc_relaxB',sp=SPf(hyper='pcB'))
LOG('NB_STRUCT_DONE 10 specifications sequential; no KH or legacy-queen branch')
