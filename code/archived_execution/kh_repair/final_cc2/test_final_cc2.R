suppressPackageStartupMessages({library(INLA);library(Matrix);library(jsonlite)})
D<-normalizePath('BRCA_재투고준비_20261104/kh_repair/final_cc2');source(file.path(D,'kh_constraints.R'));inla.setOption(num.threads='1:1')
m<-fromJSON(file.path(D,'dependency_manifest.json'));g<-inla.read.graph(file.path(D,'chosen_nb_BRCA.graph'));n<-g$n;regions<-fromJSON(file.path(D,'region_order.json'))
kh_assert_dependencies(g,file.path(D,'chosen_nb_BRCA.graph'),regions,m$panel_source,D)
stopifnot(sum(g$nnbs)==1186,identical(as.integer(lengths(kh_components(g)$nodes)),c(227L,2L)))
gmvariance<-function(Q){ee<-eigen(as.matrix(Q),symmetric=TRUE);keep<-ee$values>max(abs(ee$values))*1e-9;v<-rowSums(sweep(ee$vectors[,keep,drop=FALSE]^2,2,ee$values[keep],'/'));exp(mean(log(v)))}
records<-list()
for(ty in c('II','III','IV')){
 z<-kh_matrices(ty,g,14L)
 es<-eigen(as.matrix(z$Rs),symmetric=TRUE,only.values=TRUE)$values;et<-eigen(as.matrix(z$Rt),symmetric=TRUE,only.values=TRUE)$values
 ev<-switch(ty,II=rep(et,each=n),III=rep(es,14),IV=as.vector(outer(et,es)))
 rankq<-sum(ev>max(abs(ev))*1e-9);ea<-eigen(tcrossprod(z$A),symmetric=TRUE,only.values=TRUE)$values;ranka<-sum(ea>max(abs(ea))*1e-9)
 qa<-max(abs(z$Q%*%t(z$A)));stopifnot(rankq==z$rank,ranka==z$rankdef,qa<1e-9)
 chol<-Cholesky(forceSymmetric(z$Q+crossprod(Matrix(z$A,sparse=TRUE))),LDL=FALSE);stopifnot(all(is.finite(diag(as(chol,'sparseMatrix')))))
 if(ty=='IV'){
  Ac<-as.matrix(kronecker(Diagonal(14),z$C));err<-max(abs(Ac-Ac%*%t(z$A)%*%solve(tcrossprod(z$A))%*%z$A));stopifnot(err<1e-10)
 }
 scalecheck<-list()
 if(ty!='II'){scalecheck$spatial<-vapply(z$components$nodes,function(ii)gmvariance(z$Rs[ii,ii,drop=FALSE]),numeric(1));stopifnot(max(abs(scalecheck$spatial-1))<1e-5)}
 if(ty!='III'){scalecheck$temporal<-gmvariance(z$Rt);stopifnot(abs(scalecheck$temporal-1)<1e-5)}
 saveRDS(z,file.path(D,paste0('KH_',ty,'_matrices.rds')))
 records[[ty]]<-list(rank=rankq,rankdef=z$rankdef,rankA=ranka,max_Q_At=qa,positive_definite_after_constraints=TRUE,scaled_generalized_variance=scalecheck)
}
stopifnot(records$II$rankdef==229,records$III$rankdef==28,records$IV$rankdef==255)
panel<-read.csv(m$panel_source,check.names=FALSE);panel$area<-match(panel$fbd_region,regions);panel$time<-panel$year-2010;panel$yiid<-panel$time
cfg<-fromJSON(file.path(dirname(m$panel_source),'config_BRCA_std.json'));vars<-unique(c(cfg$sets_after_vif[[cfg$rules$main_set]],cfg$forced))
LOG<-function(...)invisible(NULL)
safe<-function(vars,data,extra,env){stopifnot(identical(data$cases,panel$cases),identical(data$E,panel$E),identical(data$year,panel$year),identical(data[,vars],panel[,vars]));list(extra=extra)}
for(ty in c('II','III','IV')){m0<-kh_fit_repaired(ty,panel,vars);stopifnot(grepl(paste0('rankdef=',records[[ty]]$rankdef),m0$extra,fixed=TRUE))}
reject<-function(expr)inherits(tryCatch({force(expr);NULL},error=function(e)e),'error')
stopifnot(reject(kh_check_panel(panel[panel$year>=2016,],n,vars)))
bad<-panel;bad$E[1]<-NA;stopifnot(reject(kh_check_panel(bad,n,vars)))
stopifnot(reject(kh_assert_dependencies(g,file.path(D,'chosen_nb_BRCA.graph'),rev(regions),m$panel_source,D)))
oldg<-inla.read.graph(file.path(dirname(D),'chosen_nb_BRCA.graph'))
stopifnot(reject(kh_assert_dependencies(oldg,file.path(dirname(D),'chosen_nb_BRCA.graph'),regions,m$panel_source,D)))
invisible(parse(file.path(D,'88k_struct_kh.R')))
write_json(list(passed=TRUE,fits_run=0,panel_hash_gate_passed=TRUE,graph_hash_gate_passed=TRUE,input_hash_gates_passed=TRUE,graph_components=c(227,2),tests=records,final_panel_mock_fit_passed=TRUE,partial_period_rejected=TRUE,wrong_graph_rejected=TRUE,wrong_region_order_rejected=TRUE,missing_E_rejected=TRUE,drop_in_parsed=TRUE),file.path(D,'validation.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('FINAL cc2 PASS: rankdefs229/28/255; final panel/hash gates; mocks only, no fits\n')
