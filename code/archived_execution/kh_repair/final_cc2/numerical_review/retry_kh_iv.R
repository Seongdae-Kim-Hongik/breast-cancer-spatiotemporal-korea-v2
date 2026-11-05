source('BRCA_재투고준비_20261104/kh_repair/final_cc2/numerical_review/common.R')
z<-readRDS(file.path(D,'KH_IV_matrices.rds'));KHQ<-inla.as.sparse(z$Q);KHC<-list(A=z$A,e=rep(0,z$rankdef));stopifnot(z$rankdef==255)
ex<-"f(st,model='generic0',Cmatrix=KHQ,rankdef=255,constr=FALSE,extraconstr=KHC,hyper=pcp,values=1:3206)"
old<-readRDS(file.path(N,'initial_snapshot/S43_fit_KH_typeIV.rds'))
ci<-list(control.vb=list(emergency=50),tolerance=1e-4,tolerance.f=1e-6,tolerance.g=1e-3,tolerance.x=1e-3,restart=1,stupid.search=TRUE)
cmA<-list(theta=old$mode$theta,x=old$mode$x,restart=TRUE,fixed=FALSE)
a<-fit_once('retry_A',file.path(N,'retry'),extra=ex,ci=ci,cm=cmA,A=z$A)
# Deliberately dispersed starting values, not a change to likelihood or priors.
th<-old$mode$theta;th[1]<-th[1]-2;th[2]<-th[2]+.15;th[3]<-th[3]-.15;th[4]<-th[4]+.2;th[5:6]<-th[5:6]-1
cmB<-list(theta=th,restart=TRUE,fixed=FALSE)
b<-fit_once('retry_B',file.path(N,'retry'),extra=ex,ci=ci,cm=cmB,A=z$A)
LOG('BOUNDED_KH_RETRIES_DONE two attempted; no adoption decision here')
