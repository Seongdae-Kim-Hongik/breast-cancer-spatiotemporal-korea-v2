# Pure matrix construction; never fits a model. Spatial order = saved graph order.
# time-major vectorization: st = (time - 1) * n + area.
kh_components <- function(g) {
  n <- g$n; nb <- lapply(g$nbs, as.integer)
  stopifnot(length(nb) == n)
  for (i in seq_len(n)) {
    stopifnot(length(nb[[i]]) == g$nnbs[i], !anyDuplicated(nb[[i]]),
              all(nb[[i]] >= 1 & nb[[i]] <= n), !i %in% nb[[i]])
    for (j in nb[[i]]) stopifnot(i %in% nb[[j]])
  }
  cc <- integer(n); k <- 0L
  for (i in seq_len(n)) if (cc[i] == 0L) {
    k <- k + 1L; todo <- i; cc[i] <- k
    while (length(todo)) {
      v <- todo[1]; todo <- todo[-1]
      new <- nb[[v]][cc[nb[[v]]] == 0L]
      cc[new] <- k; todo <- c(todo, new)
    }
  }
  list(id=cc, n=k, nodes=split(seq_len(n), cc))
}

kh_scale_connected <- function(R) {
  # Explicit per-component calls preserve installed INLA scaling convention,
  # including its numerical ridge for the constrained inverse (not added to Q).
  m <- nrow(R)
  if (m < 2L) stop("Singleton graph component requires a separately reviewed prior")
  INLA::inla.scale.model(R, constr=list(A=matrix(1, 1, m), e=0))
}

kh_matrices <- function(type, g, Tn) {
  type <- match.arg(type, c("II", "III", "IV"))
  stopifnot(length(Tn)==1L, is.finite(Tn), Tn>=2L, Tn==as.integer(Tn))
  n <- g$n; cc <- kh_components(g); cnum <- cc$n
  if (any(lengths(cc$nodes) < 2L)) stop("Singleton policy not authorized")
  C <- Matrix::sparseMatrix(i=cc$id, j=seq_len(n), x=1,
                           dims=c(cnum,n))
  B <- Matrix::sparseMatrix(i=rep(seq_len(Tn-1L), each=2L),
                           j=as.vector(rbind(seq_len(Tn-1L),2:Tn)),
                           x=rep(c(-1,1),Tn-1L), dims=c(Tn-1L,Tn))
  Rt0 <- crossprod(B)
  Rs0 <- Matrix::Diagonal(n,g$nnbs) - Matrix::sparseMatrix(
    i=rep(seq_len(n),g$nnbs), j=unlist(g$nbs), x=1, dims=c(n,n))
  Rs <- Rs0; Rt <- Rt0
  # II needs only temporal scaling; III needs only spatial scaling.
  if (type != "III") Rt <- kh_scale_connected(Rt0)
  if (type != "II") for (ii in cc$nodes)
    Rs[ii,ii] <- kh_scale_connected(Rs0[ii,ii,drop=FALSE])
  At <- kronecker(matrix(1,1,Tn), Matrix::Diagonal(n))
  Ac <- kronecker(Matrix::Diagonal(Tn), C)
  if (type == "II") {
    Q <- kronecker(Rt,Matrix::Diagonal(n)); A <- At; rd <- n
  } else if (type == "III") {
    Q <- kronecker(Matrix::Diagonal(Tn),Rs); A <- Ac; rd <- cnum*Tn
  } else {
    Q <- kronecker(Rt,Rs)
    # Omit the LAST time point for EVERY component, not one global row.
    A <- rbind(At,Ac[seq_len(cnum*(Tn-1L)),,drop=FALSE])
    rd <- n+cnum*(Tn-1L)
  }
  stopifnot(nrow(A)==rd,ncol(A)==n*Tn)
  list(Q=Q,A=as.matrix(A),rankdef=rd,rank=n*Tn-rd,
       Rs=Rs,Rt=Rt,Rs0=Rs0,Rt0=Rt0,C=C,components=cc,type=type,n=n,Tn=Tn)
}

kh_check_panel <- function(d, n, vars) {
  # All currently authorized KH/KH4x variants retain the full 2011–2024 panel.
  stopifnot(identical(sort(unique(as.integer(d$year))),2011:2024),
            all(d$time == d$year-2010), all(d$yiid == d$time),
            nrow(d)==n*14L, all(d$area %in% seq_len(n)))
  st <- (d$time-1L)*n+d$area
  stopifnot(!anyDuplicated(st),identical(sort(as.integer(st)),seq_len(n*14L)),
            all(is.finite(d$E)),all(d$E>0),all(is.finite(d$cases)))
  for (v in vars) if (!v %in% names(d) || !is.numeric(d[[v]]) ||
                      any(!is.finite(d[[v]])))
    stop("Missing/nonfinite KH covariate: ",v,"; no silent row deletion")
  invisible(TRUE)
}

kh_fit_repaired <- function(type, d, vars) {
  kh_check_panel(d,n,vars)
  z <- kh_matrices(type,g,max(d$time))
  d$st <- (d$time-1L)*n+d$area
  # Retain the original formula environment, likelihood, E, main effects,
  # covariates, fixed/PC priors and numerical retry settings.
  assign("KHQ",INLA::inla.as.sparse(z$Q),envir=.GlobalEnv)
  assign("KHC",list(A=z$A,e=rep(0,z$rankdef)),envir=.GlobalEnv)
  ex <- sprintf(paste0("f(st,model='generic0',Cmatrix=KHQ,rankdef=%d,",
    "constr=FALSE,extraconstr=KHC,hyper=pcp,values=1:%d)"),z$rankdef,n*z$Tn)
  LOG("KH %s repaired: n=%d T=%d components=%d rank=%d rankdef=%d",
      type,n,z$Tn,z$components$n,z$rank,z$rankdef)
  safe(vars,data=d,extra=ex,env=.GlobalEnv)
}

kh_assert_dependencies <- function(g, graph_file, regions, panel_file, helper_dir) {
  m <- jsonlite::fromJSON(file.path(helper_dir,"dependency_manifest.json"))
  stopifnot(digest::digest(file=graph_file,algo="sha256")==m$graph_sha256)
  ref <- INLA::inla.read.graph(file.path(helper_dir,"chosen_nb_BRCA.graph"))
  stopifnot(g$n==ref$n,all(vapply(seq_len(g$n),function(i)
    identical(sort(as.integer(g$nbs[[i]])),sort(as.integer(ref$nbs[[i]]))),logical(1))))
  expected <- jsonlite::fromJSON(file.path(helper_dir,"region_order.json"))
  stopifnot(identical(enc2utf8(as.character(regions)),enc2utf8(expected)))
  if (!isTRUE(m$panel_authorized) || is.null(m$authorized_panel_sha256))
    stop("Final panel independent verification/baseline regeneration pending; no fits")
  stopifnot(digest::digest(file=panel_file,algo="sha256")==m$authorized_panel_sha256)
  for (fname in names(m$authorized_input_hashes))
    stopifnot(digest::digest(file=file.path(dirname(panel_file),fname),algo="sha256")==m$authorized_input_hashes[[fname]])
  stopifnot(identical(as.integer(lengths(kh_components(g)$nodes)),c(227L,2L)))
  invisible(TRUE)
}
