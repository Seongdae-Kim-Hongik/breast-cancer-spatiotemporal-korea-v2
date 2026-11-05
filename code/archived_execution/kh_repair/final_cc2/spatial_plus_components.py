"""Explicit component-adjusted graph-filter sensitivity; no outcome fitting.
Definition: residualize every specified exposure/control on all component
indicators and K strictly positive low-frequency Laplacian eigenvectors,
pooled over the supplied full balanced panel; standardize residuals ddof=0.
Not year-specific residualization, not a full Spatial+ spline method.
"""
from pathlib import Path
import numpy as np

def mm(a,b):
    # Explicit accumulation avoids platform BLAS floating-point flag artifacts.
    return np.einsum("ij,jk->ik",a,b,optimize=False)

def read_graph(path):
    lines=Path(path).read_text().splitlines(); n=int(lines[0]); adj=np.zeros((n,n)); seen=set()
    for line in lines[1:]:
        if not line.strip(): continue
        t=list(map(int,line.split())); i=t[0]-1
        if i in seen or not 0<=i<n or t[1]!=len(t)-2: raise ValueError('Bad graph row')
        seen.add(i); js=np.array(t[2:],dtype=int)-1
        if len(set(js))!=len(js) or np.any(js<0) or np.any(js>=n) or i in js: raise ValueError('Bad edge')
        adj[i,js]=1
    if len(seen)!=n or not np.array_equal(adj,adj.T): raise ValueError('Incomplete/asymmetric graph')
    return adj

def component_labels(adj):
    n=len(adj); labels=np.full(n,-1,dtype=int); c=0
    for i in range(n):
        if labels[i]>=0: continue
        labels[i]=c; queue=[i]
        while queue:
            v=queue.pop()
            for j in np.flatnonzero(adj[v]):
                if labels[j]<0: labels[j]=c; queue.append(j)
        c+=1
    return labels

def spatial_basis(adj,K=None):
    n=len(adj)
    if adj.shape!=(n,n) or not np.array_equal(adj,adj.T) or np.any(np.diag(adj)) or not np.all(np.isin(adj,[0,1])): raise ValueError('Need binary undirected graph')
    labels=component_labels(adj); c=labels.max()+1
    C=np.eye(c)[labels]; C/=np.sqrt(C.sum(axis=0))
    L=np.diag(adj.sum(axis=1))-adj
    w,V=np.linalg.eigh(L); tol=1e-10*max(1.0,float(np.max(np.abs(w))))
    if np.min(w)<-tol or np.count_nonzero(np.abs(w)<=tol)!=c: raise ValueError('Nullity inconsistent with graph')
    positive=np.flatnonzero(w>tol); K=int(round(.1*n)) if K is None else K
    if not isinstance(K,(int,np.integer)) or not 1<=K<=len(positive): raise ValueError('Invalid K')
    # Do not arbitrarily split a repeated eigenspace at the requested boundary.
    if K<len(positive) and abs(w[positive[K]]-w[positive[K-1]])<=tol:
        raise ValueError('K splits tied eigenspace; review K explicitly')
    U=V[:,positive[:K]]; B=np.column_stack([C,U])
    if not np.allclose(mm(B.T,B),np.eye(c+K),atol=1e-9): raise ValueError('Basis not orthonormal')
    return B,dict(n=n,components=int(c),component_sizes=np.bincount(labels).tolist(),K=int(K),null_tolerance=tol,selected_eigenvalues=w[positive[:K]].tolist(),cutoff_gap=float(w[positive[K]]-w[positive[K-1]]) if K<len(positive) else None,rank=c+K,labels=labels.tolist())

def residualize(values,area0,basis):
    values=np.asarray(values,dtype=float); idx=np.asarray(area0)
    if values.ndim==1: values=values[:,None]
    if not np.issubdtype(idx.dtype,np.integer) or len(idx)!=len(values) or np.any(idx<0) or np.any(idx>=len(basis)): raise ValueError('Invalid graph row mapping')
    if not np.isfinite(values).all(): raise ValueError('No silent missing-value removal')
    counts=np.bincount(idx,minlength=len(basis))
    if np.any(counts==0) or len(set(counts))!=1: raise ValueError('Require balanced full-region panel')
    Ba=basis[idx]; beta,_,rank,_=np.linalg.lstsq(Ba,values,rcond=None)
    if rank!=Ba.shape[1]: raise ValueError('Projection design rank deficiency')
    raw=values-mm(Ba,beta); sd=raw.std(axis=0,ddof=0)
    if np.any(sd<=1e-12): raise ValueError('A covariate is annihilated; do not silently drop/redefine it')
    standardized=(raw-raw.mean(axis=0))/sd
    return standardized,raw,sd
