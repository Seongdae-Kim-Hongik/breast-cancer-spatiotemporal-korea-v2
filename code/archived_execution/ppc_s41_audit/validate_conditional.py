"""Conditional 2024 is a separate fit, not an S41_part analysis."""
def validate_status(j):
 issues=[]
 if j.get('analysis')!='conditional_2024' or j.get('success') is not True:issues.append('Conditional2024 missing successful fit identity')
 for k,n in [('n_rows',3206),('n_observed',2977),('n_heldout',229)]:
  if j.get(k)!=n:issues.append('Conditional2024 incorrect '+k)
 for k in ['ok','fixed_finite','hyper_finite','waic_finite']:
  if j.get(k) is not True:issues.append('Conditional2024 failed '+k)
 for k in ['mode_status','misc_mode_status']:
  if k not in j:issues.append('Conditional2024 missing '+k);continue
  v=j[k]
  if k=='mode_status' and v is None:issues.append('Conditional2024 missing primary mode status')
  if v is not None and any(x!=0 for x in (v if isinstance(v,list) else [v])):issues.append('Conditional2024 nonzero '+k)
 if 'internal_warnings' not in j or j.get('internal_warnings'):issues.append('Conditional2024 internal warnings absent or require review')
 if j.get('strategy')!='default' and not j.get('verified_isolated_retry'):issues.append('Conditional2024 nondefault strategy requires review')
 return issues
