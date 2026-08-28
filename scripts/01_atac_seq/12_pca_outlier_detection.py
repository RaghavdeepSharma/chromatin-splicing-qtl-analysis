#!/usr/bin/env python3
"""Flag extreme PCA outliers from the count matrix.
This is an audit step because Calderon mentions PCA outlier exclusion but not a fully specified algorithm.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ap=argparse.ArgumentParser()
ap.add_argument('--counts', required=True)
ap.add_argument('--sample-list', required=True)
ap.add_argument('--out-prefix', required=True)
ap.add_argument('--mad-z', type=float, default=6)
args=ap.parse_args()

# Handle variation in nucleoATAC get_count output formats.
# User may adapt first columns if needed.
df=pd.read_csv(args.counts, sep='\t', comment='#')
samples=[x.strip() for x in open(args.sample_list) if x.strip()]
missing=[s for s in samples if s not in df.columns]
if missing:
    raise SystemExit('Counts file does not contain expected sample columns: '+','.join(missing[:10]))
X=df[samples].T.astype(float)
X=np.log2(X+1)
X=StandardScaler(with_mean=True, with_std=True).fit_transform(X)
pca=PCA(n_components=2).fit_transform(X)
out=pd.DataFrame({'biosample_id':samples, 'PC1':pca[:,0], 'PC2':pca[:,1]})
for pc in ['PC1','PC2']:
    med=out[pc].median(); mad=np.median(np.abs(out[pc]-med)) or 1e-9
    out[pc+'_mad_z']=0.6745*(out[pc]-med)/mad
out['pca_outlier']=(out['PC1_mad_z'].abs()>args.mad_z) | (out['PC2_mad_z'].abs()>args.mad_z)
out.to_csv(args.out_prefix+'.pca_outliers.tsv', sep='\t', index=False)
out.loc[out.pca_outlier, ['biosample_id']].to_csv(args.out_prefix+'.pca_outlier_ids.txt', index=False, header=False)
print(out.to_string(index=False))
