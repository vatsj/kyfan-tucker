import numpy as np, pickle, time, os, sys
norb,nrows,keys,oid=pickle.load(open('z9meta.pkl','rb'))
if os.path.exists('z9state.pkl'):
    M=np.load('z9M.npy'); r,c0=pickle.load(open('z9state.pkl','rb'))
else:
    M=np.load('z9M.npy'); r,c0=0,0
t=time.time(); R=M.shape[0]; one=np.uint64(1)
for c in range(c0,norb):
    w,b=c//64,c%64
    colbits=(M[r:,w]>>np.uint64(b))&one
    nz=np.nonzero(colbits)[0]
    if len(nz)==0: continue
    p=r+nz[0]
    if p!=r: M[[r,p]]=M[[p,r]]
    mask=((M[:,w]>>np.uint64(b))&one).astype(bool); mask[r]=False
    idx=np.nonzero(mask)[0]
    if len(idx): M[idx]^=M[r]
    r+=1
    if time.time()-t>36000:
        np.save('z9M.npy',M); pickle.dump((r,c+1),open('z9state.pkl','wb'))
        print(f"checkpoint: column {c+1}/{norb}, rank so far {r}  [{time.time()-t:.0f}s]"); sys.exit(0)
rhs=(M[r:,norb//64]>>np.uint64(norb%64))&one
print(f"DONE: rank(A)={r}; system {'INCONSISTENT -> NO F_2 degree-3 certificate on S^3' if rhs.any() else 'consistent -> F_2 degree-3 certificate EXISTS on S^3'}")
np.save('z9M.npy',M); pickle.dump((r,norb),open('z9state.pkl','wb'))
