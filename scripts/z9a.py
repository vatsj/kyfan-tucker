exec(open('s3.py').read().split("ID=tuple(range(M))")[0])
import pickle
ID=(0,1,2,3); I1={1:1,2:2,3:3}; E1={1:1,2:1,3:1}
gens=[((1,2,0,3),(1,1,1,1),I1,E1), (ID,(1,1,1,1),{1:2,2:3,3:1},E1)]   # Z_3 x Z_3, odd order => averaging valid over F_2
d=3; t=time.time(); orb={}; norb=0
for S in itertools.combinations(range(n),d):
    for ls in itertools.product(labels,repeat=d):
        a=tuple(zip(S,ls))
        if a in orb or not is_viol(a): continue
        stack=[a]; orb[a]=norb
        while stack:
            b=stack.pop()
            for g in gens:
                c=act(g,b)
                if c not in orb: orb[c]=norb; stack.append(c)
        norb+=1
print(f"Z3xZ3 orbits at degree 3: {norb}  [{time.time()-t:.0f}s]"); sys.stdout.flush()
# integer keys: key(alpha) = sum over (i,l): base^k ... use (i,label_index) packed
li={l:k for k,l in enumerate(labels)}
def key(a): return sum((i*6+li[l])*(240**k) for k,(i,l) in enumerate(a))
keys=np.array(sorted(key(a) for a in orb),dtype=np.int64)
oid=np.array([orb[a] for a in sorted(orb,key=key)],dtype=np.int64)
triples=np.array(list(itertools.combinations(range(n),3)),dtype=np.int64)
rng=np.random.default_rng(11); nrows=int(norb*1.15)
def row(Lidx):   # Lidx: label indices per vertex
    k=(triples[:,0]*6+Lidx[triples[:,0]]) + (triples[:,1]*6+Lidx[triples[:,1]])*240 + (triples[:,2]*6+Lidx[triples[:,2]])*240**2
    pos=np.searchsorted(keys,k); pos[pos>=len(keys)]=0
    hit=keys[pos]==k
    return np.bincount(oid[pos[hit]],minlength=norb)&1
t=time.time()
W=(norb+1+63)//64
Mat=np.zeros((nrows,W),dtype=np.uint64)
for r in range(nrows):
    bits=row(rng.integers(0,6,n)); bits=np.append(bits,1)   # rhs=1 in last column
    Mat[r]=np.packbits(bits.astype(np.uint8),bitorder='little').view(np.uint8).tobytes().ljust(W*8,b'\0') and np.frombuffer(np.packbits(bits.astype(np.uint8),bitorder='little').tobytes().ljust(W*8,b'\0'),dtype=np.uint64)
np.save('z9M.npy',Mat); pickle.dump((norb,nrows,keys,oid),open('z9meta.pkl','wb'))
print(f"system: {nrows} rows x {norb} unknowns (+rhs), built [{time.time()-t:.0f}s]")
