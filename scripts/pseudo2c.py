exec(open('pseudo2b.py').read().split("ID=(0,1,2)")[0])
# full-group orbits of non-violating pairs, as types
orbT={}; types=[]
ID=(0,1,2); I1={1:1,2:2}; E1={1:1,2:1}
full=[((1,0,2),(1,1,1),I1,E1),((0,2,1),(1,1,1),I1,E1),(ID,(-1,1,1),I1,E1),(ID,(1,1,1),{1:2,2:1},E1),(ID,(1,1,1),I1,{1:-1,2:1})]
for S in itertools.combinations(range(n),d):
    for ls in itertools.product(labels,repeat=d):
        a=tuple(zip(S,ls))
        if a in orbT or is_viol(a,V): continue
        stack=[a]; orbT[a]=len(types); mem=[a]
        while stack:
            b_=stack.pop()
            for g in full:
                c=act(g,b_)
                if c not in orbT: orbT[c]=len(types); mem.append(c); stack.append(c)
        types.append(mem)
print(f"{len(types)} orbit types of non-violating pairs\n")
def desc(a):
    (i,li),(j,lj)=a; x,y=free[i],free[j]
    if leq(x,y) or leq(y,x): rel='comparable'
    elif leq(neg(x),y) or leq(y,neg(x)): rel='anti-comp'
    else: rel='unrelated'
    return f"{rel:10s} ranks={tuple(sorted((size(x),size(y))))} {'same-mag' if abs(li)==abs(lj) else 'diff-mag'} {'same-sign' if (li>0)==(lj>0) else 'opp-sign'}"
# knockout test: is there a solution avoiding type t entirely?
necessary=[]
for t,mem in enumerate(types):
    memset=set(mem)
    x,_=solve_restricted([],lambda a,ms=memset: a not in ms, f"exclude type {t:2d} [{desc(mem[0])}] (size {len(mem):3d})")
    if x is None: necessary.append(t)
print("\nNECESSARY types (no solution without them):")
for t in necessary: print("  ",desc(types[t][0]))
