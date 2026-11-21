import sys
from setup import instance
m, k = int(sys.argv[1]), int(sys.argv[2])
cx, N, labels, edges, flags = instance(m, k)
q = len(labels)
out = [f"#define N {N}", f"#define Q {q}", f"#define NE {len(edges)}", f"#define NF {len(flags)}", f"#define M {m}",
       "static const int LAB[Q]={" + ",".join(map(str, labels)) + "};",
       "static const int EA[NE]={" + ",".join(str(a) for (a, b) in edges) + "};",
       "static const int EB[NE]={" + ",".join(str(b) for (a, b) in edges) + "};",
       "static const int ES[NE]={" + ",".join(str(s) for s in edges.values()) + "};",
       "static const int FI[NF][M]={" + ",".join("{" + ",".join(str(i) for i, _ in f) + "}" for f in flags) + "};",
       "static const int FS[NF][M]={" + ",".join("{" + ",".join(str(s) for _, s in f) + "}" for f in flags) + "};"]
open("inst.h", "w").write("\n".join(out) + "\n")
print(N, q, len(edges), len(flags))
