import sys
from itertools import islice

block = 4
headerdic = {}

for line in open(sys.argv[1],"r"):
    headerdic[line[:-1]]=1

with open(sys.argv[3],"w") as out:
    with open(sys.argv[2], "r") as f:
        while True:
            next_n_lines = list(islice(f, 4))
            if not next_n_lines:
                break
            try:
                key = next_n_lines[0].split(" ")[0][1:]
                headerdic[key]
                out.write("".join(next_n_lines))
            except KeyError:
                continue
