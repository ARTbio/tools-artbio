import sys

F = open(sys.argv[1], 'r')
lines = F.readlines()
lines = [line[:-1] for line in lines]
for line in sorted(lines):
    print(line)
F.close()
