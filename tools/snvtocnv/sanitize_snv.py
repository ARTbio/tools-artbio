import sys

handle = open(sys.argv[1], 'r')
out = open(sys.argv[2], 'w')
for line in handle:
    if line[0] == '#':
        out.write(line)
        continue
    linelist = line[:-1].split('\t')
    refcol = linelist[0].split('chr')
    infocol = linelist[7].split('INDEL')
    if len(infocol) > 1:
        continue
    if len(refcol) > 1:
        refcol = refcol[1]
    else:
        refcol = refcol[0]
    if refcol not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                      '11', '12', '13', '14', '15', '16', '17', '18', '19',
                      '20', '21', '22', 'X', 'Y']:
        continue
    else:
        linelist[0] = refcol
        out.write('\t'.join(linelist))
