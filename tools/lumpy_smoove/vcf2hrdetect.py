import sys

handle = open(sys.argv[1], 'r')
vcfdict = dict()
tabdict = dict()
for line in handle:
    if line[0] == "#":
        continue
    else:
        tabline = line[:-1].split("\t")
        vcfdict[tabline[2]] = tabline
for id in vcfdict.keys():
    if "_1" in id:
        newid = id[:-2]
        pointbreak = vcfdict[id][4]
        if "]" in pointbreak:
            coordbreak = pointbreak.split("]")[1].split(":")[1]
            chrom = pointbreak.split("]")[1].split(":")[0]
        elif "[" in pointbreak:
            coordbreak = pointbreak.split("[")[1].split(":")[1]
            chrom = pointbreak.split("[")[1].split(":")[0]
        if vcfdict[id][0] == chrom:
            tabdict[newid] = [chrom, vcfdict[id][1], chrom, coordbreak, "INV"]
        else:
            tabdict[newid] = [vcfdict[id][0], vcfdict[id][1],
                              chrom, coordbreak, "TRA"]
for id in list(vcfdict):
    if "_" in id:
        del vcfdict[id]
for id in vcfdict.keys():  # only sv that are not of type TRA or INV
    chr1 = vcfdict[id][0]
    chr2 = vcfdict[id][0]
    pos1 = vcfdict[id][1]
    pos2 = vcfdict[id][7].split("END=")[1].split(";")[0]
    type = vcfdict[id][7].split("SVTYPE=")[1].split(";")[0]
    tabdict[id] = [chr1, pos1, chr2, pos2, type]
out = open(sys.argv[2], 'w')
out.write("chr1\tpos1\tchr2\tpos2\ttype\n")
for key in tabdict:
    line = "\t".join(tabdict[key]) + "\n"
    out.write(line)
