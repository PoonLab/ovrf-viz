"""
- how many genes per virus have overlaps?
- what is the overall coding proportion in overlaps? (nt and %)
- what is the type of overlap (0, +1, -1, etc)
- is the amino acid composition of overlapping regions different?
"""
from csv import DictReader

def get_records(file):
    reader = DictReader(open(file))
    this_accn = None
    part = []
    for row in reader:
        if row['accno'] != this_accn:
            # new record
            if this_accn is not None:
                yield this_accn, part
            this_accn = row['accno']
            part = [row]
        else:
            part.append(row)
    yield this_accn, part  # last record

def parse_coords(coords):
    res = []
    for interval in coords.split(';'):
        left, right = interval.split(':')
        res.append((int(left), int(right)))
    return(res)

def find_ovrfs(accn, records, outfile):
    """
    Examine coordinates of ORFs for overlaps.
    :param records:
    :return: list of overlapping regions, if any
    """
    n = len(records)
    for i in range(n):
        rec1 = records[i]
        c1 = parse_coords(rec1['coords'])
        for j in range(i):
            rec2 = records[j]
            c2 = parse_coords(rec2['coords'])
            for l1, r1 in c1:
                for l2, r2 in c2:
                    left = max(l1, l2)
                    right = min(r1, r2)
                    overlap = right - left
 
                    if overlap > 0:
                        outfile.write('{},"{}","{}",{},{},{}\n'.format(accn, rec1['product'], rec2['product'], overlap, rec1['strand'], rec2['strand']))


outfile = open('../data/find_ovrfs.csv', 'w')
outfile.write('accn,prod1,prod2,overlap,strand1,strand2\n')

for accn, records in get_records('../data/orfs-fixed.csv'):
    #print(accn)
    find_ovrfs(accn, records, outfile)

