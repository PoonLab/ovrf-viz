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

def find_ovrfs(records):
    """
    Examine coordinates of ORFs for overlaps.
    :param records:
    :return: list of overlapping regions, if any
    """
    n = len(records)
    for i in range(n):
        c1 = parse_coords(records[i]['coords'])
        for j in range(i):
            c2 = parse_coords(records[j]['coords'])
            for l1, r1 in c1:
                for l2, r2 in c2:
                    if r1 > l2 and l1 < r2:
                        print('{}--{}  {}--{}'.format(l1, r1, l2, r2))


for accn, records in get_records('orfs-fixed.csv'):
    find_ovrfs(records)
    break