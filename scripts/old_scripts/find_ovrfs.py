"""
- how many genes per virus have overlaps?
- what is the overall coding proportion in overlaps? (nt and %)
- what is the type of overlap (0, +1, -1, etc)
- is the amino acid composition of overlapping regions different?
"""

from csv import DictReader

def get_records(file):
    """
    Stream through CSV file, grouping rows by accession
    :param file:
    :return:
    """
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
        print(interval)
        left, right = interval.split(':')
        res.append((int(left), int(right)))
    return(res)

def find_ovrfs(accn, records, outfile):
    """
    Examine coordinates of ORFs for overlaps.
    :param records:
    :return: list of overlapping regions, if any
    """
    print(accn)
    n = len(records)
    for i in range(n):
        rec1 = records[i]
        print(rec1['coords'], type(rec1['coords']), len(rec1))         
        c1 = parse_coords(rec1['coords'])
        xl1 = min([l for l, r in c1])  # extreme left
        len1 = sum([r-l for l, r in c1])  # total CDS length
        dir1 = rec1['strand']

        for j in range(i):
            rec2 = records[j]
            c2 = parse_coords(rec2['coords'])
            xl2 = min([l for l, r in c2])
            len2 = sum([r-l for l, r in c2])
            dir2 = rec2['strand']

            for l1, r1 in c1:
                for l2, r2 in c2:
                    left = max(l1, l2)
                    right = min(r1, r2)
                    overlap = (right - left) +1

                    shift = abs(xl1-xl2) % 3
                    if dir1 != dir2:
                        shift = '-'+str(shift)
                    else:
                        shift = '+'+str(shift)

                    if overlap > 0:
                        outfile.write('{},"{}",{},{},"{}",{},{},{},{},{},{}\n'.format(
                            accn, rec1['product'], xl1, dir1, rec2['product'], xl2, dir2,
                            len1, len2, overlap, shift
                        ))


outfile = open('find_ovrfs_out.csv', 'w')
outfile.write('accn,prod1,loc1,dir1,prod2,loc2,dir2,seqlen1,seqlen2,overlap,shift\n')

for accn, records in get_records('/home/lmunoz/Projects/ovrf-review/dataset_2020/total_orfs.csv'):
    #print(accn)
    find_ovrfs(accn, records, outfile)
