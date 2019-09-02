"""
find_ovrfs.csv has duplicate entries for overlaps between reading frames that
occur in multiple segments because of "introns".
"""

from csv import DictReader, DictWriter
handle = open('../data/find_ovrfs.csv')
reader = DictReader(handle)

last_state = None
last_overlap = 0
last_shift = None

outfile = open('../data/ovrfs-reduced.csv', 'w')
writer = DictWriter(outfile, fieldnames=[
    'accn', 'prod1', 'loc1', 'dir1', 'prod2', 'loc2', 'dir2', 'seqlen1', 'seqlen2', 'overlap', 'shift'
])

# read first row
prev_row = next(reader)
prev_row['overlap'] = int(prev_row['overlap'])

test_fields = ['accn', 'loc1', 'dir1', 'seqlen1', 'loc2', 'dir2', 'seqlen2']

for row in reader:
    row['overlap'] = int(row['overlap'])

    is_different = False
    for field in test_fields:
        if row[field] != prev_row[field]:
            # current row is different
            is_different = True
            break

    if is_different:
        # output previous row
        writer.writerow(prev_row)
        prev_row = row
    else:
        # duplicate row, accumulate overlap
        prev_row['overlap'] += row['overlap']


