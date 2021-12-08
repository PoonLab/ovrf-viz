# Generate random noise on the orfs coordinates
# Modified from https://dzone.com/articles/simple-csv-transformations
    
from pathlib import Path
import csv
import random

source_path = Path("total_orfs.csv")
target_path = Path(source_path.stem + "_1").with_suffix('.csv')

def randcoord(coordinate):
    return str(random.randrange(coordinate+1))

def to_multi_coords(multi_coords_str):
    return [coords.split(":") for coords in multi_coords_str.split(";")]

def to_multi_coords_str(multi_coords):
    coords = []
    for multi_coord in multi_coords:
        coords.append(":".join([randcoord(int(coord)) for coord in multi_coord]))
    return ";".join(coords)

def transform(row, row_number):
    coords = to_multi_coords(row["coords"])
    row["coords"] = to_multi_coords_str(coords)
    return row

with source_path.open() as source_file:
    with target_path.open('w', newline='') as target_file:
        reader = csv.DictReader(source_file)
        columns = reader.fieldnames
        writer = csv.DictWriter(target_file, columns)
        writer.writeheader()
        n_rows = 451228
        percentage = 0.1
        candidates = dict.fromkeys(random.sample(range(1, n_rows), int(n_rows*percentage)))
        for i, row in enumerate(reader):
            if i % 10000 == 0:
                print(".", end="", flush=True)
            new_row = row
            if i in candidates:
                new_row = transform(row, i)

            writer.writerow(new_row)