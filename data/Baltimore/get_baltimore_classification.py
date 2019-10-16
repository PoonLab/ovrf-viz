# Clasify viruses according to their type of molecule (dsDNA, ssDNA, dsRNA, RT-RNA, (+)ss-RNA, (-)ssRNA, (circular)ssRNA)
import csv
import pandas as pd
#dir = '../data/Baltimore'

filename_list = ("baltimore_cirucular_ssRNA.txt", "baltimore_RTviruses.txt",
                 "baltimore_ssRNA_plus.txt", "baltimore_dsDNA.txt", "baltimore_ssDNA.txt",
                 "baltimore_dsRNA.txt", "baltimore_ssRNA_minus.txt")

# Create dictionary with baltimore classification as keys, and list of families as values
virus_dict = {}
for filename in filename_list:
        families = []
        with open(filename) as fp:
            line = fp.readline()
            cnt = 1
            while line:
                if 'Baltimore' in line:
                    my_line = line.strip()
                    key = my_line.split(" ")[1]
                elif 'dae' in line:
                    my_line = line.strip()
                    families.append(my_line)

                line = fp.readline()
                cnt += 1
                # elif 'unassigned' in line:
                #     print(line)
                # TODO: include viruses with unsassigned families
        virus_dict[key] = families


# Import virus table
species_file_name = "virus_full_info.txt"
#species_file = pd.read_csv(species_file_name)

species_file = pd.read_csv(species_file_name,converters={"Taxonomy": lambda x: x.strip("[]").replace("'", "").split(", ")})
print(len(species_file))
subset = species_file[0:100]
baltimore_class = []
viral_family = []


for row in species_file['Taxonomy']:
    #values = species_file['Taxonomy'][row]
    #print(values[0])
    classification = ''
    vf = ''  # Viral Family
    for virus_classification in row:
        for baltimore, family_list in virus_dict.items():
            for family in family_list:
                if virus_classification == family:
                    classification = baltimore
                    vf = family
                elif 'dae' in virus_classification:    #print("I founded this family: {} for this virus: {}, on this baltimore: {}".format(family, species_file['Genome'][i], baltimore))
                    vf = virus_classification
    if classification:
        baltimore_class.append(classification)
        viral_family.append(vf)
    elif vf:  # Taxonomy has family classification but is not on ViralZone
        baltimore_class.append("Unknown")
        viral_family.append(vf)
    else:
        baltimore_class.append('Unknown')
        viral_family.append('Unknown')

species_file["baltimore.class"] = baltimore_class
species_file["family"] = viral_family

species_file.to_csv(r'species_file_baltimore2.csv')