# Clasify viruses according to their type of molecule (dsDNA, ssDNA, dsRNA, RT-RNA, (+)ss-RNA, (-)ssRNA, (circular)ssRNA)
import os
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

