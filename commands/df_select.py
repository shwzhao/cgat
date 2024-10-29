#!/usr/bin/env python3

def parse_hierarchy(file_path):
    hierarchy = {}
    # 'clade': {}, 'order': {}, 'family': {}

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                name, rank, parent = parts
                if rank == "clade":
                    clade_name = name
                    hierarchy[clade_name] = {}
                elif rank == "order":
                    order_name = name
                    hierarchy[clade_name][order_name] = []
                elif rank == "family":
                    hierarchy[clade_name][order_name].append(name)
    return hierarchy

file_path = "/pfs/stor10/users/home/s/shwzhao/apg4.2.txt"

aa = parse_hierarchy(file_path)
print(aa)