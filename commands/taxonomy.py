#!/usr/bin/env python3

def get_phylo(species_names):
    from Bio import Entrez
    try:
        Entrez.email = "shwzhao997@gmail.com"
        results = {}
        for species_name in species_names:
            with Entrez.esearch(db="Taxonomy", term=species_name) as handle:
                record = Entrez.read(handle)
                if "IdList" in record and record["IdList"]:
                    id = record["IdList"][0]
                    with Entrez.efetch(db="Taxonomy", id=id, retmode="xml") as handle:
                        records = Entrez.read(handle)
                        if records:
                            results[species_name] = records[0]
                        else:
                            print(species_name + ": No records found")
                else:
                    print(species_name + ": No records found")
        return results
    except Exception as e:
        print("An error occurred: {}".format(e))

def get_rank_scientific_names(results):
    rank_scientific_names = {}
    for species_name, data in results.items():
        rank_scientific_names[species_name] = {}
        rank_scientific_names[species_name]['ScientificName'] = data['ScientificName']
        rank_scientific_names[species_name]['TaxId'] = data['TaxId']
        for item in data["LineageEx"]:
            rank = item.get("Rank")
            scientific_name = item.get("ScientificName")
            if rank and scientific_name:
                rank_scientific_names[species_name][rank] = scientific_name
    return rank_scientific_names


def setup_parser(parser):
    taxonomy_parser = parser.add_parser('taxonomy', help='Get species taxonomy')
    # Add command specific arguments
    taxonomy_parser.add_argument('-n', '--species_name', help='Species name. | For example: "Arabidopsis thaliana"')
    taxonomy_parser.add_argument('-f', '--species_file', help='File with Species names')
    taxonomy_parser.add_argument('-o', '--output_file', help='Path to output file')

    return taxonomy_parser

def read_species_names_from_file(file_path):
    species_names = []
    with open(file_path, 'r') as file:
        for line in file:
            species_names.append(line.strip())
    return species_names

def run(args):
    # if args.species_name and args.species_file:
    #     args.error("Options -n and -f are mutually exclusive. Please choose only one.")

    # species_names = ["Arabidopsis thaliana", "Solanum lycopersicum"]

    if args.species_file:
        species_names = read_species_names_from_file(args.species_file)
    elif args.species_name:
        species_names = []
        species_names.append(args.species_name)
    else:
        print("\nWrong species name parameter setting!\n")
        import sys
        sys.exit()

    phylo_info = get_phylo(species_names)
    if phylo_info:
        taxonomy_info = get_rank_scientific_names(phylo_info)
        import pandas as pd
        df = pd.DataFrame.from_dict(taxonomy_info, orient='index')
        if args.output_file:
            df.to_csv(args.output_file, sep='\t', index=False)
        else:
            import sys
            df.to_csv(sys.stdout, sep='\t', index=False)