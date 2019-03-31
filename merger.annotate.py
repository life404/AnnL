import sys
import re
import os

def read_annotated_gtf(gtf_file):
    annotation = {}
    with open (gtf_file, "r") as file:
        for line in file:
            if "#" in line:
                continue
            else:
                array = line.strip().split("\t")
                if array[2] == "gene":
                    annotation[array[8].split(";")[2]] = array[8].split(";")[-2]
    biotypes = []
    for transcript, biotype in annotation.items():
        biotypes.append(biotype)
    return annotation, biotypes

def read_assembly_gtf(assembly_gtf_file):
    for line in file:
        line = line.strip()
        array = line.split("\t")
        merge_annotated = {}
        merge_unannotated = {}
        gene = ""
        if array[2] == "transcript":
            if "class_code \"u\"" in array[8] or "class_code \"r\"" in array[8]:
                gene_id = re.search(r'gene_id [A-Za-z0-9".-_]+', array[8])[0].split("\"")[-2]
                transcript_id = re.search(r'transcript_id [A-Za-z0-9".-_]+', array[8])[0].split("\"")[-2]
                position_start = array[3]
                position_end = array[4]
                if gene_id != gene:
                    gene = gene_id
                    merge_unannotated[gene] = {}
                    merge_unannotated[gene][transcript_id] = []
                    merge_unannotated[gene][transcript_id].append(position_start)
                    merge_unannotated[gene][transcript_id].append(position_end)
                else:
                    merge_unannotated[gene][transcript_id] = []
                    merge_unannotated[gene][transcript_id].append(position_start)
                    merge_unannotated[gene][transcript_id].append(position_end)
            else:
                gene_name = re.search(r'gene_name [A-Za-z0-9"._-]+', array[8])[0].split("\"")[-2]
                transcript_id = re.search(r'transcript_id [A-Za-z0-9".-_]+', array[8])[0].split("\"")[-2]
                position_start = array[3]
                position_end = array[4]
#                class_code = re.search(r'class_code [=ckmnjeosmiyp"]*', array[8])[0].split("\"")[-2]
                if gene_id != gene:
                    gene = gene_id
                    merge_annotated[gene] = {}
                    merge_annotated[gene][transcript_id] = []
                    merge_annotated[gene][transcript_id].append(position_start)
                    merge_annotated[gene][transcript_id].append(position_end)
                else:
                    merge_annotated[gene][transcript_id] = []
                    merge_annotated[gene][transcript_id].append(position_start)
                    merge_annotated[gene][transcript_id].append(position_end)
    return merge_annotated, merge_unannotated

def biotype_distribution(biotypes, annotation, merge_annotated, merge_unannotated):
    biotype_distribution = {}
    for biotype in set(biotypes):
        biotype_distribution[biotype] = []
    for gene_id in merge_annotated.keys():
        biotype_distribution[annotation[gene_id]].append(gene_id)
    for biotype, number  in biotype_distribution.items():
        print("%s\t%d\n" %(biotype, len(number)))
    print("%s\t%d\n" %("unannotated", len(list(merge_unannotated.keys()))))


def main():
    gtf_file = sys.argv[1]
    assembly_gtf_file = sys.argv[2]
    annotation, biotypes = read_annotated_gtf(gtf_file)
    merge_annotated, merge_unannotated = read_assembly_gtf(assembly_gtf_file)
    biotype_distribution(biotypes, annotation, merge_annotated, merge_unannotated)

if __name__ == "__main__":
    main()
   
