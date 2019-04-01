import sys
import re
import os
import subprocess


def read_annotated_gtf(gtf_file):
    annotation = {}
    with open(gtf_file, "r") as file:
        for line in file:
            if "#" in line:
                continue
            else:
                array = line.strip().split("\t")
                if array[2] == "gene":
                    try:
                        gene_name = re.search(r'gene_name [A-Za-z0-9"._-]*', array[8]).group().split("\"")[-2]
                    except AttributeError:
                        gene_name = re.search(r'gene_id [A-Za-z0-9"._-]*', array[8]).group().split("\"")[-2]
                    gene_biotype = re.search(r'gene_biotype [A-Za-z"._-]*', array[8]).group().split("\"")[-2]
                    annotation[gene_name] = gene_biotype
    biotypes = []
    for transcript, biotype in annotation.items():
        biotypes.append(biotype)
    return annotation, biotypes


def read_assembly_gtf(assembly_gtf_file):
    merge_annotated = {}
    merge_unannotated = {}
    gene = ""
    with open(assembly_gtf_file, "r") as file:
        for line in file:
            line = line.strip()
            array = line.split("\t")
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
                        merge_unannotated[gene][transcript_id].append(
                            position_start)
                        merge_unannotated[gene][transcript_id].append(
                            position_end)
                        merge_unannotated[gene][transcript_id].append(array[0])
                    else:
                        merge_unannotated[gene][transcript_id] = []
                        merge_unannotated[gene][transcript_id].append(
                            position_start)
                        merge_unannotated[gene][transcript_id].append(
                            position_end)
                        merge_unannotated[gene][transcript_id].append(array[0])
                else:
                    gene_name = re.search(r'gene_name [A-Za-z0-9"._-]+', array[8])[0].split("\"")[-2]
                    transcript_id = re.search(r'transcript_id [A-Za-z0-9".-_]+', array[8])[0].split("\"")[-2]
                    position_start = array[3]
                    position_end = array[4]
                    #                class_code = re.search(r'class_code [=ckmnjeosmiyp"]*', array[8])[0].split("\"")[-2]
                    if gene_name != gene:
                        gene = gene_name
                        merge_annotated[gene] = {}
                        merge_annotated[gene][transcript_id] = []
                        merge_annotated[gene][transcript_id].append(position_start)
                        merge_annotated[gene][transcript_id].append(position_end)
                    else:
                        merge_annotated[gene][transcript_id] = []
                        merge_annotated[gene][transcript_id].append(position_start)
                        merge_annotated[gene][transcript_id].append(position_end)
#                    print(merge_annotated)
    return merge_annotated, merge_unannotated


def biotype_distribution(biotypes, annotation, merge_annotated, merge_unannotated):
    print("%s\t%s\t%s" % ("biotype", "annotated", "unannotated"))
#    print(annotation)
    biotype_distribution = {}
    for biotype in set(biotypes):
        biotype_distribution[biotype] = []
    for gene_id in merge_annotated.keys():
        biotype_distribution[annotation[gene_id]].append(gene_id)
    for biotype, number in biotype_distribution.items():
        print("%s\t%d" % (biotype, len(number)))
    print("%s\t%d" % ("unannotated", len(list(merge_unannotated.keys()))))


def non_protein_coding(merge_unannotated, length_threshold, work_dir):
    uncoding = []
    for gene, transcripts in merge_unannotated.items():
        for transcript in transcripts.keys():
            position_start = transcripts[transcript][0]
            position_end = transcripts[transcript][1]
            chromosome = transcripts[transcript][2]
            if int(position_end) - int(position_start) >= length_threshold:
                uncoding.append(chromosome)
                uncoding.append(transcript)
                uncoding.append(position_start)
                uncoding.append(position_end)
    #print(os.path.join(work_dir, "uncoding.bed.file"))
    bed_file = open(os.path.join(work_dir, "uncoding.bed.file"), "a")
    for seq in range(0, len(uncoding), 4):
        bed_file.write('\t'.join([uncoding[seq], str(int(uncoding[seq+2])-1), uncoding[seq+3], uncoding[seq+1]]))
        bed_file.write("\n")
    bed_file.close()
#        print('\t'.join([uncoding[seq], uncoding[seq+1], uncoding[seq+2]]))


def get_non_protein_coding_seq(work_dir, fasta):
    bed_file_path = os.path.join(work_dir, "uncoding.bed.file")
    bed_out = os.path.join(work_dir, "uncoding.fasta")
    bed_tool_command = ''.join(["bedtools getfasta", " -fi ", fasta, " -bed ", bed_file_path,
                                " -fo ", os.path.join(work_dir, "uncoding.fasta")])

    print("BEGIN geting sequence...")
    result = subprocess.getoutput(bed_tool_command)
    if "bedtools: not found" in result:
        print("ERROR, can't find bedtools")
        print("Please add the path of bedtools to envirment")
    else:
        print("FINISH")

    return bed_out


def CPC2_analysis(work_dir, bed_out):
    cpc_input = os.path.join(work_dir, bed_out)
    which_cpc = subprocess.getoutput("which CPC2.py")
    cpc_output = os.path.join(work_dir, "cpc.out")
    cpc_command = ' '.join(['python', which_cpc, '-i', cpc_input, '-o', cpc_output])
#    print(cpc_command)
    print("BEGIN CPC2 analysis")
    result = subprocess.getoutput(cpc_command)
    print(result)
    print("REMOVE coding sequence")
    uncoding_list = []
    with open(cpc_output, "r") as file:
        for line in file:
            if line.strip().split("\t")[-1] == "coding":
                continue
            else:
                uncoding_list.append(line)
    cpc_uncoding = os.path.join(work_dir, "cpc.uncoding")
    file = open(cpc_uncoding, "a")
    for i in uncoding_list:
        file.write(i)
    file.close()

    print("FINISH CPC2 analysis")


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf-file', "-G", required=True, dest="ensembl_gtf", help="the path of ensembl gtf file")
    parser.add_argument('--merge_gtf', "-g", required=True, dest="merge_gtf",
                        help="the path of merge gtf file, produced by Stringtir -merge or Cuffmerge")
    parser.add_argument('--llimit', '-l', default="200", type=int, dest="length",
                        help="the length threshold of LncRNA, default = 200")
    parser.add_argument('--output', '-o', dest="work_dir", help="the directory of results and temfile")
    parser.add_argument('--fasta', '-f', dest="fasta", required=True, help="the path of fasta file")
    parser.add_argument('--cpc2', dest="cpc", default=True,
                        help="Boolen value, if False, the programe will not perform CPC2 analysi")
    args = parser.parse_args()
    return args


def main():
    args = get_args()
#   print(args)
    gtf_file = args.ensembl_gtf
    assembly_gtf_file = args.merge_gtf
    length_threshold = args.length
    if args.work_dir:
        work_dir = args.work_dir
    else:
        work_dir = args.getcwd()
    fasta = args.fasta
    annotation, biotypes = read_annotated_gtf(gtf_file)
#   print(annotation)
    merge_annotated, merge_unannotated = read_assembly_gtf(assembly_gtf_file)
#    print(merge_unannotated)
    biotype_distribution(biotypes, annotation, merge_annotated, merge_unannotated)
    non_protein_coding(merge_unannotated, length_threshold, work_dir)
    bed_out = get_non_protein_coding_seq(work_dir, fasta)
    if args.cpc:
        CPC2_analysis(work_dir, bed_out)


if __name__ == "__main__":
    main()
