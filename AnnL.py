import sys
import re
import os
import subprocess
from Bio import SeqIO
'''
读入ensembl注释文件
'''


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


'''
读入gffcompare后的string merge的gtf文件
'''


def read_assembly_gtf(assembly_gtf_file):
    merge_annotated = {}
    merge_unannotated = {}
    gene = ""
    with open(assembly_gtf_file, "r") as file:
        for line in file:
            line = line.strip()
            array = line.split("\t")
            if array[2] == "transcript":
                '''
                class_code "i" 也存在可能
                '''
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


def calculate_exon_num(assembly_gtf_file):
    exon_num = 0
    transcript_exon = []
    transcript_exon_num = {}
    with open(assembly_gtf_file, "r") as file:
        for line in file:
            line = line.strip()
            array = line.split("\t")
            if array[2] == "transcript":
                if exon_num > 0:
                    transcript_exon.append(exon_num)
                transcript_id = re.search(r'transcript_id [A-Za-z0-9".-_]+', array[8])[0].split("\"")[-2]
                transcript_exon.append(transcript_id)
                exon_num = 0
            if array[2] == "exon":
                exon_num += 1

    for i in range(0, len(transcript_exon)-1, 2):
        transcript_exon_num[transcript_exon[i]] = int(transcript_exon[i+1])

    return transcript_exon_num


'''
对这些文件进行biotype分类，并计算数目
'''


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


'''
对于分编码蛋白，获取位置信息，小于长度阈值的会被舍弃，阈值默认为200
'''


def non_protein_coding(merge_unannotated, length_threshold, work_dir, transcript_exon_num, exonlimit):
    uncoding = []
    for gene, transcripts in merge_unannotated.items():
        for transcript in transcripts.keys():
            position_start = transcripts[transcript][0]
            position_end = transcripts[transcript][1]
            chromosome = transcripts[transcript][2]
            if int(position_end) - int(position_start) >= length_threshold and transcript_exon_num[transcript] > int(exonlimit):
                uncoding.append(chromosome)
                uncoding.append(transcript)
                uncoding.append(position_start)
                uncoding.append(position_end)
    #print(os.path.join(work_dir, "uncoding.bed.file"))
    bed_file = open(os.path.join(work_dir, "uncoding.bed"), "a")
    bed_file_name = os.path.join(work_dir, "uncoding.bed")
    for seq in range(0, len(uncoding), 4):
        bed_file.write('\t'.join([uncoding[seq], str(int(uncoding[seq+2])-1), uncoding[seq+3], uncoding[seq+1]]))
        bed_file.write("\n")
    bed_file.close()
    #print('\t'.join([uncoding[seq], uncoding[seq+1], uncoding[seq+2]]))
    return bed_file_name


'''
调用bedtools获取序列
'''


def get_non_protein_coding_seq(work_dir, fasta, bed_file_name):
    bed_file_path = bed_file_name
    fasta_name = '.'.join([os.path.splitext(bed_file_name)[0], "fasta"])
    bed_out = fasta_name
    bed_tool_command = ''.join(["bedtools getfasta", " -fi ", fasta, " -bed ", bed_file_path,
                                " -fo ", os.path.join(work_dir, fasta_name)])
    print("="*50)
    print("[INFO] BEGIN geting sequence...")
    result = subprocess.getoutput(bed_tool_command)
    if "bedtools: not found" in result:
        print("[ERROR] ERROR, can't find bedtools")
        print("[ERROR] Please add the path of bedtools to envirment")
    else:
        print("[INFO] FINISH")
    print("="*50)
    return bed_out


'''
对序列进行CPC2分析，分析编码能力
'''


def CPC2_analysis(work_dir, bed_out):
    cpc_input = os.path.join(work_dir, bed_out)
    which_cpc = subprocess.getoutput("which CPC2.py")
    cpc_output = os.path.join(work_dir, "cpc.out")
    cpc_command = ' '.join(['python', which_cpc, '-i', cpc_input, '-o', cpc_output])
    # print(cpc_command)
    print("="*50)
    print("[INFO] BEGIN CPC2 analysis")
    result = subprocess.getoutput(cpc_command)
    print(result)
    print("[INFO] REMOVE coding sequence")
    uncoding_list = []
    with open(cpc_output, "r") as file:
        for line in file:
            if "#" in line:
                continue
            elif line.strip().split("\t")[-1] == "coding":
                continue
            else:
                uncoding_list.append(line)
    cpc_uncoding = os.path.join(work_dir, "cpc.uncoding")
    file = open(cpc_uncoding, "a")
    for i in uncoding_list:
        file.write(i)
    file.close()

    print("[INFO] FINISH")
    print("="*50)


'''
对编码蛋白进行CNCI分析，进行编码能力分析
'''


def CNCI_analysis(work_dir, bed_out, threads):
    CNCI_input = os.path.join(work_dir, bed_out)
    which_CNCI = subprocess.getoutput("which CNCI.py")
    CNCI_output = os.path.join(work_dir, "CNCI")
    CNCI_command = ' '.join(["python", which_CNCI, "-f", CNCI_input, "-o", CNCI_output, "-m ve", "-p", str(threads)])
    print("="*50)
    print("[INFO] BEGIN CNCI analysis...")
    print("[INFO] ANALYSIS is running..., please be patient")
    result = subprocess.getoutput(CNCI_command)
    print(result)
    print("[INFO] FINISH")
    print("[INFO] REMOVE coding sequence")
    uncoding_list = []
    #print(os.path.join(CNCI_output, "CNCI.index"))
    with open(os.path.join(CNCI_output, "CNCI.index"), "r") as file:
        for line in file:
            array = line.strip().split("\t")
            if "Transcript ID" in line:
                continue
            if array[1] == "noncoding":
                uncoding_list.append(line)
    file = open(os.path.join(work_dir, "CNCI.uncoding"), "a")
    for i in uncoding_list:
        file.write(i)
    file.close()
    print("[INFO] FINISH")
    print("="*50)


'''
对分析的结果寻找交集，并提取交集的序列
'''


def find_intersection(work_dir):
    file_list = []
    print("="*50)
    print("[INFO] FIND INTERSECTION")
    for root, dirs, files in os.walk(work_dir):
        for file in files:
            if os.path.splitext(file)[1] == ".uncoding":
                list_name = os.path.splitext(file)[0]
                list_name = []
                with open(os.path.join(root, file), "r") as fi:
                    for line in fi:
                        array = line.strip().split("\t")
                        list_name.append(array[0])
                file_list.append(list_name)
    '''
    三个列表找交集的好方法，也可以用‘&=’
    '''
    intersection = set(file_list[0]).intersection(*file_list[1:])
#    print(intersection)
    file = open(os.path.join(work_dir, "intersection.uncoding.bed"), "a")
    bed_file_name = os.path.join(work_dir, "intersection.uncoding.bed")
    for gene in intersection:
        file.write('\t'.join([gene.split(":")[0], gene.split(":")[1].split("-")[0], gene.split(":")[1].split("-")[1]]))
        file.write("\n")
    file.close()
    print("[INFO] FINISH")
    print("="*50)
    return bed_file_name


def Pfam_detection(work_dir, bed_out, PfamA, threads):
    Pfam_input = os.path.join(work_dir, "Pfam.input.fasta")
    '''
    将DNA序列翻译成蛋白质
    '''
    file = open(Pfam_input, "a")
    for seq_record in SeqIO.parse(os.path.join(work_dir, bed_out), "fasta"):
        dna_id = seq_record.id
        DNA = seq_record.seq
        AND = DNA.reverse_complement()
        file.write(''.join([">", dna_id]))
        file.write("\n")
        file.write(str(DNA.translate()))
        file.write("\n")
        file.write(''.join([">", dna_id, "rev"]))
        file.write("\n")
        file.write(str(AND.translate()))
        file.write("\n")
    file.close()
    '''

    '''
    which_pfam = subprocess.getoutput("which pfam_scan.pl")
    which_hmmer = subprocess.getoutput("which hmmer")
    Pfam_out = os.path.join(work_dir, "Pfam.result")
    Pfam_dir = PfamA
    Pfam_command = ' '.join(["perl", which_pfam, "-fasta", Pfam_input, "-dir",
                             Pfam_dir, "-outfile", Pfam_out, "-cpu", str(threads)])
    print("="*50)
    print("[INFO] BEGIN Pfam detection...")
    print("THIS will take a LONG time")
    result = subprocess.getoutput(Pfam_command)
    print(result)
    print("[INFO] FINISH Pfam detection...")

    Pfam_annotated = []
    with open(Pfam_out, "r") as file:
        for line in file:
            array = line.strip().split("\t")
            if "#" in line:
                continue
            else:
                Pfam_annotated.append(array[0].strip("rev"))
    print("[INFO] FILTER the result...")
    print("[INFO] FINISH")

    return Pfam_annotated


def finall_result(Pfam_annotated, bed_file_name, work_dir):
    print("="*50)
    finall_result = open(os.path.join(work_dir, "AnnL.finall.out"), "a")
    with open(bed_file_name, "r") as file:
        for line in file:
            array = line.strip().split("\t")
            lncRNA = ''.join([array[0], ":", array[1], "-", array[2]])
            if lncRNA in set(Pfam_annotated):
                continue
            else:
                finall_result.write(line)
    finall_result.close()
    print("[INFO] ALL ANALYSIS FINISH")


'''
参数设置
'''


def get_args():
    import argparse
    usage = "python3 AnnL.py -G ensembl.gtf -g merge.gtf -f reference.fasta <--CPC [1,0] --CNCI [1,0] --threads num --Pfam [1,0] --PfamA /path/PfamA/ > -o /path/out/ "
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('--gtf-file', "-G", required=True, dest="ensembl_gtf", help="the path of ensembl gtf file")
    parser.add_argument('--merge_gtf', "-g", required=True, dest="merge_gtf",
                        help="the path of merge gtf file, produced by Stringtir -merge or Cuffmerge")
    parser.add_argument('--llimit', '-l', default="200", type=int, dest="length",
                        help="the length threshold of LncRNA, default = 200")
    parser.add_argument('--exonlimit', '-el', default=1, type=int, dest="exon_num",
                        help="the exon number threshold of LncRNA, default = 1")
    parser.add_argument('--output', '-o', dest="work_dir", help="the directory of results and temfile")
    parser.add_argument('--fasta', '-f', dest="fasta", required=True, help="the path of fasta file")
    parser.add_argument('--CPC', dest="cpc", default=1, type=int,
                        help="Boolen value, if False, the programe will not perform CPC2 analysi")
    parser.add_argument('--CNCI', dest="cnci", default=1, type=int,
                        help="Boolen value, if False, the programe will not perform CNCI analysi")
    parser.add_argument('--threads', '-t', dest="threads", default=4, type=int, help="threads number")
    parser.add_argument('--Pfam', dest="Pfam_d", default=1, type=int,
                        help="Boolen value, if False, the programe will not perform Pfam detection")
    parser.add_argument('--PfamA', dest="PfamA", default=" ",  help="the path of PfamA directory")
    args = parser.parse_args()
    return args


def main():
    '''
    获取变量的值
    '''
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
    threads = args.threads
    PfamA = args.PfamA

    exonlimit = args.exon_num

    '''
    调用函数
    '''
    annotation, biotypes = read_annotated_gtf(gtf_file)
#   print(annotation)
    merge_annotated, merge_unannotated = read_assembly_gtf(assembly_gtf_file)
#    print(merge_unannotated)
    biotype_distribution(biotypes, annotation, merge_annotated, merge_unannotated)
    transcript_exon_num = calculate_exon_num(assembly_gtf_file)
    bed_file_name = non_protein_coding(merge_unannotated, length_threshold, work_dir, transcript_exon_num, exonlimit)
    bed_out = get_non_protein_coding_seq(work_dir, fasta, bed_file_name)

    '''
    设定是否进行CPC和CNCI分析，1为进行分析，0为拒绝
    '''
    if args.cpc > 0:
        CPC2_analysis(work_dir, bed_out)
    if args.cnci > 0:
        CNCI_analysis(work_dir, bed_out, threads)

    '''
    交集
    '''
    intersection_bed = find_intersection(work_dir)
    intersection_fasta = get_non_protein_coding_seq(work_dir, fasta, intersection_bed)  # 提取交集的序列
    #intersection_fasta = "/home/insilicon/Desktop/aNNl/intersection.uncoding.fasta"
    if args.Pfam_d > 0:
        if args.PfamA == " ":
            print("If using Pfam detection, please provid the path of Pfam datadir")
        else:
            Pfam_annotated = Pfam_detection(work_dir, intersection_fasta, PfamA, threads)
    finall_result(Pfam_annotated, intersection_bed, work_dir)


if __name__ == "__main__":
    main()
