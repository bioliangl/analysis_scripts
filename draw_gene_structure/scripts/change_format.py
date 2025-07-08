#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import argparse

def processGFF(GFFfile_path):
    gene_dict = {}
    trans_dict = {}
    exon_dict = {}
    CDS_dict = {}
    error = []
    with open(GFFfile_path, 'r') as GFFfile_path_file:
        for line in GFFfile_path_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chr_id, _, feature, start, end, _, strand, _, attribute = line.strip().split('\t')
            start = int(start)
            end = int(end)
            if feature == "mRNA":
                RNA_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                mRNA_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                exon_count = CDS_count = 0
                if chr_id:
                    trans_dict[mRNA_id] = (RNA_parent_id, int(start), int(end), strand, chr_id)
                else:
                    error.append(line.strip())
            elif feature == "exon":
                exon_count += 1
                exon_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                try:
                    exon_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                except Exception as e:
                    exon_id = exon_count
                if strand not in ["+", "-"]:
                    error.append(line.strip())
                    continue
                if trans_dict[exon_parent_id]:
                    if not exon_parent_id in exon_dict:
                        exon_dict[exon_parent_id] = {}
                    exon_dict[exon_parent_id][str(exon_id)] = (int(start), int(end), strand)
                else:
                    error.append(line.strip())
            elif feature == "CDS":
                CDS_count += 1
                CDS_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                try:
                    CDS_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                except Exception as e:
                    CDS_id = CDS_count
                if strand not in ["+", "-"]:
                    error.append(line.strip())
                    continue
                if trans_dict[CDS_parent_id]:
                    if not CDS_parent_id in CDS_dict:
                        CDS_dict[CDS_parent_id] = {}
                    CDS_dict[CDS_parent_id][str(CDS_id)] = (int(start), int(end), strand)
                else:
                    error.append(line.strip())
    return trans_dict, exon_dict, CDS_dict, error

def processGFF_scale(GFFfile_path):
    trans_dict = {}
    exon_dict = {}
    CDS_dict = {}
    error = []
    with open(GFFfile_path, 'r') as GFFfile_path_file:
        for line in GFFfile_path_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chr_id, _, feature, start, end, _, strand, _, attribute = line.strip().split('\t')
            start = int(start)
            end = int(end)
            if feature == "mRNA":
                RNA_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                mRNA_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                exon_count = CDS_count = 0
                mrna_start = start
                mrna_end = end
                mrna_length = end - start + 1
                mrna_start_scale = 1
                mrna_end_scale = end - start + 1
                if chr_id:
                    trans_dict[mRNA_id] = (RNA_parent_id, int(mrna_start_scale), int(mrna_end_scale), strand, chr_id)
                else:
                    error.append(line.strip())
            elif feature == "exon":
                exon_count += 1
                exon_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                try:
                    exon_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                except Exception as e:
                    exon_id = exon_count
                if strand not in ["+", "-"]:
                    error.append(line.strip())
                    continue
                if trans_dict[exon_parent_id]:
                    if not exon_parent_id in exon_dict:
                        exon_dict[exon_parent_id] = {}
                    exon_start = start - mrna_start + 1
                    exon_end = mrna_length - (mrna_end - end)
                    exon_dict[exon_parent_id][str(exon_id)] = (int(exon_start), int(exon_end), strand)
                else:
                    error.append(line.strip())
            elif feature == "CDS":
                CDS_count += 1
                CDS_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                try:
                    CDS_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                except Exception as e:
                    CDS_id = CDS_count
                if strand not in ["+", "-"]:
                    error.append(line.strip())
                    continue
                if trans_dict[CDS_parent_id]:
                    if not CDS_parent_id in CDS_dict:
                        CDS_dict[CDS_parent_id] = {}
                    CDS_start = start - mrna_start + 1
                    CDS_end = mrna_length - (mrna_end - end)
                    CDS_dict[CDS_parent_id][str(CDS_id)] = (int(CDS_start), int(CDS_end), strand)
                else:
                    error.append(line.strip())
    return trans_dict, exon_dict, CDS_dict, error

def processInterpro(interpro_path, domainDB):
    interpro_dict={}
    hsp_id = 0
    with open(interpro_path, 'r') as lines:
        for line in lines:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            gene_id, type, _, start, end, evalue, _, _, attribute = line.strip().split('\t')
            start = int(start)
            end = int(end)
            if type in domainDB:
                hsp_id += 1
                domain_name = re.search(r"Name=([^\s,;]+)", attribute).group(1)
                try:
                    domain_desc = re.search(r"signature_desc=([^,;]+)", attribute).group(1)
                except Exception as e:
                    domain_desc = domain_name
                if not gene_id in interpro_dict:
                    hsp_id = 1
                    interpro_dict[gene_id] = {}
                    interpro_dict[gene_id][hsp_id] = (start,end,evalue,domain_name,domain_desc)
                else:
                    interpro_dict[gene_id][hsp_id] = (start, end, evalue, domain_name, domain_desc)
    return interpro_dict

def domainOrdChange(cdd_dict, trans_dict, CDS_dict):
    domain_dict = {}
    for query, hsp_data in cdd_dict.items():
        domain_dict[query] = {}
        for hsp, (hsp_start, hsp_end, _, _, hsp_short_name) in hsp_data.items():
            hsp_length = (hsp_end - hsp_start + 1) * 3
            hsp_start = (hsp_start - 1) * 3 + 1
            hsp_end = hsp_end * 3
            CDS_list = []
            CDS_list_length = 0
            domain_list = []
            CDS_dict_sorted = sorted(CDS_dict[query].items(), key = lambda x:x[1][0])
            for CDS_id, (CDS_start, CDS_end, CDS_strand) in CDS_dict_sorted:
                CDS_list.append([CDS_start, CDS_end, CDS_strand])
                CDS_list_length += 1
            if trans_dict[query][3] == '+':
                for i in range(1, CDS_list_length + 1):
                    head_list = CDS_list.pop(0)
                    head_start, head_end, head_strand = head_list
                    if head_end - head_start + 1 < hsp_start:
                        hsp_start = hsp_start - (head_end - head_start + 1)
                        continue
                    if head_start + hsp_start - 1 + hsp_length - 1 > head_end:
                        hsp_length = (head_start + hsp_start - 1 + hsp_length - 1) - head_end + 1
                        domain_list.append([head_start + hsp_start - 1, head_end, head_strand, hsp_short_name])
                        hsp_start = 1
                    else:
                        domain_list.append([head_start + hsp_start - 1, head_start + hsp_start - 1 + hsp_length - 1, head_strand, hsp_short_name])
                        break
            elif trans_dict[query][3] == '-':
                for i in range(1, CDS_list_length + 1):
                    head_list = CDS_list.pop(-1)
                    head_start, head_end, head_strand = head_list
                    if head_end - head_start + 1 < hsp_start:
                        hsp_start = hsp_start - (head_end - head_start + 1)
                        continue
                    if head_end - hsp_start + 1 - hsp_length + 1 < head_start:
                        hsp_length = head_start - (head_end - hsp_start + 1 - hsp_length + 1)
                        domain_list.append([head_start, head_end - hsp_start + 1, head_strand, hsp_short_name])
                        hsp_start = 1
                    else:
                        domain_list.append([head_end - hsp_start + 1 - hsp_length + 1, head_end - hsp_start + 1, head_strand, hsp_short_name])
                        break

            domain_dict[query][hsp] = domain_list
    return domain_dict

def printoutput(trans_dict, exon_dict, CDS_dict, domain_dict):
    gff_lines = []
    trans_dict_sorted = sorted(trans_dict.items(), key = lambda x: (int(re.search(r"(\d+)", x[1][4]).group(1)), x[1][1]))
    for mRNA_id, (gene_id, mRNA_start, mRNA_end, mRNA_strand, chr_id) in trans_dict_sorted:
        gff_lines.append('\t'.join([chr_id, ".", "mRNA", str(mRNA_start), str(mRNA_end), ".", mRNA_strand, ".", mRNA_id, gene_id, "."]))

        exon_dict_sorted = sorted(exon_dict[mRNA_id].items(), key = lambda x:x[1][0])
        for exon_id, (exon_start, exon_end, exon_strand) in exon_dict_sorted:
            gff_lines.append('\t'.join([chr_id, ".", "exon", str(exon_start - mRNA_start + 1), str(exon_end - mRNA_start + 1), ".", exon_strand, ".", mRNA_id, gene_id, "."]))

        CDS_dict_sorted = sorted(CDS_dict[mRNA_id].items(), key=lambda x: x[1][0])
        for CDS_id, (CDS_start, CDS_end, CDS_strand) in CDS_dict_sorted:
            gff_lines.append('\t'.join([chr_id, ".", "CDS", str(CDS_start - mRNA_start + 1), str(CDS_end - mRNA_start + 1), ".", CDS_strand, ".", mRNA_id, gene_id, "."]))

        if mRNA_id in domain_dict:
            domain_dict_sorted = sorted(domain_dict[mRNA_id].items(), key=lambda x: x[1][0])
            for domain_id, domain_line in domain_dict_sorted:
                if mRNA_strand == "-":
                    domain_line = list(reversed(domain_line))
                for domain_start, domain_end, domain_strand, short_name in domain_line:
                    gff_lines.append('\t'.join([chr_id, ".", "domain", str(domain_start - mRNA_start + 1), str(domain_end - mRNA_start + 1), ".", domain_strand, ".", mRNA_id, gene_id, short_name]))
    return gff_lines

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="aln2structureplot")

    parser.add_argument('-g', '--gff3', required=True, help='the gff3 annotation of all gene in .aln')
    parser.add_argument('-a', '--anno', required=True, help='annotation file')
    parser.add_argument('-d', '--db', required=False, default = 'Pfam', choices = ["Pfam", "CDD", "SUPERFAMILY"], help='use annotation from which database as priority')
    parser.add_argument('-o', '--out', required=True, help='ggtranscript plot gff file')
    parser.add_argument('-s', '--scale', required=False, default = 'true', choices = ["true", "false"], help='gene position normalization')
    args = parser.parse_args()

    gff3 = args.gff3
    anno = args.anno
    out = args.out
    domainDB = args.db
    scale = args.scale
    if scale == "true":
        trans_dict, exon_dict, CDS_dict, error = processGFF_scale(gff3)
    else:
        trans_dict, exon_dict, CDS_dict, error = processGFF(gff3)
        
    anno_dict = processInterpro(anno, domainDB)
    domain_dict = domainOrdChange(anno_dict, trans_dict, CDS_dict)
    gff_lines = printoutput(trans_dict, exon_dict, CDS_dict, domain_dict)

    with open(out, 'w') as f:
        for line in gff_lines:
            f.write(line + '\n')

    with open(out + "_errorlines.txt", 'w') as e:
        for line in error:
            e.write(line + '\n')

if __name__ == "__main__":
    main()
