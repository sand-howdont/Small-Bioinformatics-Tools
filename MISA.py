#!/usr/bin/env python3
import sys
import os
import re
import datetime

VERSION = "2.2 (Python Port - Deduplicated & Renumbered)"

def print_help_and_exit(script_name):
    help_text = """
DESCRIPTION: Tool for the identification and localization of
              (I)  perfect microsatellites as well as
              (II) compound microsatellites (two individual microsatellites,
                   disrupted by a certain number of bases)

SYNTAX:    misa.py <FASTA file>

    <FASTAfile>    Single file in FASTA format containing the sequence(s).

    In order to specify the search criteria, an additional file containing
    the microsatellite search parameters is required named "misa.ini". 
    If "misa.ini" is missing, the following default settings are applied:
      definition(unit_size,min_repeats):          1-10 2-6 3-5 4-5 5-5 6-5
      interruptions(max_difference_for_2_SSRs):   100
      GFF:                                        true
"""
    print(help_text)
    sys.exit(1)

def load_specifications():
    # 默认配置
    typrep = {1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5}
    amb = 100
    gff = True

    if os.path.exists("misa.ini"):
        with open("misa.ini", "r") as f:
            for line in f:
                if re.match(r'^def', line, re.IGNORECASE):
                    pairs = re.findall(r'(\d+)\s*-\s*(\d+)|(\d+)\s+(\d+)', line)
                    if pairs:
                        typrep = {}
                        for p in pairs:
                            p = [x for x in p if x] # 移除空匹配
                            typrep[int(p[0])] = int(p[1])
                elif re.match(r'^int', line, re.IGNORECASE):
                    match = re.search(r'(\d+)', line)
                    if match:
                        amb = int(match.group(1))
                elif re.match(r'^GFF', line, re.IGNORECASE):
                    if re.search(r'false', line, re.IGNORECASE):
                        gff = False
                    elif re.search(r'true', line, re.IGNORECASE):
                        gff = True
    return typrep, amb, gff

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(dna))

def main():
    if len(sys.argv) == 0 or len(sys.argv) < 2:
        print_help_and_exit(sys.argv[0])
    
    fasta_path = sys.argv[1]
    if fasta_path in ["-help", "--help", "-h"]:
        print_help_and_exit(sys.argv[0])

    if not os.path.exists(fasta_path):
        print(f"\nError: FASTA file '{fasta_path}' doesn't exist !\n")
        sys.exit(1)

    loctime = datetime.date.today().strftime('%Y-%m-%d')
    typrep, amb, gff = load_specifications()
    typ_sorted = sorted(typrep.keys())

    out_misa = None
    if not gff:
        out_misa = open(f"{fasta_path}.misa", "w")
        out_misa.write("ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n")

    max_repeats = 1
    min_repeats = min(typrep.values()) if typrep else 1000
    count_motifs = {}  
    count_motifs_repeats = {} 
    count_class = {t: 0 for t in typ_sorted}
    
    number_sequences = 0
    size_sequences = 0
    ssr_containing_seqs = {}
    ssr_in_compound = 0

    with open(fasta_path, "r") as f:
        content = f.read()
    
    entries = content.split(">")
    for entry in entries:
        if not entry.strip():
            continue
        
        lines = entry.split("\n", 1)
        id_line = lines[0]
        seq_body = lines[1] if len(lines) > 1 else ""
        
        seq = re.sub(r'[\d\s>]', '', seq_body)
        if gff:
            seq_id = id_line.split()[0] if id_line.split() else "Unknown"
        else:
            seq_id = id_line.strip().replace(" ", "_")
        
        number_sequences += 1
        size_sequences += len(seq)

        nr = 0
        start_dict, end_dict, motif_dict, repeats_dict, len_dict = {}, {}, {}, {}, {}
        
        for motiflen in typ_sorted:
            minreps = typrep[motiflen] - 1
            search_pat = re.compile(r'(([acgt]{' + str(motiflen) + r'})\2{' + str(minreps) + r',})', re.IGNORECASE)
            
            for m in search_pat.finditer(seq):
                ssr_string = m.group(1).upper()
                motif = m.group(2).upper()
                
                redundant = False
                for j in range(motiflen - 1, 0, -1):
                    if motiflen % j == 0:
                        red_pat = re.compile(r'^([ACGT]{' + str(j) + r'})\1{' + str(int(motiflen/j - 1)) + r'}$')
                        if red_pat.match(motif):
                            redundant = True
                            break
                if redundant:
                    continue
                
                nr += 1
                motif_dict[nr] = motif
                repeats_dict[nr] = len(ssr_string) / motiflen
                end_dict[nr] = m.end()
                start_dict[nr] = m.start() + 1 
                len_dict[nr] = len(seq)
                
                count_motifs[motif] = count_motifs.get(motif, 0) + 1
                if motif not in count_motifs_repeats:
                    count_motifs_repeats[motif] = {}
                rep_round = repeats_dict[nr]
                count_motifs_repeats[motif][rep_round] = count_motifs_repeats[motif].get(rep_round, 0) + 1
                count_class[motiflen] = count_class.get(motiflen, 0) + 1
                
                if max_repeats < rep_round:
                    max_repeats = rep_round

        if nr == 0:
            continue
        
        ssr_containing_seqs[nr] = ssr_containing_seqs.get(nr, 0) + 1
        order = sorted(start_dict.keys(), key=lambda x: start_dict[x])
        
        i = 0
        seq_len = len(seq)

        # 临时存储当前序列所有待写入的潜在记录
        pending_records = []
        # 用于物理位置去重的集合
        written_coordinates = set()

        while i < nr:
            space = amb + 1
            
            # 1. 最后一个或者唯一的 SSR
            if i + 1 >= nr:
                motiflen = len(motif_dict[order[i]])
                if gff:
                    if motiflen == 1:
                        ssrtype, note = "monomeric_repeat", "monomeric_repeat,"
                    else:
                        ssrtype = "microsatellite"
                        note = {2:"dinucleotide_repeat_microsatellite_feature,", 
                                3:"trinucleotide_repeat_microsatellite_feature,", 
                                4:"tetranucleotide_repeat_microsatellite_feature,"}.get(motiflen, "microsatellite,")
                    note += f"({motif_dict[order[i]]}){repeats_dict[order[i]]}"
                else:
                    ssrtype = f"p{motiflen}"
                    ssrseq = f"({motif_dict[order[i]]}){repeats_dict[order[i]]}"
                
                start, end = start_dict[order[i]], end_dict[order[i]]
                pending_records.append({"start": start, "end": end, "ssrtype": ssrtype, "note": note if gff else None, "ssrseq": None if gff else ssrseq})
                i += 1
                continue

            # 2. 与下一个 SSR 的间距大于最大允许间距 (非复合)
            if (start_dict[order[i+1]] - end_dict[order[i]]) > space:
                motiflen = len(motif_dict[order[i]])
                if gff:
                    if motiflen == 1:
                        ssrtype, note = "monomeric_repeat", "monomeric_repeat,"
                    else:
                        ssrtype = "microsatellite"
                        note = {2:"dinucleotide_repeat_microsatellite_feature,", 
                                3:"trinucleotide_repeat_microsatellite_feature,", 
                                4:"tetranucleotide_repeat_microsatellite_feature,"}.get(motiflen, "microsatellite,")
                    note += f"({motif_dict[order[i]]}){repeats_dict[order[i]]}"
                else:
                    ssrtype = f"p{motiflen}"
                    ssrseq = f"({motif_dict[order[i]]}){repeats_dict[order[i]]}"
                
                start, end = start_dict[order[i]], end_dict[order[i]]
                pending_records.append({"start": start, "end": end, "ssrtype": ssrtype, "note": note if gff else None, "ssrseq": None if gff else ssrseq})
                i += 1
                continue

            # 3. 属于复合微卫星 (Compound SSR)
            # 3a. 重叠复合微卫星 (< 1 bp)
            if (start_dict[order[i+1]] - end_dict[order[i]]) < 1:
                ssr_in_compound += 1
                if gff:
                    ssrtype = "repeat_region"
                    m_len1 = len(motif_dict[order[i]])
                    n1 = {1:"monomeric_repeat", 2:"dinucleotide_repeat_microsatellite_feature", 3:"trinucleotide_repeat_microsatellite_feature", 4:"tetranucleotide_repeat_microsatellite_feature"}.get(m_len1, "microsatellite")
                    m_len2 = len(motif_dict[order[i+1]])
                    n2 = {1:"monomeric_repeat", 2:"dinucleotide_repeat_microsatellite_feature", 3:"trinucleotide_repeat_microsatellite_feature", 4:"tetranucleotide_repeat_microsatellite_feature"}.get(m_len2, "microsatellite")
                    note = f"repeat_region,{n1},({motif_dict[order[i]]}){repeats_dict[order[i]]},{n2},({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}"
                    start, end = start_dict[order[i]], end_dict[order[i+1]]
                    ssrseq = None
                else:
                    ssrtype = 'c*'
                    m1, m2 = motif_dict[order[i]], motif_dict[order[i+1]]
                    e1, s2 = end_dict[order[i]], start_dict[order[i+1]]
                    part1 = seq[e1 - len(m1) : e1 - (e1 - s2 + 1)].upper()
                    rep1 = int(repeats_dict[order[i]] - 1)
                    overlap = seq[e1 - (e1 - s2 + 1) : e1].upper()
                    part2 = seq[s2 + (e1 - s2) : s2 + len(m2) - 1].upper()
                    rep2 = int(repeats_dict[order[i+1]] - 1)
                    ssrseq = f"({m1}){rep1}({part1}<{overlap}>{part2})({m2}){rep2}"
                    start, end = start_dict[order[i]], end_dict[order[i+1]]
                    note = None
            
            # 3b. 邻近复合微卫星 (有一定碱基间隔)
            else:
                ssr_in_compound += 1
                e1, s2 = end_dict[order[i]], start_dict[order[i+1]]
                interssr = seq[e1 : s2 - 1].lower()
                
                if gff:
                    m_len1 = len(motif_dict[order[i]])
                    ssrtype = "monomeric_repeat" if m_len1 == 1 else "microsatellite"
                    n1 = {1:"monomeric_repeat", 2:"dinucleotide_repeat_microsatellite_feature", 3:"trinucleotide_repeat_microsatellite_feature", 4:"tetranucleotide_repeat_microsatellite_feature"}.get(m_len1, "microsatellite")
                    note = f"compound_repeat,{n1},({motif_dict[order[i]]}){repeats_dict[order[i]]}"
                    start, end = start_dict[order[i]], end_dict[order[i]]
                    
                    pending_records.append({"start": start, "end": end, "ssrtype": ssrtype, "note": note, "ssrseq": None})
                    
                    m_len2 = len(motif_dict[order[i+1]])
                    ssrtype = "monomeric_repeat" if m_len2 == 1 else "microsatellite"
                    n2 = {1:"monomeric_repeat", 2:"dinucleotide_repeat_microsatellite_feature", 3:"trinucleotide_repeat_microsatellite_feature", 4:"tetranucleotide_repeat_microsatellite_feature"}.get(m_len2, "microsatellite")
                    note = f"compound_repeat,{n2},({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}"
                    start, end = start_dict[order[i+1]], end_dict[order[i+1]]
                    ssrseq = None
                else:
                    ssrseq = f"({motif_dict[order[i]]}){repeats_dict[order[i]]}{interssr}({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}"
                    ssrtype = 'c'
                    start, end = start_dict[order[i]], end_dict[order[i+1]]
                    note = None

            # 级联延伸检查
            while (i + 2 < nr) and ((start_dict[order[i+2]] - end_dict[order[i+1]]) <= space):
                i += 1
                ssr_in_compound += 1
                e_next, s_next = end_dict[order[i]], start_dict[order[i+1]]
                
                if (s_next - e_next) < 1:
                    if gff:
                        ssrtype = "repeat_region"
                        note = note.replace("compound_repeat", "repeat_region") + ","
                        m_len_x = len(motif_dict[order[i+1]])
                        nx = {1:"monomeric_repeat", 2:"dinucleotide_repeat_microsatellite_feature", 3:"trinucleotide_repeat_microsatellite_feature", 4:"tetranucleotide_repeat_microsatellite_feature"}.get(m_len_x, "microsatellite")
                        note += f"{nx},({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}"
                        end = end_dict[order[i+1]]
                    else:
                        ssrseq += f"({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}*"
                        ssrtype = 'c*'
                        end = end_dict[order[i+1]]
                else:
                    interssr_next = seq[e_next : s_next - 1].lower()
                    if gff:
                        pending_records.append({"start": start, "end": end, "ssrtype": ssrtype, "note": note, "ssrseq": None})
                        
                        m_len_x = len(motif_dict[order[i+1]])
                        ssrtype = "monomeric_repeat" if m_len_x == 1 else "microsatellite"
                        nx = {1:"monomeric_repeat", 2:"dinucleotide_repeat_microsatellite_feature", 3:"trinucleotide_repeat_microsatellite_feature", 4:"tetranucleotide_repeat_microsatellite_feature"}.get(m_len_x, "microsatellite")
                        note = f"compound_repeat,{nx},({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}"
                        start, end = start_dict[order[i+1]], end_dict[order[i+1]]
                    else:
                        ssrseq += f"{interssr_next}({motif_dict[order[i+1]]}){repeats_dict[order[i+1]]}"
                        end = end_dict[order[i+1]]
            
            pending_records.append({"start": start, "end": end, "ssrtype": ssrtype, "note": note, "ssrseq": ssrseq})
            i += 1

        # ---- 执行去重、重新排序与物理编号 ----
        final_records = []
        for r in pending_records:
            # 使用 start_end 坐标串组合作为排重键值
            coord_key = f"{r['start']}_{r['end']}"
            if coord_key not in written_coordinates:
                written_coordinates.add(coord_key)
                final_records.append(r)
        
        # 按照起始坐标物理位置升序重排
        final_records.sort(key=lambda x: x["start"])

        # 开始输出写入并严格生成递增编号
        count_seq = 0
        if gff:
            if len(final_records) > 0:
                out_gff = open(f"{seq_id}.gff", "w")
                out_gff.write(f"##gff-version 3\n##sequence-region {seq_id} 1 {seq_len}\n#!Date {loctime}\n#!Type DNA\n#!Source-version MISA {VERSION}\n")
                out_gff.write(f"{seq_id}\tMISA\tregion\t1\t{seq_len}\t.\t.\t.\tID={seq_id}.1\n")
                
                for r in final_records:
                    count_seq += 1
                    # ID 严格顺延：从 2 开始递增 (ID=seq_id.2, ID=seq_id.3 ...)
                    out_gff.write(f"{seq_id}\tMISA\t{r['ssrtype']}\t{r['start']}\t{r['end']}\t.\t.\t.\tNote={r['note']};ID={seq_id}.{count_seq+1}\n")
                out_gff.close()
        else:
            for r in final_records:
                count_seq += 1
                out_misa.write(f"{seq_id}\t{count_seq}\t{r['ssrtype']}\t{r['ssrseq']}\t{r['end'] - r['start'] + 1}\t{r['start']}\t{r['end']}\n")

    if out_misa:
        out_misa.close()

    # 5. 生成统计报表 (.statistics)
    with open(f"{fasta_path}.statistics", "w") as out_stat:
        out_stat.write("Specifications\n==============\n\n")
        out_stat.write(f"Sequence source file: \"{fasta_path}\"\n\n")
        out_stat.write("Definement of microsatellites (unit size / minimum number of repeats):\n")
        for t in typ_sorted:
            out_stat.write(f"({t}/{typrep[t]}) ")
        out_stat.write("\n")
        if amb > 0:
            out_stat.write(f"\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  {amb}\n")
        out_stat.write("\n\n\n")

        total_ssrs = sum(count_class.values())
        ssr_seqs_vals = list(ssr_containing_seqs.values())
        total_ssr_containing_seqs = sum(ssr_seqs_vals)
        more_than_one = total_ssr_containing_seqs - ssr_containing_seqs.get(1, 0)

        out_stat.write("RESULTS OF MICROSATELLITE SEARCH\n================================\n\n")
        out_stat.write(f"Total number of sequences examined:              {number_sequences}\n")
        out_stat.write(f"Total size of examined sequences (bp):           {size_sequences}\n")
        out_stat.write(f"Total number of identified SSRs:                 {total_ssrs}\n")
        out_stat.write(f"Number of SSR containing sequences:              {total_ssr_containing_seqs}\n")
        out_stat.write(f"Number of sequences containing more than 1 SSR:  {more_than_one}\n")
        out_stat.write(f"Number of SSRs present in compound formation:    {ssr_in_compound}\n\n\n")

        out_stat.write("Distribution to different repeat type classes\n---------------------------------------------\n\n")
        out_stat.write("Unit size\tNumber of SSRs\n")
        for t in typ_sorted:
            out_stat.write(f"{t}\t{count_class.get(t, 0)}\n")
        out_stat.write("\n")

        out_stat.write("Frequency of identified SSR motifs\n----------------------------------\n\nRepeats")
        max_rep_int = int(max_repeats)
        for r in range(min_repeats, max_rep_int + 1):
            out_stat.write(f"\t{r}")
        out_stat.write("\ttotal\n")

        sorted_motifs = sorted(count_motifs.keys(), key=lambda x: (len(x), x))
        for motif in sorted_motifs:
            mlen = len(motif)
            out_stat.write(motif)
            for r in range(min_repeats, max_rep_int + 1):
                if r < typrep.get(mlen, 0):
                    out_stat.write("\t-")
                else:
                    val = count_motifs_repeats.get(motif, {}).get(r, "")
                    out_stat.write(f"\t{val}")
            out_stat.write(f"\t{count_motifs[motif]}\n")
        out_stat.write("\n")

        out_stat.write("Frequency of classified repeat types (considering sequence complementary)\n-------------------------------------------------------------------------\n\nRepeats")
        for r in range(min_repeats, max_rep_int + 1):
            out_stat.write(f"\t{r}")
        out_stat.write("\ttotal\n")

        red_rev_groups = {}
        visited_motifs = set()

        for motif in sorted_motifs:
            if motif in visited_motifs:
                continue
            
            group = set()
            cur = motif
            rev_comp = reverse_complement(motif)
            cur_rc = rev_comp
            
            for _ in range(len(motif)):
                group.add(cur)
                group.add(cur_rc)
                cur = cur[1:] + cur[0]
                cur_rc = cur_rc[1:] + cur_rc[0]

            repr_normal = min(group)
            rc_group = {reverse_complement(m) for m in group}
            repr_rc = min(rc_group)

            if repr_normal <= repr_rc:
                group_name = f"{repr_normal}/{repr_rc}"
            else:
                group_name = f"{repr_rc}/{repr_normal}"

            if group_name not in red_rev_groups:
                red_rev_groups[group_name] = {"total": 0, "repeats": {r: 0 for r in range(min_repeats, max_rep_int + 1)}}

            for m in group:
                if m in count_motifs:
                    red_rev_groups[group_name]["total"] += count_motifs[m]
                    for r in range(min_repeats, max_rep_int + 1):
                        red_rev_groups[group_name]["repeats"][r] += count_motifs_repeats.get(m, {}).get(r, 0)
                    visited_motifs.add(m)

        sorted_groups = sorted(red_rev_groups.keys(), key=lambda x: (len(x.split('/')[0]), x))
        for gname in sorted_groups:
            unit_len = len(gname.split('/')[0])
            out_stat.write(gname)
            for r in range(min_repeats, max_rep_int + 1):
                if r < typrep.get(unit_len, 0):
                    out_stat.write("\t-")
                else:
                    v = red_rev_groups[gname]["repeats"][r]
                    out_stat.write(f"\t{v}" if v > 0 else "\t")
            out_stat.write(f"\t{red_rev_groups[gname]['total']}\n")

if __name__ == "__main__":
    main()
