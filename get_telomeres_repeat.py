#!/usr/bin/env python3
"""
Telomere Repeat Analysis Pipeline - v3.5 (With Track Plot & Auto Output Dir)
==============================================================================
Robust detection using:
  1. Exact motif counting (TTTAGGG + CCCTAAA) in sliding windows
  2. 7-bp sequence periodicity analysis (seq[i] == seq[i+7])
  3. Proper reverse-complement handling for both chromosome ends
  4. Fix for circular permutations / phase shifts in variant calling
  5. Context-aware validation for Variant types (2+ repeats OR flanked by Typical)
  6. Automatic 7-panel visualization suite (outputs both PNG and PDF)
"""

import sys
import os
import argparse
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

# 设置 matplotlib 后端为 Agg（非交互式，适用于服务器群集群环境）
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

CANONICAL_FWD = "TTTAGGG"
CANONICAL_REV = "CCCTAAA"

def get_circular_permutations(motif):
    permutations = set()
    for i in range(len(motif)):
        permutations.add(motif[i:] + motif[:i])
    return permutations

CIRCULAR_FWD_SET = get_circular_permutations(CANONICAL_FWD)
CIRCULAR_REV_SET = get_circular_permutations(CANONICAL_REV)

class DualLogger:
    """同时将 stdout 输出到屏幕和指定文件的记录器"""
    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.log = open(filepath, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

def count_motifs_in_window(seq, start, end, motif):
    """Count occurrences of a motif in seq[start:end] (overlapping allowed)."""
    count = 0
    for i in range(start, min(end, len(seq) - len(motif) + 1)):
        if seq[i:i+len(motif)] == motif:
            count += 1
    return count

def find_telomere_boundary(seq, window=500, step=50, density_threshold=0.15):
    """Scan from position 0 inward to find the telomere boundary."""
    scan_len = len(seq)
    n_windows = (scan_len - window) // step

    profile = []
    for wi in range(n_windows):
        start = wi * step
        end = start + window

        n_fwd = count_motifs_in_window(seq, start, end, CANONICAL_FWD)
        n_rev = count_motifs_in_window(seq, start, end, CANONICAL_REV)
        n_total = n_fwd + n_rev

        max_motifs = window / 7.0
        motif_density = n_total / max_motifs if max_motifs > 0 else 0

        periodic_count = 0
        total_positions = 0
        for i in range(start, min(end - 7, scan_len - 7)):
            if seq[i] == seq[i+7]:
                periodic_count += 1
            total_positions += 1
        periodicity = periodic_count / max(1, total_positions)

        combined = motif_density * 0.5 + periodicity * 0.5

        profile.append({
            'position': start + window // 2,
            'motif_density': motif_density,
            'periodicity': periodicity,
            'combined': combined,
            'n_motifs': n_total,
            'motif_type': 'fwd' if n_fwd > n_rev else 'rev',
            'dominant_motif': CANONICAL_FWD if n_fwd > n_rev else CANONICAL_REV,
        })

    boundary = scan_len
    for p in profile:
        if p['combined'] < density_threshold:
            boundary = p['position']
            break

    return boundary, profile

def detailed_telomere_analysis(seq, telo_end, label):
    """
    Detailed characterization of the telomeric region with Context-Aware Variant Filtering.
    """
    telo_seq = seq[:telo_end]

    result = {
        'label': label,
        'telo_length': telo_end,
        'n_exact_fwd': 0,
        'n_exact_rev': 0,
        'n_exact_total': 0,
        'n_positions_7bp_periodic': 0,
        'n_positions_total': 0,
        'variant_motifs': Counter(),
        'all_7mers': Counter(),
        'non_canonical_at_motif_sites': [],
        'canonical_motif_sites': [],
        'raw_seq': seq # 保留原始扫描切片序列用于图表渲染
    }

    # 1. 基础滑动窗口统计所有标准的 7-mers 
    for i in range(len(telo_seq) - 6):
        kmer = telo_seq[i:i+7]
        result['all_7mers'][kmer] += 1
        if kmer == CANONICAL_FWD:
            result['n_exact_fwd'] += 1
        elif kmer == CANONICAL_REV:
            result['n_exact_rev'] += 1

    dominant_motif = CANONICAL_FWD if result['n_exact_fwd'] > result['n_exact_rev'] else CANONICAL_REV

    # 2. 7-bp 周期性统计
    for i in range(len(telo_seq) - 7):
        if telo_seq[i] == telo_seq[i+7]:
            result['n_positions_7bp_periodic'] += 1
        result['n_positions_total'] += 1

    # 3. 第一轮处理：顺序抽取 7-bp 块，标记初步状态
    raw_blocks = []
    # 建立用于寻找最优相位(Phase Shift)的简易统计
    phase_scores = [0] * 7
    for p in range(7):
        for i in range(p, len(telo_seq) - 6, 7):
            if telo_seq[i:i+7] in CIRCULAR_FWD_SET or telo_seq[i:i+7] in CIRCULAR_REV_SET:
                phase_scores[p] += 1
    best_phase = np.argmax(phase_scores) if len(telo_seq) > 0 else 0
    result['best_phase'] = best_phase

    for i in range(best_phase, telo_end - 6, 7):
        variant = telo_seq[i:i+7]
        if variant in CIRCULAR_FWD_SET:
            raw_blocks.append({'pos': i, 'kmer': variant, 'type': 'Typical', 'canonical': CANONICAL_FWD})
        elif variant in CIRCULAR_REV_SET:
            raw_blocks.append({'pos': i, 'kmer': variant, 'type': 'Typical', 'canonical': CANONICAL_REV})
        else:
            raw_blocks.append({'pos': i, 'kmer': variant, 'type': 'Pending_Variant'})

    # 4. 第二轮处理：根据上下文环境（Context）深度过滤变异类型
    num_blocks = len(raw_blocks)
    for idx, block in enumerate(raw_blocks):
        if block['type'] == 'Typical':
            result['canonical_motif_sites'].append((block['pos'], block['canonical']))
            result['n_exact_total'] += 1
        elif block['type'] == 'Pending_Variant':
            current_kmer = block['kmer']
            is_valid_variant = False
            
            has_prev = (idx > 0)
            has_next = (idx < num_blocks - 1)
            prev_block = raw_blocks[idx - 1] if has_prev else None
            next_block = raw_blocks[idx + 1] if has_next else None
            
            is_same_as_prev = (has_prev and prev_block['type'] == 'Pending_Variant' and prev_block['kmer'] == current_kmer)
            is_same_as_next = (has_next and next_block['type'] == 'Pending_Variant' and next_block['kmer'] == current_kmer)
            
            if is_same_as_prev or is_same_as_next:
                is_valid_variant = True
            else:
                flanked_by_typical_prev = (has_prev and prev_block['type'] == 'Typical')
                flanked_by_typical_next = (has_next and next_block['type'] == 'Typical')
                if flanked_by_typical_prev and flanked_by_typical_next:
                    is_valid_variant = True
            
            if is_valid_variant:
                result['variant_motifs'][current_kmer] += 1
                result['non_canonical_at_motif_sites'].append((block['pos'], current_kmer))

    # 计算该末端真正的优势7-mer基序（从Phased序列中提取频数最高的）
    phased_kmers = Counter([telo_seq[i:i+7] for i in range(best_phase, len(telo_seq)-6, 7)])
    result['dominant_unit'] = phased_kmers.most_common(1)[0][0] if phased_kmers else (CANONICAL_REV if label.endswith('right') else CANONICAL_FWD)
    
    # 计算优势基序与标准基序的相似碱基数
    ref_canon = CANONICAL_REV if label.endswith('right') else CANONICAL_FWD
    result['canon_match'] = sum(1 for a, b in zip(result['dominant_unit'], ref_canon))

    return result, dominant_motif


def save_dual_format(fig, out_dir, base_name):
    """辅助函数：一键同时输出高清 PNG 图像与矢量 PDF 图像"""
    png_path = os.path.join(out_dir, f"{base_name}.png")
    pdf_path = os.path.join(out_dir, f"{base_name}.pdf")
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
    print(f"  → 💾 Saved: {base_name}.png & {base_name}.pdf")


def main():
    parser = argparse.ArgumentParser(description='Telomere Repeat Analysis & Visual Suite v3.5')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output file prefix (controls output naming)')
    parser.add_argument('--max_scan', type=int, default=100000, help='Max bases to scan per end')
    parser.add_argument('--window', type=int, default=500, help='Sliding window size for boundary detection')
    parser.add_argument('--threshold', type=float, default=0.15, help='Combined density threshold for boundary')
    args = parser.parse_args()

    # 动态生成图片输出文件夹 (基于 --output 参数)
    fig_dir = f"{args.output}.telomeres_figures"
    os.makedirs(fig_dir, exist_ok=True)

    summary_info_file = f"{args.output}.summary.info"
    sys.stdout = DualLogger(summary_info_file)

    print("=" * 90)
    print("🔬 TELOMERE REPEAT ANALYSIS & VISUALIZATION v3.5 — Dual Output Engine")
    print("=" * 90)
    print(f"\n📁 Input FASTA:  {args.input}")
    print(f"📁 Output Pref:  {args.output}")
    print(f"🎨 Figure Dir:  {fig_dir}/ (PNG + PDF dual formats)")

    print(f"\n[1/3] Loading genome...")
    records = list(SeqIO.parse(args.input, 'fasta'))
    chroms = sorted([rec.id for rec in records])
    print(f"  → {len(records)} chromosomes loaded.")

    print(f"\n[2/3] Detecting telomere arrays & calculating metrics...")
    all_results = []
    # 建立快捷索引结构方便后续图表快速定位查找
    results_dict = defaultdict(dict)
    step = max(1, args.window // 10)

    for rec in records:
        seq_full = str(rec.seq).upper()
        name = rec.id
        chr_len = len(seq_full)

        # Left end (5')
        left_seq = seq_full[:args.max_scan]
        left_boundary, left_profile = find_telomere_boundary(left_seq, window=args.window, step=step, density_threshold=args.threshold)
        left_detail, _ = detailed_telomere_analysis(left_seq, left_boundary, f"{name}_left")
        left_detail.update({'chromosome': name, 'end_label': 'left', 'chr_len': chr_len, 'telo_extended': (left_boundary >= args.max_scan - args.window), 'profile': left_profile})
        all_results.append(left_detail)
        results_dict[name]['left'] = left_detail

        # Right end (3') - reverse complement scanning
        right_seq = seq_full[-args.max_scan:]
        right_seq_rc = str(Seq(right_seq).reverse_complement())
        right_boundary, right_profile = find_telomere_boundary(right_seq_rc, window=args.window, step=step, density_threshold=args.threshold)
        right_detail, _ = detailed_telomere_analysis(right_seq_rc, right_boundary, f"{name}_right")
        right_detail.update({'chromosome': name, 'end_label': 'right', 'chr_len': chr_len, 'telo_extended': (right_boundary >= args.max_scan - args.window), 'profile': right_profile})
        all_results.append(right_detail)
        results_dict[name]['right'] = right_detail

    # 打印控制台表格汇总数据
    print(f"\n[3/3] Analysis completed. Generation of Text Summaries...\n")
    print("=" * 120)
    print("📊 TELOMERE ARRAY SUMMARY")
    print("=" * 120)
    print(f"{'Chr':<8} {'End':<6} {'TeloLen(bp)':<14} {'ExactMotifs':<13} "
          f"{'MotifDens%':<12} {'7bpPeriod%':<12} {'Dominant':<18} "
          f"{'#VariantTypes':<14} {'Status'}")
    print("-" * 120)

    for r in all_results:
        max_motifs = max(1, r['telo_length']) / 7.0
        motif_dens = r['n_exact_total'] / max_motifs * 100
        period_pct = (r['n_positions_7bp_periodic'] / max(1, r['n_positions_total']) * 100)
        status = "⚡>scan_limit" if r['telo_extended'] else "✅bounded"
        n_var = len(r['variant_motifs'])
        print(f"{r['chromosome']:<8} {r['end_label']:<6} {r['telo_length']:<14,} "
              f"{r['n_exact_total']:<13,} {motif_dens:<11.1f}% {period_pct:<11.1f}% "
              f"{r['dominant_unit']:<18} {n_var:<14} {status}")
    print("-" * 120)

    # 文本数据存盘 (TSV & BED)
    tsv_file = f"{args.output}.telomeres_region.tsv"
    with open(tsv_file, 'w') as f:
        f.write("chromosome\tend\ttelo_length_bp\tn_exact_motifs\tmotif_density_pct\tperiodicity_7bp_pct\tdominant_motif\tn_variant_types\textended\n")
        for r in all_results:
            max_m = max(1, r['telo_length']) / 7.0
            f.write(f"{r['chromosome']}\t{r['end_label']}\t{r['telo_length']}\t{r['n_exact_total']}\t{r['n_exact_total']/max_m*100:.1f}\t"
                    f"{r['n_positions_7bp_periodic']/max(1,r['n_positions_total'])*100:.1f}\t{r['dominant_unit']}\t{len(r['variant_motifs'])}\t{r['telo_extended']}\n")

    var_file = f"{args.output}.variants.tsv"
    with open(var_file, 'w') as f:
        f.write("chromosome\tend\tposition\tkmer\tdominant_motif\n")
        for r in all_results:
            for pos, kmer in r['non_canonical_at_motif_sites']:
                f.write(f"{r['chromosome']}\t{r['end_label']}\t{pos}\t{kmer}\t{r['dominant_unit']}\n")

    bed_file = f"{args.output}.telomeres_repeat.bed"
    with open(bed_file, 'w') as f:
        for r in all_results:
            chrom = r['chromosome']
            end_label = r['end_label']
            chr_len = r['chr_len']
            all_repeats_with_gpos = []
            for pos, kmer in r['canonical_motif_sites']:
                g_start = pos if end_label == 'left' else chr_len - pos - 7
                all_repeats_with_gpos.append((g_start, g_start+7, kmer, "Typical", "+" if end_label=='left' else "-"))
            for pos, kmer in r['non_canonical_at_motif_sites']:
                g_start = pos if end_label == 'left' else chr_len - pos - 7
                all_repeats_with_gpos.append((g_start, g_start+7, kmer, "Variant", "+" if end_label=='left' else "-"))
            for g_start, g_end, kmer, r_type, strand in sorted(all_repeats_with_gpos, key=lambda x: x[0]):
                f.write(f"{chrom}\t{g_start}\t{g_end}\t{kmer}_{r_type}\t0\t{strand}\n")

    # ======================================================================
    # 开始生成科学出版级别图形 (Matplotlib & Seaborn 渲染引擎)
    # ======================================================================
    print(f"\n🎨 [4/4] Activating Visualization Suite...")
    sns.set_style("whitegrid")
    
    # 调色板定义
    COLOR_LEFT = '#4393C3'
    COLOR_RIGHT = '#F4A582'
    COLOR_CANON = '#2166AC'
    COLOR_DIVERGED = '#F4A582'
    COLOR_NONCANON = '#B2182B'

    # ══════════════════════════════════════════════════════════════════════
    # 1️⃣ Figure 1: Telomere Length Bar Chart
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 1: Telomere Length Bar Chart...")
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    x_indices = np.arange(len(chroms))
    bar_width = 0.35
    
    left_lens = [results_dict[c]['left']['telo_length'] for c in chroms]
    right_lens = [results_dict[c]['right']['telo_length'] for c in chroms]
    
    bars_l = ax1.bar(x_indices - bar_width/2, left_lens, bar_width, label="Left end (5')", color=COLOR_LEFT, edgecolor='black', linewidth=0.5)
    bars_r = ax1.bar(x_indices + bar_width/2, right_lens, bar_width, label="Right end (3')", color=COLOR_RIGHT, edgecolor='black', linewidth=0.5)
    
    for i, chrom in enumerate(chroms):
        if results_dict[chrom]['left']['canon_match'] <= 4:
            bars_l[i].set_hatch('//')
            bars_l[i].set_edgecolor(COLOR_NONCANON)
        if results_dict[chrom]['right']['canon_match'] <= 4:
            bars_r[i].set_hatch('//')
            bars_r[i].set_edgecolor(COLOR_NONCANON)

    ax1.set_xlabel('Chromosome', fontweight='bold')
    ax1.set_ylabel('Telomere Length (bp)', fontweight='bold')
    ax1.set_title('Figure 1: Telomere Length per Chromosome End\n(Hatched bar = non-canonical / diverged repeats)', fontweight='bold')
    ax1.set_xticks(x_indices)
    ax1.set_xticklabels(chroms, rotation=45)
    ax1.legend()
    plt.tight_layout()
    save_dual_format(fig1, fig_dir, 'fig1_telomere_lengths')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════
    # 2️⃣ Figure 2: Motif Density & Periodicity Profiles (24-panel Matrix)
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 2: Density & Periodicity Profiles...")
    n_plots = len(chroms) * 2
    n_cols = 4
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig2, axes = plt.subplots(n_rows, n_cols, figsize=(20, 3.5 * n_rows))
    axes = axes.flatten()
    
    plot_idx = 0
    for chrom in chroms:
        for end in ['left', 'right']:
            ax = axes[plot_idx]
            r_data = results_dict[chrom][end]
            prof = r_data['profile']
            
            if prof:
                pos_kb = np.array([p['position'] for p in prof]) / 1000.0
                dens = [p['motif_density'] for p in prof]
                per = [p['periodicity'] for p in prof]
                
                ax.plot(pos_kb, dens, color='#377EB8', linewidth=1.5, label='Motif density')
                ax.plot(pos_kb, per, color='#4DAF4A', linewidth=1.2, label='7-bp periodicity')
                ax.axvline(x=r_data['telo_length']/1000.0, color='red', linestyle='--', alpha=0.7)
            
            match_score = r_data['canon_match']
            bg_color = COLOR_CANON if match_score >= 6 else (COLOR_DIVERGED if match_score >= 4 else COLOR_NONCANON)
            
            ax.set_title(f"{chrom}_{end} ({r_data['dominant_unit']})", fontsize=10, fontweight='bold', 
                         bbox=dict(boxstyle='round,pad=0.2', facecolor=bg_color, alpha=0.15))
            ax.set_ylim(-0.05, 1.05)
            ax.tick_params(labelsize=8)
            plot_idx += 1
            
    for i in range(plot_idx, len(axes)):
        fig2.delaxes(axes[i])
        
    fig2.suptitle('Figure 2: Motif Density & 7-bp Periodicity Profiles across Genome Ends', fontweight='bold', fontsize=14, y=1.01)
    plt.tight_layout()
    save_dual_format(fig2, fig_dir, 'fig2_density_profiles')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════
    # 3️⃣ Figure 3: Variant Repeat Heatmap
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 3: Variant Repeat Heatmap...")
    global_variants = Counter()
    for r in all_results:
        global_variants.update(r['variant_motifs'])
        
    top_20_vars = [k for k, _ in global_variants.most_common(20)]
    
    if len(top_20_vars) > 0:
        heatmap_matrix = []
        row_labels = []
        for chrom in chroms:
            for end in ['left', 'right']:
                row_labels.append(f"{chrom}_{end[0].upper()}")
                r_data = results_dict[chrom][end]
                total_v = sum(r_data['variant_motifs'].values())
                row_pct = [ (r_data['variant_motifs'][v] / total_v * 100) if total_v > 0 else 0 for v in top_20_vars ]
                heatmap_matrix.append(row_pct)
                
        fig3, ax3 = plt.subplots(figsize=(14, 8))
        sns.heatmap(heatmap_matrix, xticklabels=top_20_vars, yticklabels=row_labels, cmap='YlOrRd', ax=ax3, cbar_kws={'label': '% of Variant Population'})
        ax3.set_title('Figure 3: Variant Repeat Composition per Chromosome End (Top 20 Variants)', fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_dual_format(fig3, fig_dir, 'fig3_variant_heatmap')
        plt.close()
    else:
        print("  → ⚠️ Skipping Figure 3: No valid variants found.")

    # ══════════════════════════════════════════════════════════════════════
    # 4️⃣ Figure 4: Deep Analysis Multi-Panel
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 4: Deep Analysis Multi-Panel...")
    fig4 = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 3, figure=fig4, hspace=0.4, wspace=0.3)
    
    target_c10 = 'Chr10' if 'Chr10' in results_dict else chroms[0]
    target_c09 = 'Chr09' if 'Chr09' in results_dict else (chroms[1] if len(chroms)>1 else chroms[0])
    target_c06 = 'Chr06' if 'Chr06' in results_dict else chroms[-1]
    
    # A
    ax_a = fig4.add_subplot(gs[0, :])
    ax_a.set_xlim(0, 70); ax_a.set_ylim(-4.5, 0.5); ax_a.axis('off')
    c_map = {'A': '#4CAF50', 'C': '#2196F3', 'G': '#FFC107', 'T': '#F44336'}
    sample_seq = results_dict[target_c10]['left']['raw_seq'][:280]
    for idx, base in enumerate(sample_seq):
        row, col = idx // 70, idx % 70
        ax_a.text(col, -row, base, color=c_map.get(base, 'grey'), fontfamily='monospace', fontweight='bold', fontsize=9, ha='center', va='center')
    ax_a.set_title(f"A. {target_c10} Left End Sequence (First 280 bp)", fontweight='bold', loc='left')
    
    # B
    ax_b = fig4.add_subplot(gs[1, 0])
    dom_u = results_dict[target_c10]['left']['dominant_unit']
    ax_b.bar(np.arange(7)-0.2, [1]*7, 0.4, color=[c_map.get(b, 'grey') for b in CANONICAL_FWD], label=f'Ref: {CANONICAL_FWD}', alpha=0.5, edgecolor='black')
    ax_b.bar(np.arange(7)+0.2, [1]*7, 0.4, color=[c_map.get(b, 'grey') for b in dom_u], label=f'Actual: {dom_u}', edgecolor='black')
    ax_b.set_xticks(range(7)); ax_b.set_xticklabels([f'p{i+1}' for i in range(7)])
    ax_b.set_title(f"B. Base Comparison ({results_dict[target_c10]['left']['canon_match']}/7 Match)")
    ax_b.legend(fontsize=8)
    
    # C
    ax_c = fig4.add_subplot(gs[1, 1])
    p_phase = results_dict[target_c10]['left']['best_phase']
    t_seq = results_dict[target_c10]['left']['raw_seq'][:5000]
    p_kmers = Counter([t_seq[i:i+7] for i in range(p_phase, len(t_seq)-6, 7)])
    top_p_kmers = p_kmers.most_common(5)
    if top_p_kmers:
        ax_c.barh([k[0] for k in top_p_kmers], [k[1] for k in top_p_kmers], color=COLOR_LEFT, edgecolor='black', height=0.6)
    ax_c.set_title("C. Phased 7-mer Top Distribution (First 5 kb)")
    
    # D & E
    for t_chr, t_end, panel_title, gs_loc in [(target_c10, 'left', f'D. {target_c10} Left Profile', gs[1, 2]), (target_c09, 'left', f'E. {target_c09} Left Profile', gs[2, 0])]:
        ax_de = fig4.add_subplot(gs_loc)
        p_data = results_dict[t_chr][t_end]['profile']
        if p_data:
            ax_de.plot([p['position']/1000.0 for p in p_data], [p['motif_density'] for p in p_data], color='#377EB8', label='Density')
            ax_de.plot([p['position']/1000.0 for p in p_data], [p['periodicity'] for p in p_data], color='#4DAF4A', label='Periodicity')
        ax_de.set_title(panel_title); ax_de.set_ylim(-0.05, 1.05); ax_de.legend(fontsize=8)
        
    # F
    ax_f = fig4.add_subplot(gs[2, 1:])
    for tc, te, t_lab in [(target_c10, 'left', 'Target-L'), (target_c09, 'left', 'Canonical-L'), (target_c06, 'right', 'Contrast-R')]:
        curr_seq = results_dict[tc][te]['raw_seq']
        curr_phase = results_dict[tc][te]['best_phase']
        scores = Counter()
        for i in range(curr_phase, len(curr_seq)-6, 7):
            scores[sum(1 for a, b in zip(curr_seq[i:i+7], CANONICAL_REV))] += 1
        tot = sum(scores.values()) if sum(scores.values()) > 0 else 1
        ax_f.plot(range(8), [scores[sc]/tot*100 for sc in range(8)], 'o-', label=f"{tc}_{te} ({t_lab})")
    ax_f.set_xlabel('Match Score (out of 7)'); ax_f.set_ylabel('% Percentage')
    ax_f.set_title('F. Match Score Distribution Comparison'); ax_f.legend(fontsize=8)
    
    fig4.suptitle('Figure 4: Case Deep Dive and Structural Dissection Suite', fontweight='bold', fontsize=14)
    plt.tight_layout()
    save_dual_format(fig4, fig_dir, 'fig4_chr10_deep_analysis')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════
    # 5️⃣ Figure 5: Classification Summary
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 5: Classification Summary Plot...")
    fig5, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    
    c_lens, d_lens, n_lens = [], [], []
    for r in all_results:
        s = r['canon_match']; val_kb = r['telo_length'] / 1000.0
        if s >= 6: c_lens.append(val_kb)
        elif s >= 4: d_lens.append(val_kb)
        else: n_lens.append(val_kb)
        
    sizes = [len(c_lens), len(d_lens), len(n_lens)]
    axes[0].pie(sizes, labels=[f'Canon\nn={sizes[0]}', f'Diverged\nn={sizes[1]}', f'Non-Canon\nn={sizes[2]}'], 
                colors=[COLOR_CANON, COLOR_DIVERGED, COLOR_NONCANON], autopct='%1.1f%%', startangle=90, textprops={'fontweight':'bold'})
    axes[0].set_title('A. Telomere Category Proportions', fontweight='bold')
    
    bp = axes[1].boxplot([c_lens, d_lens, n_lens], patch_artist=True, widths=0.4)
    for patch, color in zip(bp['boxes'], [COLOR_CANON, COLOR_DIVERGED, COLOR_NONCANON]):
        patch.set_facecolor(color); patch.set_alpha(0.5)
    for i, data_g in enumerate([c_lens, d_lens, n_lens]):
        if len(data_g) > 0:
            axes[1].scatter([i+1]*len(data_g) + np.random.normal(0, 0.04, len(data_g)), data_g, color='black', alpha=0.6, s=25, zorder=3)
    axes[1].set_xticklabels(['Canonical', 'Diverged', 'Non-Canon'])
    axes[1].set_ylabel('Length (kb)', fontweight='bold')
    axes[1].set_title('B. Array Length Distribution By Type', fontweight='bold')
    
    for r in all_results:
        prof = r['profile']
        if prof:
            avg_d = np.mean([p['motif_density'] for p in prof[:50]]) 
            avg_p = np.mean([p['periodicity'] for p in prof[:50]])
            s = r['canon_match']
            col = COLOR_CANON if s >= 6 else (COLOR_DIVERGED if s >= 4 else COLOR_NONCANON)
            axes[2].scatter(avg_d, avg_p, color=col, edgecolors='black', s=50, alpha=0.8)
    axes[2].set_xlabel('Mean Motif Density'); axes[2].set_ylabel('Mean 7-bp Periodicity')
    axes[2].set_xlim(-0.05, 1.05); axes[2].set_ylim(-0.05, 1.05)
    axes[2].set_title('C. Density vs Periodicity Clustering Mapping', fontweight='bold')
    
    plt.tight_layout()
    save_dual_format(fig5, fig_dir, 'fig5_classification_summary')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════
    # 6️⃣ Figure 6: Combined Multi-Panel Overview
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 6: Mega-Panel Structural Overview Map...")
    fig6 = plt.figure(figsize=(20, 15))
    gs_main = GridSpec(2, 2, figure=fig6, hspace=0.3, wspace=0.25)
    
    ax_a = fig6.add_subplot(gs_main[0, 0])
    matrix_data = np.zeros((len(chroms), 2)); annot_labels = []
    for i, chrom in enumerate(chroms):
        row_lbls = []
        for j, end in enumerate(['left', 'right']):
            val_kb = results_dict[chrom][end]['telo_length'] / 1000.0
            matrix_data[i, j] = val_kb
            row_lbls.append(f"{val_kb:.1f} kb")
        annot_labels.append(row_lbls)
    sns.heatmap(matrix_data, annot=np.array(annot_labels), fmt="", cmap='YlOrRd', 
                xticklabels=['Left (5\')', 'Right (3\')'], yticklabels=chroms, ax=ax_a, cbar_kws={'label': 'Length (kb)'})
    ax_a.set_title('A. Whole-Genome Telomere Length Matrix Map', fontweight='bold')
    
    ax_b = fig6.add_subplot(gs_main[0, 1])
    bar_y, bar_colors = [], []
    for idx, r in enumerate(all_results):
        bar_y.append(f"{r['chromosome']}_{r['end_label'][0].upper()}")
        s = r['canon_match']
        bar_colors.append(COLOR_CANON if s >= 6 else (COLOR_DIVERGED if s >= 4 else COLOR_NONCANON))
    ax_b.barh(bar_y, [1]*len(bar_y), color=bar_colors, edgecolor='grey', height=0.7)
    ax_b.set_xticks([])
    ax_b.set_title('B. Domain Classification Strip Chart', fontsize=11, fontweight='bold')
    ax_b.tick_params(labelsize=8)
    
    ax_c = fig6.add_subplot(gs_main[1, 0])
    for idx, r in enumerate(all_results):
        prof = r['profile']
        if prof and (r['canon_match'] < 6 or idx % 6 == 0):
            ax_c.plot([p['position']/1000.0 for p in prof], [p['periodicity'] for p in prof], 
                      linewidth=1.5 if r['canon_match']<6 else 0.6, 
                      alpha=0.9 if r['canon_match']<6 else 0.3,
                      color=COLOR_NONCANON if r['canon_match']<=4 else COLOR_CANON)
    ax_c.set_xlim(0, min(40, args.max_scan/1000.0))
    ax_c.set_xlabel('Genomic Distance from Chromosome End (kb)')
    ax_c.set_ylabel('7-bp Periodicity Score')
    ax_c.set_title('C. 7-bp Periodicity Profile Overlay (Red Highlight=Non-Canonical)', fontweight='bold')
    
    ax_d = fig6.add_subplot(gs_main[1, 1])
    ax_d.axis('off')
    lines = ["TELOMERE REPEAT REPERTOIRE CATALOG SUMMARY", "="*58, ""]
    for r in all_results:
        flag = "[C]" if r['canon_match']>=6 else ("[D]" if r['canon_match']>=4 else "[N]")
        lines.append(f" {flag} {r['chromosome']:<6} {r['end_label']:<5} | Dom: {r['dominant_unit']} | Score: {r['canon_match']}/7 | Len: {r['telo_length']:,} bp")
    ax_d.text(0.01, 0.99, "\n".join(lines), transform=ax_d.transAxes, fontfamily='monospace', fontsize=8, va='top', ha='left',
              bbox=dict(boxstyle='round', facecolor='#F8F9FA', edgecolor='#DEE2E6', alpha=0.95))
    ax_d.set_title('D. Detailed Structural Repeat Catalog Box', fontweight='bold')
    
    fig6.suptitle('Figure 6: Comprehensive Chromosome-End Architecture Map', fontweight='bold', fontsize=15)
    plt.tight_layout()
    save_dual_format(fig6, fig_dir, 'fig6_comprehensive_overview')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════
    # 7️⃣ Figure 7: Telomere Repeat Sequence Distribution Tracks (NEW)
    # ══════════════════════════════════════════════════════════════════════
    print("Generating Figure 7: Telomere Repeat Sequence Distribution Tracks...")
    
    # 动态适应画布高度
    fig7, ax7 = plt.subplots(figsize=(16, max(6, len(chroms) * 0.7)))

    # 定义颜色和标签映射字典
    color_map = {}
    label_map = {}
    
    # 1. 典型序列（Typical）颜色
    canon_color = '#2166AC' 
    color_map['CANONICAL'] = canon_color
    label_map[canon_color] = 'Typical (Canonical)'
    
    # 2. 截取频率 Top 5 的变异序列赋予特异性颜色
    top_5_variants = [k for k, _ in global_variants.most_common(5)]
    variant_palette = ['#E69F00', '#009E73', '#CC79A7', '#D55E00', '#56B4E9']
    
    for i, var in enumerate(top_5_variants):
        c = variant_palette[i % len(variant_palette)]
        color_map[var] = c
        label_map[c] = f'Variant: {var}'
        
    # 3. 其余非典型变异颜色
    other_color = '#999999'
    label_map[other_color] = 'Other Variants'

    # 为避免坐标点重绘循环导致渲染缓慢，采用批量数据收集模式
    plot_data = defaultdict(lambda: {'x': [], 'y': []})
    
    y_pos = 0
    y_ticks = []
    y_labels = []

    for chrom in chroms:
        for end in ['left', 'right']:
            r_data = results_dict[chrom][end]
            
            # 画一条浅色的线段代表端粒区域底色背景（代表长度范围）
            telo_len_kb = r_data['telo_length'] / 1000.0
            ax7.plot([0, telo_len_kb], [y_pos, y_pos], color='#EEEEEE', linewidth=4, zorder=1)

            # 收集 Canonical 点位
            for pos, kmer in r_data['canonical_motif_sites']:
                plot_data[canon_color]['x'].append(pos / 1000.0)
                plot_data[canon_color]['y'].append(y_pos)

            # 收集 Variant 点位
            for pos, kmer in r_data['non_canonical_at_motif_sites']:
                c = color_map.get(kmer, other_color)
                plot_data[c]['x'].append(pos / 1000.0)
                plot_data[c]['y'].append(y_pos)

            y_ticks.append(y_pos)
            y_labels.append(f"{chrom}_{end[0].upper()}")
            y_pos -= 1

    # 批量执行 Scatter 散点图绘制出密集条带样式
    for c, coords in plot_data.items():
        if coords['x']:  # 只绘制存在数据的颜色种类
            ax7.scatter(coords['x'], coords['y'], color=c, marker='|', s=80, label=label_map[c], zorder=2)

    # 轴和修饰配置
    ax7.set_yticks(y_ticks)
    ax7.set_yticklabels(y_labels, fontfamily='monospace', fontsize=10)
    ax7.set_xlabel('Distance from Chromosome End (kb)', fontweight='bold', fontsize=12)
    ax7.set_title('Figure 7: Telomere Repeat Sequence Distribution Tracks', fontweight='bold', fontsize=15, pad=15)
    
    # 获取唯一的图例句柄以防止冗余
    handles, labels = ax7.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax7.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=11, frameon=True)

    plt.tight_layout()
    save_dual_format(fig7, fig_dir, 'fig7_telomere_tracks')
    plt.close()

    print("\n[✔] 🚀 Pipeline fully finished! Look inside your directory for data and both PNG/PDF figures!")

if __name__ == "__main__":
    main()
