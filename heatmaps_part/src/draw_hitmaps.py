import os
import sys
import subprocess
import argparse
import itertools
import pandas as pd
from Bio import SeqIO, AlignIO, SearchIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import seaborn as sns
import matplotlib.pyplot as plt


def drop_duplicates(recs):
    unique_ids = []
    unique_records = []
    for r in recs:
        if r.id not in unique_ids:
            unique_ids.append(r.id)
            unique_records.append(r)
    return tuple(unique_records)


def get_seq_pairs(recs):
    # число возможных пар сочетаний без повторений
    c = tuple(itertools.combinations(recs, 2))
    return c


# def make_muscle_pairwise(pairs):
#     alignments = []
#     for p in pairs:
#         muscle_cline = MuscleCommandline(clwstrict=True)
#         child = subprocess.Popen(str(muscle_cline),
#                                  stdin=subprocess.PIPE,
#                                  stdout=subprocess.PIPE,
#                                  stderr=subprocess.PIPE,
#                                  universal_newlines=True,
#                                  shell=(sys.platform != "win32"))
#         SeqIO.write(p, child.stdin, 'fasta')
#         child.stdin.close()
#         align = AlignIO.read(child.stdout, 'clustal')
#         alignments.append(align)
#     return tuple(alignments)


# def process_muscle_alignmet(alignments, param):
#     result = []
#     for aln in alignments:
#         ident = 0
#         subst = 0
#         gaps = 0
#         pair = [rec for rec in aln]
#         length = aln.get_alignment_length()
#         for let1, let2 in zip(pair[0].seq, pair[1].seq):
#             if let1 == let2:
#                 ident += 1
#             else:
#                 if let1 == '-' or let2 == '-':
#                     gaps += 1
#                 else:
#                     subst += 1
#         if param == 'pident':
#             # процент идентичности 60-100%
#             result.append((pair[0].id, pair[1].id, round(ident / length * 100, 3)))
#         elif param == 'corr_pident':
#             # sqrt(1-pident) 1-0, 0=max_indent
#             result.append((pair[0].id, pair[1].id, round(((1 - (ident / length)) ** 0.5), 3)))
#         elif param == 'subst':
#             # процент замен
#             result.append((pair[0].id, pair[1].id, round(subst  / length * 100, 3)))
#     return tuple(result)


def make_blast_pairwise(pairs):
    alignments = []
    for p in pairs:
        p[0].seq, p[1].seq = p[0].seq.upper(), p[1].seq.upper()
        SeqIO.write(p[0], 'query.fa', 'fasta')
        SeqIO.write(p[1], 'subject.fa', 'fasta')
        blast_cline = NcbiblastnCommandline(cmd='blastn',
                                            query='query.fa',
                                            subject='subject.fa',
                                            out='result.xml',
                                            # dust='no',
                                            # soft_masking='false',
                                            outfmt='5')
        stdout, stderr = blast_cline()
        result = SearchIO.read('result.xml', 'blast-xml')
        print(result)
        alignments.append(result)
        [os.remove(f) for f in ['query.fa', 'subject.fa', 'result.xml']]
    return tuple(alignments)


def process_blast_result(alignments, param):
    data = []
    for aln in alignments:
        hsp = aln[0][0]  # first hit, first hsp
        ident = 0
        subst = 0
        gaps = 0
        pair = [rec for rec in hsp.aln]
        length = hsp.aln.get_alignment_length()
        for let1, let2 in zip(pair[0].seq, pair[1].seq):
            if let1 == let2:
                ident += 1
            else:
                if let1 == '-' or let2 == '-':
                    gaps += 1
                elif let1 == 'N' or let2 == 'N':
                    length -= 1
                else:
                    subst += 1
        if param == 'pident':
            # процент идентичности 60-100%
            data.append((pair[0].id, pair[1].description, round(ident / length * 100, 3)))
        elif param == 'corr_pident':
            # sqrt(1-pident) 1-0, 0=max_indent
            data.append((pair[0].id, pair[1].description, round(((1 - (ident / length)) ** 0.5), 3)))
        elif param == 'subst':
            # процент замен
            data.append((pair[0].id, pair[1].description, round(subst / length * 100, 3)))
    return tuple(data)


def process_mafft_aln(pairs, param):
    result = []
    for pair in pairs:
        ident = 0
        subst = 0
        gaps = 0
        length = len(pair[0])
        for let1, let2 in zip(pair[0].seq, pair[1].seq):
            if let1 == let2:
                ident += 1
            else:
                if let1 == '-' or let2 == '-':
                    gaps += 1
                elif let1 == 'N' or let2 == 'N':
                    length -= 1
                else:
                    subst += 1
        if param == 'pident':
            # процент идентичности 60-100%
            result.append((pair[0].id, pair[1].id, round(ident / length * 100, 3)))
        elif param == 'corr_pident':
            # sqrt(1-pident) 1-0, 0=max_indent
            result.append((pair[0].id, pair[1].id, round(((1 - (ident / length)) ** 0.5), 3)))
        elif param == 'subst':
            # процент замен
            result.append((pair[0].id, pair[1].id, round(subst  / length * 100, 3)))
    return tuple(result)


def write_to_df(results, names, param, prefix, t):
    if param == 'pident':
        df = pd.DataFrame(index=names, columns=names).fillna(100)
    elif param == 'corr_pident' or param == 'subst':
        df = pd.DataFrame(index=names, columns=names).fillna(0)
    for res in results:
        row_name = res[0]
        col_name = res[1]
        val = res[2]
        df.loc[row_name, col_name] = val
        df.loc[col_name, row_name] = val
    # df.loc['Pmult_BNB-2015', 'Pmult_BNB-2015'] = 'ПИСЬКА'
    print(param)
    print(df)
    df.to_csv(f"./results/{t}/csv/{prefix}_{param}.csv")
    return df


def draw_heatmap(df, param, prefix, t):
    mi = min(df.min())
    ma = max(df.max())
    sns.clustermap(df, vmin=mi, vmax=ma, cmap='rocket_r')
    plt.savefig(f"./results/{t}/hm/{prefix}_{param}.png", format='png')


def draw_from_megax(dfs, prefix, t):
    df_R = pd.read_csv(dfs[1], index_col=0)
    df_R.fillna(0, inplace=True)
    df_L = pd.read_csv(dfs[0], index_col=0)
    df_L.fillna(0, inplace=True)
    df = df_R + df_L
    print(df)
    mi = min(df.min())
    ma = max(df.max())
    df.to_csv(f"./results/{t}/csv/{prefix}_megax_dm.csv")
    sns.clustermap(df, vmin=mi, vmax=ma, cmap='rocket_r')
    plt.savefig(f"./results/{t}/hm/{prefix}_megax_dm.png", format='png')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Draw heatmaps')
    parser.add_argument('-s', '--set', type=str, help='Path to sequences set in *.fasta format')
    parser.add_argument('-m', '--megax', nargs=2, help='Path to Left and Right MegaX output in *.csv format')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix for output images and dataframes')
    parser.add_argument('-t', '--type', type=str, help='Type of marker for folder distribution')
    args = parser.parse_args()

    records = drop_duplicates([r for r in SeqIO.parse(args.set, 'fasta')])
    # overwrite deduplicated set
    with open(args.set, 'w') as outf:
        SeqIO.write(records, outf, 'fasta')
    col_row_names = [r.name for r in records]
    df_LR = args.megax

    params = ('pident', 'corr_pident', 'subst')

    combinations = get_seq_pairs(records)
    alignments = make_blast_pairwise(combinations)
    for p in params:
        results = process_blast_result(alignments, param=p)
        # results = process_mafft_aln(pairs=combinations, param=p)
        df = write_to_df(results, col_row_names, param=p, prefix=args.prefix, t=args.type)
        draw_heatmap(df, param=p, prefix=args.prefix, t=args.type)
    draw_from_megax(df_LR, prefix=args.prefix, t=args.type)
