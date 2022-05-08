#!/usr/bin/env python3

## Adapted from https://github.com/BrooksLabUCSC/flair/blob/master/bin/call_diffsplice_events.py

import sys, csv, os

try:
    psl = open(sys.argv[1])
    outfilenamebase = sys.argv[2]
except:
    sys.stderr.write('usage: script.py .psl|.bed out.tsv \n')
    sys.exit(1)


def get_potential_utrs(starts, sizes):
    utrs = []
    for b in range(len(starts)):
        utrs += [(starts[b] + 1, starts[b] + sizes[b])]
    return utrs


def parse_iso_id(iso_gene):
    if '_' not in iso_gene:
        return iso_gene
    if '_chr' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_chr')]
    elif '_XM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XM')]
    elif '_XR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XR')]
    elif '_NM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NM')]
    elif '_NR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NR')]
    elif '_R2_' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_R2_')]
    else:
        iso = iso_gene[:iso_gene.rfind('_')]
    return iso


def update_utr_dict(jdict, exon_start, exon_end):
    if exon_start not in jdict[gene]:
        jdict[gene][exon_start] = {}  # exon_start anchor
    if exon_end not in jdict[gene][exon_start]:
        jdict[gene][exon_start][exon_end] = {}  # exon_end anchor
        jdict[gene][exon_start][exon_end]['isos'] = []  # isoform list for these UTRs
    jdict[gene][exon_start][exon_end]['isos'] += [name]
    return jdict


def find_alt(alljuncs, writer, search_threeutr=True):
    for gene in alljuncs:
        for exon_start in alljuncs[gene]:
            all_tp = list(alljuncs[gene][exon_start].keys())
            n = 0
            strand, chrom, gene_clean = gene.rstrip().split('_')
            if len(all_tp) == 1:  # if there is only one 5' UTR
                continue
            for tp1 in all_tp:
                for tp2 in all_tp:
                    if tp1 == tp2 or tp1 < tp2:
                        # two sites are the same, too close together, or have already been tested in another order
                        continue
                    if strand == '+':
                        inclusion = tp1
                        exclusion = tp2
                    else:
                        inclusion = tp2
                        exclusion = tp1
                    event = chrom + ':' + str(exon_start) + '-' + str(inclusion) + ':' + str(exon_start) + '-' + str(exclusion) + ':' + strand
                    writer.writerow([chrom] + [gene_clean] + [event] + [','.join(alljuncs[gene][exon_start][inclusion]['isos'])] + \
                            [','.join(alljuncs[gene][exon_start][exclusion]['isos'])])
                    n += 1


isoforms = {}
a3_junctions = {}  # alt 3' utr detection
a5_junctions = {}  # alt 5' utr detection
for line in psl:
    line = line.rstrip().split('\t')

    chrom, name, start, end, strand = line[0], line[3], int(line[1]), int(line[2]), line[5]

    blockstarts = [int(n) + start for n in line[11].split(',')[:-1]]
    blocksizes = [int(n) for n in line[10].split(',')[:-1]]

    gene = name.rsplit('_', 1)[1]
    name = parse_iso_id(name)

    gene = strand + '_' + chrom + '_' + gene  # stranded comparisons
    if gene not in isoforms:
        isoforms[gene] = {}
        a3_junctions[gene] = {}
        a5_junctions[gene] = {}


    isoforms[gene][name] = {}
    isoforms[gene][name]['sizes'] = blocksizes
    isoforms[gene][name]['starts'] = blockstarts
    isoforms[gene][name]['range'] = start, end

    these_exons = get_potential_utrs(blockstarts, blocksizes)

    if strand == '+':
        five_UTR_start, five_UTR_end = these_exons[0][0], these_exons[0][1]
        three_UTR_start, three_UTR_end = these_exons[-1][1], these_exons[-1][0]
    else:
        five_UTR_start, five_UTR_end = these_exons[-1][0], these_exons[-1][1]
        three_UTR_start, three_UTR_end = these_exons[0][1], these_exons[0][0]

    print(gene)
    print(five_UTR_start)
    print(three_UTR_start)
    print(five_UTR_end)
    print(three_UTR_end)

    a5_junctions = update_utr_dict(a5_junctions, five_UTR_start, five_UTR_end)
    a3_junctions = update_utr_dict(a3_junctions, three_UTR_start, three_UTR_end)

#print(a5_junctions)
with open(outfilenamebase + '.alt5utr.events.tsv', 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    writer.writerow(['seqnames'] + ['feature_id'] + ['coordinate'] + ['included_isoform'] + ['excluded_isoform'])
    find_alt(a5_junctions, writer, search_threeutr=False)

#print(a3_junctions)
with open(outfilenamebase + '.alt3utr.events.tsv', 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    writer.writerow(['seqnames'] + ['feature_id'] + ['coordinate'] + ['included_isoform'] + ['excluded_isoform'])
    find_alt(a3_junctions, writer, search_threeutr=True)
