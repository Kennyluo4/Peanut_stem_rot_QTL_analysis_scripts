'''search SNPs within QTL region'''
res = []
for lines in open('27_QTL_flankingMarker.txt'):
    items = lines.replace('\n', '').split('\t')
    QTL, start, end, chr = items[0], items[3], items[5], items[1]
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        continue
    qtl_info = "\t".join(items[0:7])
    if start > end:
        mid = start
        start = end
        end = mid
    for lns in open('combinedSNPs.txt'):
        itms = lns.replace('\n', '').split('\t')
        marker, chromosome, mk_pos = itms[0], itms[1], itms[2]
        try:
            mk_pos = int(mk_pos)
        except ValueError:
            continue
        if chr == chromosome and start <= mk_pos <= end:
            res.append(qtl_info + '\t' + lns)
        else:
            continue
with open('SNP_in_candidate_gene_result.txt', 'w') as handle:
    handle.writelines(res)
