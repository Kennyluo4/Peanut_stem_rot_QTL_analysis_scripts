#The QTL file (e.g. "consensusQTLlist.txt") should have 1st column as QTL name, 3rd column as chromosome, 5th column as start, 6th column as end.
#The gene file (e.g. "WMmap.txt") should have 3rd as chromosome, 4th column as gene start, 5th as gene end
result = []

for lines in open('consensusQTLlist.txt'):
      cols = lines.split('\t')
      qtl, chr, left, right = cols[0],cols[2], cols[4], cols[5]
      left = float(left)
      right = float(right)
      qtl_info = lines.replace("\n", "\t")
      for lns in open('WMmap.txt'):
            its = lns.strip().split('\t')
            gene_start, gene_end, chromosome= its[3], its[4],its[2]
            gene_start = int(gene_start)
            gene_end = int(gene_end)
            # if left < gene_start < right and chr == chromosome:
            #       result.append(qtl_info + '\t' + lns)
            # else:
            #       continue
            if left < gene_start < right and chr == chromosome:
                  result.append(qtl_info + '\t' + lns)
            elif left < gene_end < right and chr == chromosome:
                  result.append(qtl_info + '\t' + lns)
            else:
                  continue

with open('B_genome_qtl_gene.txt','w') as handle:
      handle.writelines(result)
