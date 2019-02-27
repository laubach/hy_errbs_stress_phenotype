rrbs <- read.table('~/Desktop/terminal_code/output_cort_ewas/cort_ewas_hy_n29.assoc.txt',
                   +                    header = T, row.names = NULL)

p = p.adjust(rrbs$pvalue, method = "BH", length(rrbs$pvalue))

test <- cbind(p, rrbs)