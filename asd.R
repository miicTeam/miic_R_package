library(miic)
setwd('/home/mribeirodantas/dev/miic_webserver/test')
inp <- read.table('simulation_all.tsv', sep='\t', header=T)
res <- miic(inp)
setwd('/home/mribeirodantas/dev/miic_r_package/')
#saveRDS(res, "res_w_newstyle.obj")
#saveRDS(res, "res_no_codingstyle.obj")
saveRDS(res, "res_tip_codingstyle.obj")

