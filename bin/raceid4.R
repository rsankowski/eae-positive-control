#RaceID4 
#create folders
library(tidyverse)
library(viridis)
library(RaceID)
library(Seurat)
library(Matrix)

date = Sys.Date()

micr <- unlist(read_csv("data/microglia-cell-ids.csv")[2])
non_micr <- unlist(read_csv("data/non-microglia-cell-ids.csv")[2])

#add hypothalamic microglia
load("data/prdata.Robj")

prdata <- prdata[, micr]
prdata <- prdata[, !colnames(prdata) %in% non_micr]

#RaceID
sc <- SCseq(prdata)

# filtering of expression data
sc <- filterdata(sc, 
                 mintotal=500,
)

sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('Jun',
                           'Fos',
                           'Zfp36',
                           'Atf3',
                           'Hspa1a|Hspa1b',
                           'Dusp1',
                           'Egr1',
                           'Malat1'))

    sc <- compdist(sc,metric="pearson")
    sc <- clustexp(sc) 
    sc <- findoutliers(sc)
    sc <- comptsne(sc)
    sc <- compfr(sc,knn=10)

#Save sc file
save(sc, file = 'data/datasets/sc-with-mt-genes.Robj')

#analysis without mt genes
prdata = prdata[!grepl("mt-", rownames(prdata)),]


    #RaceID
    sc <- SCseq(prdata)
    
    # filtering of expression data
    sc <- filterdata(sc, 
                     mintotal=500,
    )
    
    sc <- CCcorrect(sc, 
                    dimR = T, 
                    nComp = 20,
                    CGenes = c('Jun',
                               'Fos',
                               'Zfp36',
                               'Atf3',
                               'Hspa1a|Hspa1b',
                               'Dusp1',
                               'Egr1',
                               'Malat1'))
    
    sc <- compdist(sc,metric="pearson")
    sc <- clustexp(sc) 
    sc <- findoutliers(sc)
    sc <- comptsne(sc)
    sc <- compfr(sc,knn=10)
    
    #Save sc file
    save(sc, file = 'data/datasets/sc-without-mt-genes.Robj')
    
plotexpmap(sc,"Mrc1",logsc=F,fr=F)
plotexpmap(sc,"Lyve1",logsc=F,fr=F)
plotexpmap(sc,"Cd163",logsc=F,fr=F)
plotexpmap(sc,"Tmem119",logsc=F,fr=F)
plotexpmap(sc,"Cx3cr1",logsc=F,fr=F)
plotexpmap(sc,"Hexb",logsc=F,fr=F)
plotexpmap(sc,"Ptprc",logsc=F,fr=F)
plotexpmap(sc,"Cd3e",logsc=F,fr=F)
plotexpmap(sc,"Itgam",logsc=F,fr=F)
plotexpmap(sc,"Cd8a",logsc=F,fr=F)
plotexpmap(sc,"Cd4",logsc=F,fr=F)
plotexpmap(sc,"H2-Aa",logsc=F,fr=F)
plotexpmap(sc,"Zbtb46",logsc=F,fr=F)
plotexpmap(sc,"Mog",logsc=F,fr=F)
plotexpmap(sc,"Mbp",logsc=F,fr=F)
plotexpmap(sc,"Gfap",logsc=F,fr=F)
plotexpmap(sc,"Wfdc17",logsc=F,fr=F)
plotexpmap(sc,"Cd79a",logsc=F,fr=F)
plotexpmap(sc,"Cst3",logsc=F,fr=F)
plotexpmap(sc,"Nkg7",logsc=F,fr=F)
plotexpmap(sc,"Rho",logsc=F,fr=F)
plotexpmap(sc,"S100a11",logsc=F,fr=F)
plotexpmap(sc,"Fn1",logsc=F,fr=F)
plotexpmap(sc,"Gfap",logsc=F,fr=F)
plotexpmap(sc,"Cd209a",logsc=F,fr=F)
plotexpmap(sc,"Rho",logsc=F,fr=F)

plotexpmap(sc,"Mrc1",logsc=F,fr=T)
plotexpmap(sc,"Lyve1",logsc=F,fr=T)
plotexpmap(sc,"Cd163",logsc=F,fr=T)
plotexpmap(sc,"Tmem119",logsc=F,fr=T)
plotexpmap(sc,"Cx3cr1",logsc=F,fr=T)
plotexpmap(sc,"Ptprc",logsc=F,fr=T)
plotexpmap(sc,"Cd3e",logsc=F,fr=T)
plotexpmap(sc,"Itgam",logsc=F,fr=T)
plotexpmap(sc,"Cd8a",logsc=F,fr=T)
plotexpmap(sc,"H2-Aa",logsc=F,fr=T)
plotexpmap(sc,"Zbtb46",logsc=F,fr=T)
plotexpmap(sc,"Ly6c2",logsc=F,fr=T)
plotexpmap(sc,"Cd177",logsc=F,fr=T)
plotexpmap(sc,"Igkc",logsc=F,fr=T)
plotexpmap(sc,"Wfdc17",logsc=F,fr=T)
plotexpmap(sc,"Plp1",logsc=F,fr=T)
plotexpmap(sc,"Mog",logsc=F,fr=T)
#plot marker genes


#identify myeloid cells
write.csv(data.frame(names(sc@cpart[sc@cpart %in% c(7,5,4,10,2)])), 'data/microglia-cell-ids.csv')
write.csv(data.frame(names(sc@cpart[sc@cpart %in% c(8,4)])), 'data/non-microglia-cell-ids.csv')
