

enrichGOSep <- function(gene, output, OrgDb, name="GO",
                         keyType = 'ENTREZID',
                         pval = 0.1, qval = 0.1){
  
  for (ont.i in c("ALL","BP", "CC", "MF")){
    ego <- enrichGO(
                    gene          = gene,
                    keyType       = keyType,
                    OrgDb         = OrgDb,
                    ont           = ont.i,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pval,
                    qvalueCutoff  = qval,
                    readable      = TRUE,
                    pool          = TRUE)
    if (dim(as.data.frame(ego))[1] == 0){
      next()
    }
    print(ont.i)
    p_title <- name
    outputf <- file.path(output, paste(name, ".%s.qval%s_pval%s.pdf", sep = ""))  

    goplot.dir <- file.path(output, "goplot")
    if (!dir.exists(goplot.dir)) {
      dir.create(goplot.dir, recursive = T)
    }

    simgo.dir <- file.path(output, "simgo")
    if (!dir.exists(simgo.dir)) {
      dir.create(simgo.dir, recursive = T)
    }

    out.csv <- file.path(output, paste(name, ".%s.qval%s_pval%s.csv", sep = ""))
    
    write.csv(as.data.frame(ego), file = sprintf(out.csv, ont.i, qval, pval))

    out.pdf <- sprintf(outputf, ont.i, qval, pval)


    
    if (ont.i == "ALL"){
      p <- dotplot(ego, split="ONTOLOGY", showCategory = 15, title = p_title) + 
                facet_grid(ONTOLOGY~., scale="free") + 
                scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) 
    }else {
       p <- dotplot(ego, showCategory=30, title = p_title) + 
                scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) 

       simple.ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

       tryCatch({
        gop <- goplot(simple.ego, max.overlaps= 5, showCategory = 15)
        go.out <- file.path(goplot.dir, paste(name, ".%s_goplot.pdf", sep = ""))
        pdf(sprintf(go.out, ont.i), width = 15, height = 15)
        print(gop)
        dev.off()        
      }, error = function(e) {
         print(e)
      }
      )
      outputsc.pdf <- file.path(simgo.dir, paste(name, "_%s_simple.pdf", sep = ""))
      outputsc.csv <- file.path(simgo.dir, paste(name, "_%s_simple.csv", sep = ""))
      emat  <-  simplifyEnrichment::GO_similarity(ego$ID, ont.i)

      pdf(sprintf(outputsc.pdf, ont.i) ,width = 12, height = 10)
      df = tryCatch({
        df <- simplifyEnrichment::simplifyGO(emat)
      }, warning = function(w) {
          print(w)
      }, error = function(e) {
        file.remove(sprintf(outputsc.pdf, ont.i))
         return(F)
      }
      )

      dev.off()
      write.csv(df, file = sprintf(outputsc.csv, ont.i))
    }
    
    pdf(out.pdf, width = 12, height = 10)
    fit <- try(print(p))
    if ("try-error" %in% class(fit)){
      print(sprintf("%s not enrich", ont.i))
    }
    dev.off()
  }
}


compareClusterSep <- function(gene_ls, output, OrgDb, name="GO",
                         keyType = 'ENSEMBL',
                         pval = 0.1, qval = 0.1){
        for (ont.i in c("BP", "CC", "MF")){
          cpk <- compareCluster(gene_ls, fun="enrichGO", 
                                  pvalueCutoff = pval,
                                  qvalueCutoff = qval,
                                  OrgDb = OrgDb,
                                  ont = ont.i,
                                  keyType = keyType,
                                  readable      = TRUE) 

          if (dim(as.data.frame(cpk))[1] == 0){
            next()
          }
          outputf <- file.path(output, paste(name, ".cmp.%s.qval%s_pval%s.pdf", sep = "")) 
          out.pdf <- sprintf(outputf, ont.i, qval, pval)
          
          p <- dotplot(cpk, title = sprintf("%s_%s",name, ont.i)) + 
                scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                          guide=guide_colorbar(reverse=TRUE, order=1)) +
                guides(size = guide_legend(override.aes=list(shape=1))) +
                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                      panel.grid.major.x = element_blank())

          pdf(out.pdf, width = 12, height = 10)
          fit <- try(print(p))
          if ("try-error" %in% class(fit)){
            print(sprintf("%s not enrich", ont.i))
          }
          dev.off()
  }
}



compareKEGGCluster <- function(gene_ls, output, OrgDb, 
                         organism, name="KEGG",
                         keyType = 'ENSEMBL',
                         pval = 0.1, qval = 0.1){
    if (keyType == 'ENTREZID'){
        gene_ls  <-  gene_ls
    } else {
      gene_ls <- lapply(gene_ls, function(genes){
        genes <- bitr(genes, fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
        genes
      })
       
    }
      cpk <- compareCluster(gene_ls, fun="my_enrichKEGG", 
                              pvalueCutoff = pval,
                              qvalueCutoff = qval,
                              organism = organism,
                             ) 

      if (dim(as.data.frame(cpk))[1] == 0){
        next()
      }
      outputf <- file.path(output, paste(name, ".cmp.qval%s_pval%s.pdf", sep = "")) 
      out.pdf <- sprintf(outputf, qval, pval)
      
      p <- dotplot(cpk, title = sprintf("%s",name)) 
      pdf(out.pdf, width = 12, height = 10)
      fit <- try(print(p))
      if ("try-error" %in% class(fit)){
        print(sprintf("%s not enrich", ont.i))
      }
      dev.off()

}


enrichKEGGSep <- function(gene, out.dir, OrgDb, name = "KEGG",
                        org = "mmu", keytype = 'ENTREZID' ,
                        pval = 0.1, qval = 0.1) {
    if (keytype == 'ENTREZID'){
        gene  <-  gene
    } else {
       gene <- bitr(gene, fromType = keytype, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    }
    
    kk <- my_enrichKEGG(gene     = gene,
                    organism     = org, #hsa
                    pvalueCutoff = pval,
                    qvalueCutoff = qval) #

    if (!dir.exists(out.dir)){
        dir.create(out.dir, recursive = T)
    }

    out.pdf.fmt <- file.path(out.dir, paste(name, ".qval%s_pval%s.pdf", sep = ""))  
    out.csv.fmt <- file.path(out.dir, paste(name, ".qval%s_pval%s.csv", sep = ""))

    out.pdf <- sprintf(out.pdf.fmt, qval, pval)

    tryCatch({
        kp <- dotplot(kk, showCategory=30, title = name)
        pdf(out.pdf, width = 12, height = 10)
        print(kp)
        kk <- setReadable(kk, OrgDb = OrgDb, keyType="ENTREZID")
      }, error = function(e) {
        if (file.exists(out.pdf)){
          file.remove(out.pdf)
        }
        
      }, finally = {
        if (!is.null(dev.list()))
          {
              dev.off()
          }
          
          write.csv(as.data.frame(kk), file = sprintf(out.csv.fmt, qval, pval))
      })
}


my_enrichKEGG <- function(gene, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, 
    pAdjustMethod = "BH", universe, minGSSize = 10, maxGSSize = 500, 
    qvalueCutoff = 0.2, use_internal_data = FALSE) {
    species <- clusterProfiler:::organismMapper(organism)
    if (use_internal_data) {
        KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db(species)
    }
    else {
        kegg_data.file <- paste(organism, lubridate::today(),"kegg","RData", sep = ".")
        if (file.exists(kegg_data.file)){
            print("load local catch")
            load(kegg_data.file)
        }else {
          print("online downloading... ")
          KEGG_DATA <- clusterProfiler:::prepare_KEGG(species, "KEGG", keyType)
          save(KEGG_DATA, file = kegg_data.file) 
        }
        
    }
    res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff, 
        pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
        maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = KEGG_DATA)
    if (is.null(res)) 
        return(res)
    res@ontology <- "KEGG"
    res@organism <- species
    res@keytype <- keyType
    return(res)
}

download_keggdata <- function(organism, ){
    kegg_data.file <- paste(organism, lubridate::today(),"kegg","RData", sep = ".")
    if (file.exists(kegg_data.file)){
        print("load local catch")
        load(kegg_data.file)
    }else {
      print("online downloading... ")
      KEGG_DATA <- clusterProfiler:::prepare_KEGG(species, "KEGG", keyType)
      save(KEGG_DATA, file = kegg_data.file) 
    }
}

load("/home/huangshenghui/project/asla/hsa.2021-09-30.kegg.RData")
load("/home/huangshenghui/project/asla/mmu.2021-09-02.kegg.RData")
load("/home/huangshenghui/project/asla/hsa.2021-10-05.kegg.RData")

my_enrichKEGG

names(KEGG_DATA)

object.size(KEGG_DATA)

object.size(KEGG_DATA$EXTID2PATHID)


my_enrichKEGG()
