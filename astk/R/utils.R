

prepare_KEGG <- function(species, KEGG_Type = "KEGG", keyType = "kegg"){
    astk_dir <- file.path(path.expand('~'), ".astk", "kegg")  
    if (!dir.exists(astk_dir)){
      dir.create(astk_dir, recursive=T)
    }
    kegg_data_name <- paste(species, lubridate::today(),keyType,"RData", sep = ".")
    kegg_data_path <- file.path(astk_dir, kegg_data_name)

    if (file.exists(kegg_data_path)){
        print("load local catch")
        load(kegg_data_path)
    }else {
      print("online downloading... ")
      kegg <- clusterProfiler::download_KEGG(species, KEGG_Type, keyType)
      clusterProfiler:::build_Anno(kegg$KEGGPATHID2EXTID, kegg$KEGGPATHID2NAME)
      save(KEGG_DATA, file = kegg_data_path) 
    }
    return(KEGG_DATA)
}


enrichGOSep <- function(
                        gene, 
                        output, 
                        OrgDb, 
                        name    = "GO",
                        ont     = "BP",
                        keyType = 'ENTREZID',
                        pval    = 0.1, 
                        qval    = 0.1,
                        simple  = F,
                        width   = 12, 
                        height  = 10,
                        format  = "pdf") {
  
    ego <- enrichGO(
                    gene          = gene,
                    keyType       = keyType,
                    OrgDb         = OrgDb,
                    ont           = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pval,
                    qvalueCutoff  = qval,
                    readable      = TRUE,
                    pool          = FALSE)

    outputf <- file.path(output, paste(name, ".%s.qval%s_pval%s.%s", sep = "")) 
    out.pdf <- sprintf(outputf, ont, qval, pval, format)

    if (dim(as.data.frame(ego))[1] == 0){
      print("NO enrichment")
      return()
    }
    print(ont)
    p_title <- ont
     
    out.csv <- file.path(output, paste(name, ".%s.qval%s_pval%s.csv", sep = ""))
    
    write.csv(as.data.frame(ego), file = sprintf(out.csv, ont, qval, pval))

    if (ont == "ALL"){
      p <- dotplot(ego, split="ONTOLOGY", showCategory = 15, title = p_title) + 
                facet_grid(ONTOLOGY~., scale="free") + 
                scale_color_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide   = guide_colorbar(reverse=TRUE, order=1)) +
                guides(size = guide_legend(override.aes=list(shape=1))) +
                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                      panel.grid.major.x = element_blank())                                      
    }else {
       p <- dotplot(ego, showCategory=30, title = p_title) + 
                scale_color_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide   = guide_colorbar(reverse=TRUE, order=1)) +
                guides(size = guide_legend(override.aes=list(shape=1))) +
                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                      panel.grid.major.x = element_blank())                          

      if (simple){
          simple.ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
          simple.ego <- enrichplot::pairwise_termsim(simple.ego)           
          simgo.dir <- file.path(output, "simgo")
          if (!dir.exists(simgo.dir)) {
            dir.create(simgo.dir, recursive = T)
          }            
          outputsc.fig <- file.path(simgo.dir, paste(name, "_%s_simple.", format, sep = ""))
          outputsc.csv <- file.path(simgo.dir, paste(name, "_%s_simple.csv", sep = ""))
          emat  <-  simplifyEnrichment::GO_similarity(ego$ID, ont)
          df = tryCatch({
            df <- simplifyEnrichment::simplifyGO(emat)            
            title = sprintf("%s GO terms clustered by binary_cut", nrow(emat))
            ps <- simplifyEnrichment::ht_clusters(emat, df$cluster, column_title = title)
            
            save_fig(ps, 
                sprintf(outputsc.fig, ont), 
                format = format,
                width  = 12, 
                height = 10, 
                units  = "in")
            write.csv(df, file = sprintf(outputsc.csv, ont)) 

          }, warning = function(w) {
              print(w)
          }, error = function(e) {
            print(e)
            return(F)
          }
          )                 
      }
    }

    save_fig(p, 
        out.pdf, 
        format = format,
        width  = width, 
        height = height, 
        units  = "in")
}


compareClusterSep <- function(gene_ls, 
                              output, 
                              OrgDb, 
                              ont     = "BP",
                              name    = "GO",
                              keyType = 'ENSEMBL',
                              pval    = 0.1, 
                              qval    = 0.1,
                              width   = 12, 
                              height  = 10,
                              format  = "pdf") {

        if (ont == "ALL"){
          onts <- c("BP", "CC", "MF")
        } else {
          onts <- c(ont)
        }
        
        for (ont.i in onts){
          outputf <- file.path(output, paste(name, ".cmp.%s.qval%s_pval%s.%s", sep = "")) 
          out.pdf <- sprintf(outputf, ont.i, qval, pval, format)

          df = tryCatch({
              cpk <- compareCluster(gene_ls, 
                                    fun          = "enrichGO", 
                                    pvalueCutoff = pval,
                                    qvalueCutoff = qval,
                                    OrgDb        = OrgDb,
                                    ont          = ont.i,
                                    keyType      = keyType,
                                    readable     = TRUE) 

              }, warning = function(w) {
                  print(w)
              }, error = function(e) {

                pdf(out.pdf ,width = 12, height = 10)
                dev.off()            
                cpk  <- data.frame()
                return(cpk)
              }
              )
          if (class(df) == "data.frame"){
            cpk <- df
          }

          if (dim(as.data.frame(cpk))[1] == 0){

            pdf(out.pdf ,width = 12, height = 10)
            dev.off()
            next()           
          }

          p <- dotplot(cpk, title = sprintf("%s_%s",name, ont.i)) + 
                scale_color_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide   = guide_colorbar(reverse=TRUE, order=1)) +
                guides(size = guide_legend(override.aes=list(shape=1))) +
                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                      panel.grid.major.x = element_blank())

          save_fig(p, 
              out.pdf, 
              format = format,
              width  = width, 
              height = height, 
              units  = "in")
        }
}


comparePathwayCluster <- function(gene_ls, 
                               output, 
                               OrgDb, 
                               organism, 
                               database,
                               name    = "KEGG",
                               keyType = 'ENSEMBL',
                               pval    = 0.1, 
                               qval    = 0.1,
                               width   = 12, 
                               height  = 10,
                               format  = "pdf") {
    if (keyType == 'ENTREZID'){
        gene_ls  <-  gene_ls
    } else {
      gene_ls <- lapply(gene_ls, function(genes){
        genes <- bitr(genes, fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
        genes
      })
       
    }
      if (database == "KEGG"){
        enrich_func = "enrichKEGG"
      } else if (database == "Reactome"){
        enrich_func = ReactomePA::enrichPathway
      }
      cpk <- compareCluster(gene_ls, 
                            fun          = enrich_func, 
                            pvalueCutoff = pval,
                            qvalueCutoff = qval,
                            organism     = organism)

      if (dim(as.data.frame(cpk))[1] == 0){
        next()
      }
      outputf <- file.path(output, paste(name, ".cmp.qval%s_pval%s.pdf", sep = "")) 
      out.pdf <- sprintf(outputf, qval, pval)
      
      p <- dotplot(cpk, title = sprintf("%s",name))  + 
            scale_color_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                                  guide   = guide_colorbar(reverse=TRUE, order=1)) +
            guides(size = guide_legend(override.aes=list(shape=1))) +
            theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                  panel.grid.major.x = element_blank())
      save_fig(p, 
          out.pdf, 
          format = format,
          width  = width, 
          height = height, 
          units  = "in")
}


enrichKEGGSep <- function(gene, 
                          out.dir, 
                          OrgDb, 
                          name    = "KEGG",
                          org     = "mmu", 
                          keytype = 'ENTREZID' ,
                          pval    = 0.1, 
                          qval    = 0.1,
                          width   = 12, 
                          height  = 10,
                          format  = "pdf") {

    if (keytype == 'ENTREZID'){
        gene  <-  gene
    } else {
       gene <- bitr(gene, fromType = keytype, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    }
    
    kk <- enrichKEGG(gene         = gene,
                     organism     = org, 
                     pvalueCutoff = pval,
                     qvalueCutoff = qval)

    if (!dir.exists(out.dir)){
        dir.create(out.dir, recursive = T)
    }

    if (dim(as.data.frame(kk))[1] == 0){
      print("nothing enrichment")  
      return()
    }                 

    out.fig.fmt <- file.path(out.dir, paste(name, ".qval%s_pval%s.%s", sep = ""))  
    out.csv.fmt <- file.path(out.dir, paste(name, ".qval%s_pval%s.csv", sep = ""))

    out.fig <- sprintf(out.fig.fmt, qval, pval, format)

    tryCatch({
        kp <- dotplot(kk, showCategory=30, title = name) + 
                scale_color_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide   = guide_colorbar(reverse=TRUE, order=1)) +
                guides(size = guide_legend(override.aes=list(shape=1))) +
                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                      panel.grid.major.x = element_blank())
        save_fig(kp, 
            out.fig, 
            format = format,
            width  = width, 
            height = height, 
            units  = "in")
        kk <- setReadable(kk, OrgDb = OrgDb, keyType="ENTREZID")
        write.csv(as.data.frame(kk), file = sprintf(out.csv.fmt, qval, pval))
      }, error = function(e) {
        print(e)
      })
}


enrichReactomeSep <- function(
                          gene, 
                          out.dir, 
                          OrgDb, 
                          name    = "Reactome",
                          org     = "mmu", 
                          keytype = 'ENTREZID' ,
                          pval    = 0.1, 
                          qval    = 0.1,
                          width   = 12, 
                          height  = 10,
                          format  = "pdf") {

    if (keytype == 'ENTREZID'){
        gene  <-  gene
    } else {
       gene <- bitr(gene, fromType = keytype, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    }
    
    kk <- ReactomePA::enrichPathway(
                     gene         = gene,
                     organism     = org, 
                     pvalueCutoff = pval,
                     qvalueCutoff = qval,
                     readable     = FALSE
                     )
    if (!dir.exists(out.dir)){
        dir.create(out.dir, recursive = T)
    }
    if (dim(as.data.frame(kk))[1] == 0){
      print("nothing enrichment")  
      return()
    }                 


    out.fig.fmt <- file.path(out.dir, paste(name, ".qval%s_pval%s.%s", sep = ""))  
    out.csv.fmt <- file.path(out.dir, paste(name, ".qval%s_pval%s.csv", sep = ""))

    out.fig <- sprintf(out.fig.fmt, qval, pval, format)

    tryCatch({
        kp <- dotplot(kk, showCategory=30, title = name) + 
                scale_color_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide   = guide_colorbar(reverse=TRUE, order=1)) +
                guides(size = guide_legend(override.aes=list(shape=1))) +
                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                      panel.grid.major.x = element_blank())
        save_fig(kp, 
            out.fig, 
            format = format,
            width  = width, 
            height = height, 
            units  = "in")
        write.csv(as.data.frame(kk), file = sprintf(out.csv.fmt, qval, pval))
      }, error = function(e) {
        print(e)
      })
}


## from  GOSemSim
load_OrgDb <- function(OrgDb) {
    if (is(OrgDb, "character")) {
        suppressMessages(require(OrgDb, character.only = TRUE))
        OrgDb <- eval(parse(text=OrgDb))
    }
    return(OrgDb)
}    


save_fig <- function(plot, 
                     filename, 
                     format = "png",
                     width  = 6, 
                     height = 6, 
                     units  = "in",
                     res    = 72){
  if (format == "pptx"){
      eoffice::topptx(plot, filename, width = width, height = height)
      cat(sprintf("Saving %s x %s in image\n", width, height))
  } else if (format == "png") {
      png(filename, width = width, height = height, units=units, res = res)
      print(plot)
      dev.off()
      cat(sprintf("Saving %s x %s in image\n", width, height))
  } else if (format == "pdf"){
      pdf(filename, width=width, height=height)
      print(plot)
      dev.off()
      cat(sprintf("Saving %s x %s in image\n", width, height))
  }
}


fig_cmd_parser <- function(){
    parser <- argparse::ArgumentParser()
    parser$add_argument("--output", help="file path")
    parser$add_argument("--fmt", help="figure format")
    parser$add_argument("--resolution", type="integer", help="resolution")
    parser$add_argument("--width", type="double", help="figure width")
    parser$add_argument("--height", type="double", help="figure height")
    return(parser)
}

readFasta <- function(file){
    dna <- Biostrings::readDNAStringSet(file)
    rna <- Biostrings::BStringSet(Biostrings::RNAStringSet(dna))
    return(rna)
}