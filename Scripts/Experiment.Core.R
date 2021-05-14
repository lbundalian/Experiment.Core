# Author : Linnaeus Bundalian


# Guidelines for maintenance ----------------------------------------------
# 
# 1. Naming Convention :
#   a. Functions - lower camel case ( i.e yourMethod)
#   b. Variable - period.separated ( i.e your.variable)
# 2. Use of S4, RefClass or R6 to encapsulate the methods and properties
# 3. Use sections to organize the code
# 4. Matrices are preferable
# 


# Packages ----------------------------------------------------------------
source("Scripts/Packages.R")

# RNASeq Class ------------------------------------------------------------
ExperimentCore <- R6Class("ExperimentCore",
                  public = list(
                    # Properties
                    dge = NULL,
                    data = NULL,
                    assay = NULL,
                    controls = NULL,
                    design = NULL,
                    meta = NULL,
                    stat = NULL,
                    subset = NULL,
                    plots = NULL,
                    fit = NULL,
                    test = NULL,
                    DE = NULL,
                    contrast = NULL,
                    ranked.genes = NULL,
                    GO.Database = NULL,
                    GO.Result = NULL,
                    GO.Enrichment = NULL,
                    GO.Analysis = NULL,
                    GSEA.Database = NULL,
                    GSEA.Result = NULL,
                    GSEA.Enrichment = NULL,
                    GSEA.Analysis = NULL,
                    initialize = function(createDGE = TRUE, meta = NULL, assay = NULL, controls = NULL, 
                                          gene.keys = 'ENSEMBL', gene.ids = c('ENTREZID','SYMBOL')){
                      
                      self$assay <- assay
                      self$meta <- meta
                      self$controls <- controls
                      
                      
                      if (createDGE) {
                        private$initDGE()
                        gene.lookup <- rownames(assay) %>% unlist
                        for (ids in gene.ids) {
                          ids <- casefold(ids,upper=TRUE)
                          #print(ids)
                          self$dge$genes[[ids]] <-  mapIds( org.Hs.eg.db, keytype=gene.keys, 
                                                                  column = ids, keys = gene.lookup )
                          
                        }
                        self$dge <- self$dge[!is.na(self$dge$genes$SYMBOL), ]
                        nrow(self$dge$genes)
                      }
                     
                    },
                    designExperiment = function(remove.lowcounts = TRUE, custom.design = NULL){
                      if (is.null(custom.design)) {
                        self$design <- model.matrix(~0+self$meta$Group)
                        colnames(self$design) <- levels(self$meta$Group)
                        keep <- filterByExpr(self$dge, self$design)
                      } else {
                        
                        self$design <- custom.design

                        
                      }
                      #print(keep)
                      self$dge <- self$dge[keep,, keep.lib.sizes=FALSE]
                      #nrow(self$dge$genes)
                      
                      
                    },
                    
                    # Normalization 
                    normalizeData = function(log2 = TRUE, lib.size = TRUE, prior.count = 1, method = "TMM"){
                      self$dge <-  calcNormFactors(self$dge)
                      self$data <- cpm(self$dge, log=log2, normalized.lib.sizes = lib.size,prior.count = prior.count)
                      self$subset <- self$data
                    },
                    
                    
                    # Statistics
                    calcStatistics = function(q = 20, set = "data"){
                      
                      
                      if (set == "data") {
                        
                        sd <- apply(self$data,1,sd)
                        mean <- apply(self$data,1,mean)
                        median <- apply(self$data,1,median)
                        quant <- apply(self$data,1,quantile,prob=(q/100))
                        var <- apply(self$data,1,var)
                        self$stat <- cbind(sd,mean,median,quant,var)
                        
                        colnames(self$stat) <- c("SD","MEAN","MEDIAN",paste0("Q",q),"VAR")
                        
                      } else if (set == "subset") {
                        
                        sd <- apply(self$subset,1,sd)
                        mean <- apply(self$subset,1,mean)
                        median <- apply(self$subset,1,median)
                        quant <- apply(self$subset,1,quantile,prob=(q/100))
                        var <- apply(self$subset,1,var)
                        self$stat <- cbind(sd,mean,median,quant, var)
                        colnames(self$stat) <- c("SD","MEAN","MEDIAN",paste0("Q",q),"VAR")
                        
                      }
                      
                      
                      
                    },
                    
                    # Dispersion Estimation
                    
                    estimateDisperion = function(robust = TRUE){
                      self$dge <- estimateDisp(self$dge, self$design, robust=robust)
                      self$plots[['BCV']] <- plotBCV(self$dge)
                      self$fit <- glmQLFit(self$dge, self$design, robust=robust)
                      self$plots[['QL']] <- plotQLDisp(self$fit)
                    },
                    
                    contrastVariables = function(variable1 = NULL, variable2 = NULL) {
                      x <- paste0(variable1,'-',variable2)
                      #print(x)
                      self$contrast <- makeContrasts(contrasts = x ,
                                                     levels= self$design)
                    },
                    
                    performDifferential = function(threshold = 1.5){
                      
                      self$test[['QL']] <- glmQLFTest(self$fit, contrast=self$contrast)
                      self$DE[['QL']] <- decideTestsDGE(self$test[['QL']])
                      self$plots[['DE']] <- plotMD(self$test[['QL']], status=self$DE[['QL']])
                      
                      
                      self$test[['LRT']] <- glmTreat(self$fit, contrast=self$contrast, 
                                                    lfc=log2(threshold))
                      
                      self$DE[['LRT']] <- decideTestsDGE(self$test[['LRT']])
                      
                      self$plots[['DE.LRT']] <- plotMD(self$test[['LRT']], status=self$DE[['LRT']])
                      
                    },
                    
                    volcanoPlot = function(test = 'LRT',p.threshold = 0.05, logfc.threshold = 1.5){
                      
                      vp <- ggplot(self$test[[test]]$table %>% mutate(P.adj = p.adjust(.$PValue,method = 'BH'))) + 
                              geom_point(aes(x = logFC, y = -log10(P.adj), color = (P.adj <= -log10(p.threshold) & abs(logFC) >= logfc.threshold )), size = 1) + 
                              geom_vline(xintercept = c(logfc.threshold,-1*logfc.threshold), linetype = 'dashed', color = 'red') +
                              geom_hline(yintercept = c(-log10(p.threshold)), linetype = 'dashed', color = 'red') +
                              ggtitle("Differential Expression") +
                              xlab("log2 Fold Change") + 
                              ylab("-log10 Adj PVal") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(size = rel(1.5), hjust = 0.5), 
                              axis.title = element_text(size = rel(1.25)))
                      i <- length(self$plots)
                      if (is.null(i)) {
                        i <- 1
                      } else {
                        i <- i+1
                      }
                      
                      self$plots[[paste0('VolcanoPlot',i)]] <- vp
                    },
                    
                    createRankedList = function(test = 'LRT',stat = 'logFC'){
                      
                      gene.list <- self$test[[test]]$table %>% dplyr::select(c(stat)) %>% rownames_to_column(.,var="genes") %>%
                                      inner_join(.,self$dge$genes, by= "genes") %>% dplyr::select(c(SYMBOL,stat)) %>%
                                      column_to_rownames("SYMBOL") 
                      
                      gene.stats <- gene.list[[stat]]
                      names(gene.stats) <- rownames(gene.list)
                      self$ranked.genes <- gene.stats
                      
                    },
                    
                    #init GO db 
                    initializeGO = function(test = 'LRT', terms = c('BP','MF','CC'), perm = 10000){
                      
                      self$GO.Enrichment[[test]] <- goana(self$test[[test]] ,geneid = self$dge$genes$ENTREZID
                                             , species = "Hs")
                      self$GO.Enrichment[[test]] <-  self$GO.Enrichment[[test]] %>% rownames_to_column(var = "pathway")
                      go.terms <- terms
                            
                      for (term in go.terms) {
                        #print(term)
                        self$GO.Database[[term]] <- annFUN.org(whichOnto = term,feasibleGenes = NULL,
                                                        mapping = "org.Hs.eg.db",ID = "symbol")
                        #gene.ontodb[[term]] <- select(GO.db, keys=names(gene.onto[[term]]), columns="TERM")
                        self$GO.Result[[term]] <- fgsea(pathways=self$GO.Database[[term]], 
                                                         stats=sort(self$ranked.genes, decreasing = T), nPermSimple = perm)
                      }
                    },
                    
                    runGO = function(test = 'LRT', terms = 'BP', perm = 10000){
                      
                      self$GO.Enrichment[[test]] <- goana(self$test[[test]] ,geneid = self$dge$genes$ENTREZID
                                             , species = "Hs")
                      self$GO.Enrichment[[test]] <-  self$GO.Enrichment[[test]] %>% rownames_to_column(var = "pathway")
                      go.terms <- terms
                            
                      for (term in go.terms) {
                        #print(term)
                        self$GO.Database[[term]] <- annFUN.org(whichOnto = term,feasibleGenes = NULL,
                                                        mapping = "org.Hs.eg.db",ID = "symbol")
                        #gene.ontodb[[term]] <- select(GO.db, keys=names(gene.onto[[term]]), columns="TERM")
                        self$GO.Result[[term]] <- fgsea(pathways=self$GO.Database[[term]], 
                                                         stats=sort(self$ranked.genes, decreasing = T), nPermSimple = perm)
                      }
                    },
                      
                    enrichGO = function(GO.Term = 'BP',n=15, test = 'LRT'){
                                  
                        
                                self$GO.Analysis[[GO.Term]] <- inner_join(self$GO.Result[[GO.Term]],self$GO.Enrichment[[test]],by='pathway')

                                topGoUp <- self$GO.Analysis[[GO.Term]] %>%  filter(padj > 0.05) %>% head(order(pval), n=n)
                                topGoDown <- self$GO.Analysis[[GO.Term]] %>% filter(padj < 0.05) %>% head(order(pval), n=n)

                                topGeneOnto <- rbind(topGoUp, rev(topGoDown))

                                goPlot <- ggplot(topGeneOnto, aes(reorder(Term, NES), NES)) +
                                          geom_col(aes(fill=padj<0.05)) +
                                          coord_flip() +
                                          labs(x="GO Terms", y="Normalized Enrichment Score",
                                               title=paste0("Gene Ontology : ",GO.Term)) + 
                                          theme_minimal()

                                i <- length(self$plots)
                                if (is.null(i)) {
                                    i <- 1
                                } else {
                                    i <- i+1
                                }
                        
                                self$plots[[paste0('GoPlot',i)]] <- goPlot
                                
                                
                                gotop <- self$GO.Analysis[[GO.Term]] %>% filter(padj < 0.05) %>% arrange(desc(size,NES)) %>% head(1)
                                gotop <- self$GO.Analysis[[GO.Term]] %>% filter(padj < 0.05) %>% arrange(desc(NES)) %>% head(1)

                                enrichmentPlot <- plotEnrichment(gotop$leadingEdge  %>% 
                                                  unlist,sort(self$ranked.genes, decreasing = T)) +
                                  xlab("Ranked Genes") + 
                                  ylab("Enrichment Scores") + 
                                  ggtitle(casefold(gotop$Term,upper = TRUE)) 
                                
                                i <- length(self$plots)
                                if (is.null(i)) {
                                    i <- 1
                                } else {
                                    i <- i+1
                                }
                                self$plots[[paste0('Enrichment',i)]] <- enrichmentPlot
  
                    },
                    
                    plotExpression = function(test = 'LRT', n.genes = 5){
                        # prepares data from TMM + log-cpm normalized data to be used in plots
                        norm.counts <- data.frame(self$data) %>% 
                          rownames_to_column(var = "genes") 
                        norm.counts <- gather(norm.counts,key = 'samples', value = 'Norm', 2:7)

                        norm.counts <- inner_join(norm.counts,self$dge$genes,by='genes')
                        norm.counts <- inner_join(norm.counts,self$dge$samples, by = 'samples') %>%
                          inner_join(topTags(self$test[[test]],n.genes)$table,.,by="SYMBOL") %>%
                          dplyr::select(SYMBOL,samples,group,Norm) 

                        i <- length(self$plots)
                        if (is.null(i)) {
                            i <- 1
                        } else {
                            i <- i+1
                        }
                        
                        
                        # Expression plot of Normalized counts of the DE genes (top or desired number of
                        # genes from the DE genes)
                        expressionPlot1 <- ggplot(norm.counts) + geom_point(aes(x=SYMBOL,y=Norm,color = group)) +
                                          scale_y_log10() +
                                          xlab("Genes") + 
                                          ylab("Normalized Counts") + 
                                          ggtitle("DE Genes") +
                                          theme_bw() +
                                          theme(axis.text.x = element_text(angle=45,hjust=1)) +
                                          theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

                        self$plots[[paste0('Expression',i)]] <- expressionPlot1
                        
                        i <- length(self$plots)
                        if (is.null(i)) {
                            i <- 1
                        } else {
                            i <- i+1
                        }
                        
                        # Bar plot of Normalized counts of the DE genes (top or desired number of
                        # genes from the DE genes)
                        expressionPlot2 <- ggplot(norm.counts, aes(x=SYMBOL, y=Norm, fill=group)) + 
                                          geom_bar(stat = "identity",
                                                   position = "dodge") +
                                          xlab("Genes") + 
                                          ylab("Normalized Counts") + 
                                          ggtitle("DE Genes") +
                                          theme_bw() +
                                          theme(axis.text.x = element_text(angle=45,hjust=1)) +
                                          theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
                        self$plots[[paste0('Expression',i)]] <- expressionPlot2
                    },
                    
                    addGOPlot = function(GO.Term = 'BP', test = 'LRT', GO.Pathway = NULL ){
                                
                                
                        
                                if(!is.null(GO.Pathway)){
                                    gotop <- self$GO.Analysis[[GO.Term]] %>% filter(Term == GO.Pathway)
                                } else {
                                    gotop <- self$GO.Analysis[[GO.Term]] %>% filter(padj < 0.05) %>% arrange(desc(size,NES)) %>% head(1)
                                    gotop <- self$GO.Analysis[[GO.Term]] %>% filter(padj < 0.05) %>% arrange(desc(NES)) %>% head(1)
                                }

                                enrichmentPlot <- plotEnrichment(gotop$leadingEdge  %>% 
                                                  unlist,sort(self$ranked.genes, decreasing = T)) +
                                                  xlab("Ranked Genes") + 
                                                  ylab("Enrichment Scores") + 
                                                  ggtitle(casefold(gotop$Term,upper = TRUE)) 
                                
                                i <- length(self$plots)
                                if (is.null(i)) {
                                    i <- 1
                                } else {
                                    i <- i+1
                                }
                        
                                self$plots[[paste0('Enrichment',i)]] <- enrichmentPlot
  
                    },                    
                    
                    initializeGSEA = function(name = "Hallmark",pathway.db = "Data/h.all.v7.4.symbols.gmt" ){
                        self$GSEA.Database[[name]] <- gmtPathways(pathway.db)
                    },
                    
                    performGSEA = function(name = "Hallmark",perm =10000){
                        self$GSEA.Result[[name]] <- fgsea(pathways=self$GSEA.Database[[name]], 
                                                  stats=sort(self$ranked.genes, decreasing = T), nPermSimple = perm)
                        self$GSEA.Result[[name]] <- self$GSEA.Result[[name]] %>% arrange(desc(NES))
                    },
                    
                    plotEnrichedPathways = function(name = "Hallmark",n = 15){
                        
                            gsea.result <- self$GSEA.Result[[name]]

                            topPathwaysUp <- gsea.result[ES > 0][head(order(pval), n = n), pathway]
                            topPathwaysDown <- gsea.result[ES < 0][head(order(pval), n = n), pathway]
                            topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

                            gseaPlot <- ggplot(self$GSEA.Result[[name]] %>% filter(pathway %in% topPathways), aes(reorder(pathway, NES), NES)) +
                                              geom_col(aes(fill=padj<0.05)) +
                                              coord_flip() +
                                              labs(x="Pathway", y="Normalized Enrichment Score",
                                                   title="Gene Set Enrichment Analysis" ) + 
                                              theme_minimal()

                            i <- length(self$plots)
                            if (is.null(i)) {
                                i <- 1
                            } else {
                                i <- i+1
                            }
                            self$plots[[paste0('Enrichment',i)]] <- gseaPlot
                        
                        
                    },
                    
                    plotGSEA = function(database = "Hallmark",pathway = 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', perm = 10000){
                        
                            selected.pathway <- pathway
                            selected.database <- self$GSEA.Database[[database]]
                        
                            self$GSEA.Analysis[[pathway]] <- fgseaSimple(selected.database[selected.pathway], 
                                                                         sort(self$ranked.genes, decreasing = T), nperm=perm)
                            enrichedPlot <- plotEnrichment(selected.database[selected.pathway] %>% 
                                             unlist,sort(self$ranked.genes, decreasing = T)) +
                                                geom_text(aes(label=paste0("pval = ",formatC(self$GSEA.Analysis[[pathway]]$pval,3)," ", 
                                                             "ES = ",formatC(self$GSEA.Analysis[[pathway]]$ES,3)),
                                                x = -Inf, y = Inf),hjust = -2, vjust = 5) +
                                                xlab("Ranked Genes") + 
                                                ylab("Enrichment Scores") + 
                                                ggtitle(selected.pathway)
                            
                            i <- length(self$plots)
                            if (is.null(i)) {
                                i <- 1
                            } else {
                                i <- i+1
                            }
                            self$plots[[paste0('Enrichment',i)]] <- enrichedPlot
                        
                        
                    },
                      
                    # Filtering
                    filterData =  function(sd.threshold = 1, q.threshold = 5,sd.only =FALSE,q.only = FALSE){
                      
                      sub <- self$stat %>% as.data.frame
                      
                      filter.sd <- sub %>% filter(sub[1]>sd.threshold) %>% rownames
                      filter.quantile <- sub %>% filter(sub[4]>q.threshold) %>% rownames
                      
                      
                      if(sd.only){
                        
                        self$subset <- self$subset[filter.sd,]
                        
                      }
                      else if (q.only) {
                        self$subset <- self$subset[filter.quantile,]
                      }
                      else {
                        
                        filter <- c(filter.sd,filter.quantile)
                        self$subset <- self$data[unique(filter),]
                        
                      }
                      
                      
                      
                    }
                    
                  ),
                  private = list(
                    file = NULL,
                    samples = NULL,
                    counts = NULL,
                    genes = NULL,
                    groups = NULL,
                    # Initialize DGE 
                    initDGE = function(...){
                      
                      
                      private$samples <- colnames(self$assay)
                      private$genes <- rownames(self$assay)
                      self$assay <- self$assay %>% as.matrix
                      private$counts <- self$assay
                      private$counts[is.na(private$counts)] <- 0
                      private$groups = self$controls
                      self$dge <- DGEList(counts=private$counts, samples = private$samples, 
                                          genes = private$genes, group = private$groups)
                      
                    }
                    
                  )
                  
                  
                  
)





