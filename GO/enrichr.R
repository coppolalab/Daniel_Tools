library(httr)
library(readr)
library(openxlsx)
library(plyr)
library(stringr)
library(jsonlite)
library(purrr)
library(reshape2)
library(magrittr)
library(igraph)
library(STRINGdb)
library(Cairo)
library(TeachingDemos)
library(dplyr)
library(BayesFactor)

enrichr.names <- list()
enrichr.clean <- list()

gobiol <- read_lines("~/Documents/_Research/code/GO/GO_Biological_Process_2015.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
enrichr.names[["GO Biological Process"]] <- map_chr(gobiol, magrittr::extract, 1)
enrichr.clean[["GO Biological Process"]] <- map(gobiol, head, -1) %>% map(tail, -2) %>% map(toupper)

gomole <- read_lines("~/Documents/_Research/code/GO/GO_Molecular_Function_2015.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
enrichr.names[["GO Molecular Function"]] <- map_chr(gomole, magrittr::extract, 1)
enrichr.clean[["GO Molecular Function"]] <- map(gomole, head, -1) %>% map(tail, -2) %>% map(toupper)

gocell <- read_lines("~/Documents/_Research/code/GO/GO_Cellular_Component_2015.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
enrichr.names[["GO Cellular Component"]] <- map_chr(gocell, magrittr::extract, 1)
enrichr.clean[["GO Cellular Component"]] <- map(gocell, head, -1) %>% map(tail, -2) %>% map(toupper)

kegg <- read_lines("~/Documents/_Research/code/GO/KEGG_2016.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
enrichr.names[["KEGG"]] <- map_chr(kegg, magrittr::extract, 1)
enrichr.clean[["KEGG"]] <- map(kegg, head, -1) %>% map(tail, -2) %>% map(str_replace, ",.*$", "") %>% map(toupper)

reactome <- read_lines("~/Documents/_Research/code/GO/Reactome_2016.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
enrichr.names[["Reactome"]] <- map_chr(reactome, magrittr::extract, 1)
enrichr.clean[["Reactome"]] <- map(reactome, head, -1) %>% map(tail, -2) %>% map(str_replace, ",.*$", "") %>% map(toupper)

#gtex.up <- read_lines("~/Documents/_Research/code/GO/GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
#enrichr.names[["GTEx Up"]] <- map_chr(gtex.up, magrittr::extract, 1)
#enrichr.clean[["GTEx Up"]] <- map(gtex.up, head, -1) %>% map(tail, -2) %>% map(str_replace, ",.*$", "") %>% map(toupper)

#gtex.down <- read_lines("~/Documents/_Research/code/GO/GTEx_Tissue_Sample_Gene_Expression_Profiles_down.txt") %>% map(str_split, "\t") %>% map(extract2, 1) 
#enrichr.names[["GTEx Down"]] <- map_chr(gtex.down, magrittr::extract, 1)
#enrichr.clean[["GTEx Down"]] <- map(gtex.down, head, -1) %>% map(tail, -2) %>% map(str_replace, ",.*$", "") %>% map(toupper)

GetHyper <- function(database, gene.sig, all.genes) {
    enrichr.name <- enrichr.names[[database]]
    enrichr.cleaned <- enrichr.clean[[database]]

    sig.intersect <- map(enrichr.cleaned, intersect, gene.sig)
    list.diff <- map_int(sig.intersect, length)

    all.intersect <- map(enrichr.cleaned, intersect, all.genes) 
    list.only <- map_int(all.intersect, length)

    list.notdiff <- list.only - list.diff
    notlist.diff <- length(gene.sig) - list.diff

    notlist.notdiff <- length(all.genes) - length(gene.sig) - notlist.diff + list.diff
    column1 <- map2(list.diff, notlist.diff, c)
    column2 <- map2(list.notdiff, notlist.notdiff, c)
    table.list <- map2(column1, column2, cbind)
    bf.tables <- map(table.list, contingencyTableBF, sampleType = "hypergeom")
    bf.tables.extract <- map(bf.tables, extractBF) %>% map_dbl(extract2, "bf")

    intersect.format <- map(sig.intersect, str_c, collapse = ",") %>% reduce(c)
    bf.df <- data.frame(Term = enrichr.name, Num.Genes = list.diff, Bayes.Factor = bf.tables.extract) %>% filter(Num.Genes > 0)
    bf.df$Genes <- intersect.format
    bf.filter <- filter(bf.df, Bayes.Factor > 3) %>% filter(Num.Genes > 4) %>% arrange(desc(Bayes.Factor))
    bf.filter
}

GetEnrichrData <- function(database, gene.df, use.weights = FALSE) {
    mainurl <- "http://amp.pharm.mssm.edu/Enrichr"
    if (use.weights == TRUE) {
        gene.list <- select(gene.df, Symbol, Weight)
        gene.list.combined <- paste(gene.list$Symbol, gene.list$Weight, sep = ",") 
        gene.list.format <- paste(gene.list.combined, collapse = "\n")
    } else {
        gene.list <- select(gene.df, Symbol) %>% unlist
        gene.list.format <- paste(gene.list, collapse = "\n")
    }

    post.request <- POST(url = paste(mainurl, "addList", sep = "/"), body = list(list = gene.list.format, description = ""))

    #This was changed because Enrichr updated their API so that POST requests to addList return both a shortID and the full length userListID.  The code has been adjusted so that it properly extracts the userListID only, since the uniqueness of the shortIDs is unclear.
    userlist.raw <- content(post.request, "text") 
    userlist <- str_split(userlist.raw, ",")[[1]][2] %>% str_extract("[0-9]+")
    
    get.request <- GET(url = paste(mainurl, "enrich", sep = "/"), query = list(userListId = userlist, backgroundType = database))
    get.content <- content(get.request)[[1]]

    content.df <- map(get.content, ReshapeData) %>% reduce(rbind) %>% data.frame

    if (ncol(content.df) == 9) {
        colnames(content.df) <- c("Index", "Term", "P.value", "Z.score", "Combined.Score", "Genes", "Adj.P.value", "Old.P.value", "Old.Adj.P.value")
        content.df$P.value %<>% as.numeric
        content.df$Z.score %<>% as.numeric
        content.df$Combined.Score %<>% as.numeric
        content.df$Adj.P.value %<>% as.numeric
        content.df$Index %<>% as.numeric
        content.df$Old.P.value %<>% as.numeric
        content.df$Old.Adj.P.value %<>% as.numeric

        content.df$Term %<>% as.character
        content.df$Genes %<>% as.character

        content.df %<>% select(Index:P.value, Adj.P.value, Z.score:Genes)
        return(content.df)
    } else {
        print(paste(database, "returned no results"))
        return(NA)
    }
}

ReshapeData <- function(orig.list) {
    orig.list[[6]] %<>% paste(collapse = ",")
    return(orig.list)
}

GetStringDB <- function(symbols.df, plot.name, prefix = ".", edge.threshold = 0, species.id = 9606) {
    string_db <- STRINGdb(species = species.id, version = "10")
    symbols.mapped <- string_db$map(symbols.df, "Symbol", removeUnmappedRows = TRUE) %>% select(Symbol, STRING_id) #Get STRINGDB protein ids for all symbols
    symbols.filtered <- filter(symbols.mapped, !duplicated(STRING_id)) #Remove redundant proteins
    symbols.subnet <- string_db$get_subnetwork(symbols.filtered$STRING_id) #Get the interaction network for the given protein
    symbols.names <- data.frame(STRING_id = vertex_attr(symbols.subnet, "name")) #Get the exact order of the protein ids
    symbols.names$STRING_id %<>% as.character #Convert protein 
    symbols.joined <- join(symbols.names, symbols.filtered)
    vertex_attr(symbols.subnet, "name") <- symbols.joined$Symbol

    edge.weights <- edge_attr(symbols.subnet, "combined_score")
    pruned.subnet <- delete.edges(symbols.subnet, which(edge.weights < edge.threshold))
    num.edges <- map(1:vcount(pruned.subnet), incident, graph = pruned.subnet) %>% map_dbl(length) 
    pruned.subnet2 <- delete.vertices(pruned.subnet, which(num.edges == 0))
    vertex.colors <- rainbow(vcount(pruned.subnet2))
    V(pruned.subnet2)$color <- vertex.colors
    edge.df <- data.frame(edge_attr(symbols.subnet))
    edge.thickness <- edge.df$combined_score / 200

    filepath <- file.path(prefix, plot.name)

    CairoPDF(filepath, width = 30, height = 30)
    plot.igraph(pruned.subnet2, vertex.size = 2, vertex.label.dist = 0.12, vertex.label.degree = -pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.width = edge.thickness, edge.color = "#0000FF99")
    dev.off()
    return(symbols.subnet)
}

