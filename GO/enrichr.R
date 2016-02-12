library(httr)
library(readr)
library(openxlsx)
library(plyr)
library(dplyr)
library(stringr)
library(jsonlite)
library(purrr)
library(reshape2)
library(magrittr)
library(igraph)
library(STRINGdb)
library(Cairo)
library(TeachingDemos)

get.enrichrdata <- function(database, gene.df, use.weights = FALSE)
{
    mainurl <- "http://amp.pharm.mssm.edu/Enrichr"
    if (use.weights == TRUE)
    {
        gene.list <- select(gene.df, Symbol, Weight)
        gene.list.combined <- paste(gene.list$Symbol, gene.list$Weight, sep = ",") 
        gene.list.format <- paste(gene.list.combined, collapse = "\n")
    }
    else
    {
        gene.list <- select(gene.df, Symbol) %>% as.matrix %>% as.vector
        gene.list.format <- paste(gene.list, collapse = "\n")
    }

    post.request <- POST(url = paste(mainurl, "addList", sep = "/"), body = list(list = gene.list.format, description = ""))

    #This was changed because Enrichr updated their API so that POST requests to addList return both a shortID and the full length userListID.  The code has been adjusted so that it properly extracts the userListID only, since the uniqueness of the shortIDs is unclear.
    userlist.raw <- content(post.request, "text") 
    userlist <- str_split(userlist.raw, ",")[[1]][2] %>% str_extract("[0-9]+")
    
    get.request <- GET(url = paste(mainurl, "enrich", sep = "/"), query = list(backgroundType = database, userListId = userlist))
    get.content <- content(get.request)[[1]]

    content.df <- lapply(get.content, reshapedata) %>% reduce(rbind) %>% data.frame
    if (ncol(content.df) == 7)
    {
        colnames(content.df) <- c("Index", "GO.Term", "P.value", "Z.score", "Combined.Score", "Genes", "Adj.P.value")
        content.df$P.value %<>% as.numeric
        content.df$Z.score %<>% as.numeric
        content.df$Combined.Score %<>% as.numeric
        content.df$Adj.P.value %<>% as.numeric
        content.df$Index %<>% as.numeric

        content.df$GO.Term %<>% as.character
        content.df$Genes %<>% as.character

        content.df %<>% select(Index:P.value, Adj.P.value, Z.score:Genes)
        return(content.df)
    }    
    else
    {
        print(paste(database, "returned no results"))
        return(NA)
    }
}

reshapedata <- function(orig.list)
{
    orig.list[[6]] %<>% paste(collapse = ",")
    return(orig.list)
}

get.stringdb <- function(symbols.df, plot.name, prefix = "./", edge.threshold = 0)
{
    string_db <- STRINGdb(species = 9606, version = "10")
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
}

