
##### This function takes three files as input as follows:
# 1. sequences is the fasta format of DNA sequences of one or more chromosom fragments
# 2. intervals is a table containing four columns named chr, start, stop, and id. the intervals of the 
# genes to be translated are defined in this file.
# 3. codons is a table containing the genetic codons for RNA translation to amino acids. 

## if no codons file provided the function will go through the R built-in translation function
# and finish the task. 

## The output of the function will be file of standard fasta format of translated amino acids
# for defined intervals, which will be saved in working directory.
# The function will also print the results on the console.

## The required packages are defined in the function, which need to be installed before running the function.


get_ORF_aa = function(sequences, intervals, codons = NULL){
  required_packages <- c("seqinr", "tidyverse", "Biostrings", "foreach", "progress")
  suppressPackageStartupMessages(lapply(required_packages, require, character.only = T))
  
  
  
  sequence <- seqinr::read.fasta(sequences)
  
  interval <- read.delim(intervals)
  
  
  ##### subsetting the genes from chromosomes, transforming to reverse complement to get the coding DNA
  
  
  ## from template
  genes_l_f <- list()
  for(i in 1:nrow(interval)){
    for(j in 1:length(sequence)){
      if(interval[i, "chr"] == names(sequence)[j]){
        genes_l_f[[i]] <- sequence[[j]][interval[i,"start"]:interval[i,"stop"]] %>% c2s()
        names(genes_l_f)[i] <- interval[i, "id"]
      }
    }
  }
  
  
  ## from coding
  genes_l_r <- list()
  for(i in 1:nrow(interval)){
    for(j in 1:length(sequence)){
      if(interval[i, "chr"] == names(sequence)[j]){
        genes_l_r[[i]] <- rev(comp(sequence[[j]][interval[i,"start"]:interval[i,"stop"]])) %>% c2s()
        names(genes_l_r)[i] <- paste0(interval[i, "id"], "_rev", collapse = "")
      }
    }
  }
  
  ## combine both lists of genes in one list
  genes_l <- c(genes_l_r, genes_l_f)
  
  
  
  ##### find the potential start and stop codons
  
  codon_pattern <- c("atg", "taa", "tag", "tga")
  potential_codon <- list()
  genes_codon <- list()
  
  pb <- progress_bar$new(total = length(genes_l))
  for(j in 1:length(genes_l)){
    pb$tick()
    Sys.sleep(1/length(genes_l))
    while(j < (length(genes_l)+1)){
      for(i in 1:4){
        while(i < 5){
          codon_select = codon_pattern[i]
          
          positions = Biostrings::matchPattern(codon_select, genes_l[[j]])@ranges@start
          if(i == 1 & length(positions) == 0){
            print("There is no ATG detected!")
          }
          else if(length(positions) >0){
            potential_codon[[i]] <- positions
            names(potential_codon)[i] <- codon_pattern[i]
          }
          i = i + 1
        }
      }
      genes_codon[[j]] <- potential_codon
      names(genes_codon)[j] <- names(genes_l[j])
      j = j + 1
    }
  }
  
  
  ## create a data frame of start and stop codons for each gene
  
  # Starts
  n <- list()
  Rfs <- data.frame(matrix(ncol = 2)) 
  colnames(Rfs) <- c("Start", "gene")
  
  pb <- progress_bar$new(total = length(genes_codon))
  
  for(j in 1:length(genes_codon)){
    pb$tick()
    Sys.sleep(1/length(genes_codon))
    k = 0
    for(i in 2:length(genes_codon[[j]])){
      k = k + length(genes_codon[[j]][[i]])
      n[[j]] <- k
      row_names <- data.frame(matrix(rep(names(genes_codon[j]), n[[j]]*length(genes_codon[[j]]$atg)), ncol=1))
      colnames(row_names) <- "X"
    }
    
    df <- data.frame(matrix(rep(genes_codon[[j]]$atg, n[[j]]), ncol=1))
    colnames(df) <- "Start"
    df$gene <- row_names$X
    Rfs <- rbind(Rfs, df)
  }
  Rfs <- Rfs[-1,c(2,1)]
  
  # Stops
  Stop <- data.frame(matrix(ncol = 1))
  colnames(Stop) <- "Stop"
  
  pb <- progress_bar$new(total = length(genes_codon))
  for(j in 1:length(genes_codon)){
    pb$tick()
    Sys.sleep(1/length(genes_codon))
    for(i in 2:length(genes_codon[[j]])){
      
      df <- data.frame(matrix(rep(genes_codon[[j]][[i]], each = length(genes_codon[[j]]$atg), ncol=1)))
      colnames(df) <- "Stop"
      Stop <- rbind(Stop, df)
    }
  }
  Stop <- Stop[-1,]
  
  # cbind Stops to Rfs
  Rfs$Stop <- Stop
  
  
  ## select the biologically meaningful ORFs
  Rfs$Size <- Rfs$Stop - Rfs$Start
  
  
  ## filter the reading frames that are happening by chance and the ones that are not in the same frame
  ORFs <- Rfs %>% filter(Size > 100) %>% filter(((Size/3) %% 1) == 0) 
  
  #### getting the seq of ORFs
  ORFs$seq <- NA
  pb <- progress_bar$new(total = length(genes_l))
  for(i in 1:length(genes_l)){
    pb$tick()
    Sys.sleep(1/length(genes_l))
    while(i < (length(genes_l)+1)){
      for(j in 1:nrow(ORFs)){
        while(j < (nrow(ORFs)+1)){
          if(names(genes_l[i]) == ORFs[j,1]){
            ORFs[j,"seq"] <- c2s(s2c(genes_l[[i]])[ORFs[j,"Start"]:ORFs[j,"Stop"]])
          }
          j = j + 1
        }
      }
      i = i + 1
    }
  }
  
  # rename genes in group
  ORFs$gene <- gsub("_rev", "", ORFs$gene)
  
  
  ##### translate ORFs to aa with R built in function
  if(is.null(codons)){ 
  ORFs$aa <- NA
   
  pb <- progress_bar$new(total = nrow(ORFs))
   for(i in 1:nrow(ORFs)){
     
     pb$tick()
     Sys.sleep(1/nrow(ORFs))
     ORFs[i, "aa"] <- c2s(seqinr::translate(s2c(ORFs[i,"seq"])))
   }
   
   ## select the longest orf for each of the genes by ignoring stop codons
   ORFs$n.pattern <- NA
   for(i in 1:nrow(ORFs)){
     ORFs[i,"n.pattern"] <- length(matchPattern("*", ORFs[i,"aa"])@ranges@width)
   }
   r_translated <- ORFs %>% filter(n.pattern == 0) %>% arrange(gene, desc(Size))%>% group_by(gene) %>%
     filter(Size == max(Size))
   
   
   #### write to fasta format
   fasLine <- c()
   for(i in 1:nrow(r_translated)){
     fasLine <- c(fasLine, as.character(paste0(">", r_translated[i, "gene"], collapse = "")))
     fasLine <- c(fasLine, as.character(r_translated[i, "aa"]))
   }
   conn <- file("results.fasta")
   writeLines(fasLine, conn)
   close(conn)
   
   print(fasLine)
  }else{
    codon <- read.csv(codons)
    codon <- separate(codon,col = na..aa,into = c("na", "aa"), sep = "[[:space:]]", remove = T)
  ######### translate with template table
  ORFs$peptide <- NA
  
  pb <- progress_bar$new(total = nrow(ORFs))
  for(i in 1:nrow(ORFs)){
    pb$tick()
    Sys.sleep(1/nrow(ORFs))
    X1 <- NULL
    X2 <- NULL
    for(j in 1:nrow(codon)){
      
      ## define empty vectors to collect the position of aa in sequence
      
      patterns = tolower(codon[j,1])
      positions = matchPattern(patterns, ORFs[i,"seq"]) # find codon in sequence
      if(length(positions) > 0){
        X1 <-base::append(X1,  positions@ranges@start)  # position of codons in sequence
        X2 <- base::append(X2, rep(codon[j,2], length(positions@ranges@start))) # aa of corresponding codon
      }
    }
    result <- data.frame(cbind(X1, X2)) 
    result$X1 <- as.numeric(result$X1)
    result <- result %>% dplyr::arrange(X1) # sort aa by appearance in order
    index <- seq(from = 1, to = nrow(result), by = 3) # select aa in triple frame of codons
    result <- result[index,]
    ORFs[i,"peptide"] <- c2s(result$X2) # write aa in df
  }
  
  
  ## select the longest orf for each of the genes by ignoring stop codons
  ORFs$n.stop <- NA
  for(i in 1:nrow(ORFs)){
    ORFs[i,"n.stop"] <- length(matchPattern("*", ORFs[i,"peptide"])@ranges@width)
  }
  translated <- ORFs %>% filter(n.stop == 0) %>% arrange(gene, desc(Size))%>% group_by(gene) %>%
    filter(Size == max(Size))
  
  
  #### write to fasta format
  fasLine <- c()
  for(i in 1:nrow(translated)){
    fasLine <- c(fasLine, as.character(paste0(">", translated[i, "gene"], collapse = "")))
    fasLine <- c(fasLine, as.character(translated[i, "peptide"]))
  }
  conn <- file("results.fasta")
  writeLines(fasLine, conn)
  close(conn)
  
  
  print(fasLine)
  
  }
}
