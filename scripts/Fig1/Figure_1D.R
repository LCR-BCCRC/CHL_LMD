#########################################################
############## Custom functions ##############
##########################################################
create_enrichment_plot (input = For_Forest,
                        SAMPLE = "res_id",
                        Categorical_variables_OI = Categorical_variables_OI,
                        directory_save = ".",
                        min_freq = 5)

######### custom main function:
create_enrichment_plots <- function(input,
                                    SAMPLE = "SAMPLE",
                                    GENE = "GENE",
                                    min_freq = 3,
                                    Categorical_variables_OI,
                                    directory_save = "."){

for(i in Categorical_variables_OI){
   tryCatch({  
  # define parameters:
    print(i)
  filtered_in <- input[!is.na(input[[i]]),]
  # This code is made for binary comparisons, so make sure there are only 2 variables:
  number_of_variables <- length(unique(filtered_in[[i]]))
  if(number_of_variables != 2){
    print("not a binary comparison!")
  }
  ############################################################################################  
  ################### Filter based on mutational frequencies:
  temp_binary <- create_binary_MUT_file(filtered_in,
                                       SAMPLE = SAMPLE,
                                       add_MUT_label = "yes")

  ################################
  indx <- which(colnames(temp_binary) != SAMPLE)
  Freqs_all <- colSums(temp_binary[,indx],na.rm = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(("Feature"))
  colnames(Freqs_all) <- c("Feature","Freq")
  ############################################################################################  
  Features_low_freq <- Freqs_all %>%
    dplyr::filter(Freq < min_freq )
  print(paste("Removing the following features because they are mutated in less than ", min_freq, " samples: ", Features_low_freq$Feature,
              sep = ""))
  filtered_in <- filtered_in %>%
    dplyr::filter(GENE %!in% str_replace(Features_low_freq$Feature,"_MUT",""))
 print(filtered_in)
  
   ############################################################################################  
  ## (Calculating the mutational frequencies)
  filtered_in[[i]][filtered_in[[i]] == "1"] <- "POS"
  filtered_in[[i]][filtered_in[[i]] == "0"] <- "NEG"
  variables_OI <- unique(filtered_in[[i]])
  variables_OI <- variables_OI[order(dplyr::desc(variables_OI))]
  group_1 <- variables_OI[1]
  group_2 <- variables_OI[2]
  size_group_1 <- length(unique(filtered_in[[SAMPLE]][filtered_in[[i]] == group_1]))
  size_group_2 <- length(unique(filtered_in[[SAMPLE]][filtered_in[[i]] == group_2]))
  Number_combined <- Obtain_Mut_FREQs_per_subset(input = filtered_in,
                                                 SAMPLE = SAMPLE, 
                                                 GENE = GENE,
                                                 binary_annotation = i,
                                                 pos_value = group_1,
                                                 neg_value = group_2)
   ############################################################################################  
  dir.create(directory_save)
  name_pdf <- paste(directory_save,
                    "/",
                    i,
                    "_Freq.pdf",
                    sep = "")
  
  Stats_plot <- Calculate_mutation_stats_pdf(Number_combined,
                                  Group_name_1 = group_1,
                                         Group_name_2 = group_2,
                                     p_value_cut_off = 0.05,
                                     p_adj_method = "none",
                                  p_value_label = 0.05,
                                    cut_off_x_label = 1,
                                  name_pdf = name_pdf,
                                  title = "")
},
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

########## custom helper function 1: create binary mutation file
create_binary_MUT_file <- function(input,
                                   SAMPLE = "SAMPLE",
                                   GENE = "GENE",
                                   add_MUT_label = "no"){
  edited_file <- input[,c(SAMPLE,GENE)]
  if(add_MUT_label == "yes"){
  edited_file <- edited_file %>%
    mutate(GENE = paste(GENE,"MUT",sep= "_"))}
  edited_file <- edited_file %>%
    add_column(Binary_column = 1)
  # Remove duplicated rows (multiple hits per gene per sample):
  edited_file <- edited_file[!duplicated(edited_file),]
  edited_file <- edited_file %>%
    spread(GENE, Binary_column)
  edited_file[is.na(edited_file)] <- 0
  return(edited_file)
}
####### custom helper function 2: infer mutational frequencies per subset
Obtain_Mut_FREQs_per_subset <- function(input,
                                        binary_annotation,
                                        SAMPLE = "Tumor_ID",
                                        GENE = "GENE",
                                        pos_value = 1,
                                        neg_value = 0){

  
  ## POS samples:
  if(pos_value != 1){
    input[[binary_annotation]][input[[binary_annotation]] == pos_value] <- "1"
    input[[binary_annotation]][input[[binary_annotation]] == neg_value] <- "0"
    input[[binary_annotation]] <- as.numeric(input[[binary_annotation]])
  }
  indx <- which(input[[binary_annotation]] > 0)
  POS_subset <- input[indx,] %>%
    as.data.frame()
  POS_subset <- POS_subset[,c(SAMPLE,GENE)] 
  POS_subset <- POS_subset %>%
    add_column(Binary_column =1) %>%
    distinct() %>%
    arrange(GENE) 
  POS_subset <- POS_subset %>%
    spread(GENE,Binary_column)
  POS_subset[is.na(POS_subset)] <- 0
  Group_size_POS <- length(unique(POS_subset[[SAMPLE]]))
  print(paste("The positive subset includes ", Group_size_POS, " samples"))
  # Mutational frequencies POS:
  indx <- which(colnames(POS_subset) == SAMPLE)
  Mut_Frequencies_POS <- colSums(POS_subset[,-indx]) %>%
    as.data.frame()
  Mut_Frequencies_POS <- Mut_Frequencies_POS %>%
    rownames_to_column("GENE")
  annotation_POS_MUT <- paste(binary_annotation,"POS_MUT",sep = "_")
  colnames(Mut_Frequencies_POS) <- c("GENE",annotation_POS_MUT)
  ## NEG samples:
  indx <- which(input[[binary_annotation]] == 0)
  NEG_subset <- input[indx,]
  NEG_subset <- NEG_subset[,c(SAMPLE,GENE)] 
  NEG_subset <- NEG_subset %>%
    add_column(Binary_column =1) %>%
    distinct() %>%
    arrange(GENE) 
  NEG_subset <- NEG_subset %>%
    spread(GENE,Binary_column)
  NEG_subset[is.na(NEG_subset)] <- 0
  Group_size_NEG <- length(unique(NEG_subset[[SAMPLE]]))
  print(paste("The negative subset includes ", Group_size_NEG, " samples"))
  # Mutational frequencies NEG:
  indx <- which(colnames(NEG_subset) == SAMPLE)
  Mut_Frequencies_NEG <- colSums(NEG_subset[,-indx]) %>%
    as.data.frame()
  Mut_Frequencies_NEG <- Mut_Frequencies_NEG %>%
    rownames_to_column("GENE")
  annotation_NEG_MUT <- paste(binary_annotation,"NEG_MUT",sep = "_")
  colnames(Mut_Frequencies_NEG) <- c("GENE",annotation_NEG_MUT)
  ## Combining Mutational Frequencies Negative and Positive subsets:
  Mut_Frequencies_combined <- Mut_Frequencies_POS %>%
    full_join(Mut_Frequencies_NEG)
  Mut_Frequencies_combined[is.na(Mut_Frequencies_combined)] <- 0 # Replacing NAs with zeros...
  ## Adding new columns: 
  annotation_POS_WT <- paste(binary_annotation,"POS_WT",sep = "_")
  annotation_NEG_WT <- paste(binary_annotation,"NEG_WT",sep = "_")
  #
  Mut_Frequencies_combined[[annotation_POS_WT]] <- Group_size_POS - Mut_Frequencies_combined[[annotation_POS_MUT]]
  Mut_Frequencies_combined[[annotation_NEG_WT]] <- Group_size_NEG - Mut_Frequencies_combined[[annotation_NEG_MUT]]
  # Put columns into right order for follow-up statistical analysis:
  Mut_Frequencies_combined <- Mut_Frequencies_combined %>%
    subset(select = c("GENE",
                      annotation_POS_MUT,
                      annotation_POS_WT,
                      annotation_NEG_MUT,
                      annotation_NEG_WT))
  return(Mut_Frequencies_combined)
}
####### custom helper function 3: perform Fisher test and plot enrichment
Calculate_mutation_stats_pdf <- function(input_data,
                                         Group_name_1 = NULL,
                                         Group_name_2 = NULL,
                                         p_value_cut_off = 0.01,
                                         p_adj_method = NULL,  # Options are "bf",bonferroni", "holm", "hommel","BH","BY",or "hochberg".
                                     name_pdf = Sys.Date(),
                                     title = NULL,
                                     subtitle = NULL,
                                     additional_subtitle_text = NULL,
                                     xlab = NULL,
                                     ylab = NULL,
                                     bubble_legend_size = FALSE,
                                     bubble_legend_alpha = FALSE,
                                     width = 8,
                                     height = 6,
                                     p_value_label = p_value_cut_off,
                                     cut_off_x_label = 0.4,
                                     cut_off_y_label = cut_off_x_label,
                                     max_lim_x = 1, 
                                     max_lim_y = max_lim_x,
                                     Genes_OI =NULL, 
                                     ...
                                     ) 
{
  colnames(input_data) <- c("Feature","Group_1_POS","Group_1_NEG",
                            "Group_2_POS","Group_2_NEG")

  Mutation_stats <- input_data %>%
    mutate(Group_1_total = Group_1_POS + Group_1_NEG,
           Group_2_total = Group_2_POS + Group_2_NEG,
           Group_1_fraction = Group_1_POS/Group_1_total,
           Group_2_fraction = Group_2_POS/Group_2_total,
           Enrichment = Group_1_fraction/Group_2_fraction)
  # Changing order:  
  Mutation_stats <- Mutation_stats %>%
    dplyr::select("Feature",
                  "Group_1_NEG","Group_1_POS",
                  "Group_2_NEG","Group_2_POS",
                  "Group_1_fraction",
                  "Group_2_fraction",
                  "Enrichment") 
# Performing Fisher's exact test: 
For_Fisher <- Mutation_stats[,1:5] %>%
as.data.frame()
rownames(For_Fisher) <- For_Fisher$Feature
For_Fisher <- For_Fisher[,-1] # removing redundant Feature column...
Fisher_result <- as.data.frame(t(apply(For_Fisher, 1, row_fisher)))
Fisher_result <- Fisher_result %>%
  rownames_to_column("Feature")
Fisher_result <- Fisher_result[,-c(2:5)] # removing columns that are already present in the Mutation_stats table (so I can use "left_join")
# Adding Fisher test results to Mutation_stats table:
Mutation_stats <- Mutation_stats %>%
  left_join(Fisher_result)
Mutation_stats <- Mutation_stats %>%
  mutate(Threshold = case_when(p_val>p_value_cut_off ~ paste(">",p_value_cut_off,sep = " "),
                               p_val<=p_value_cut_off ~ paste("<=",p_value_cut_off,sep = " ")),
         Inverse_P_value = -log10(p_val))
# Re-editing column names:  
colnames(Mutation_stats) <- str_replace(colnames(Mutation_stats),"Group_1",Group_name_1)
colnames(Mutation_stats) <- str_replace(colnames(Mutation_stats),"Group_2",Group_name_2)
# Arrange by p-value: 
Mutation_stats <- Mutation_stats %>%
  arrange(p_val)
# Multiple testing: 
if(!is.null(p_adj_method)) {
  Number_of_tests <- length(unique(Mutation_stats$Feature))
  Mutation_stats <- Mutation_stats %>%
    mutate(p_value_adjusted = p.adjust(p_val, method = p_adj_method, n= Number_of_tests),
           Threshold = case_when(p_value_adjusted>p_value_cut_off ~ paste(">",p_value_cut_off,sep = " "),
                                 p_value_adjusted<=p_value_cut_off ~ paste("<=",p_value_cut_off,sep = " ")))
}
print(Mutation_stats[1:10,])
# Enrichment plot
#
column_OI_1 <- paste(Group_name_1,"fraction",sep = "_")
column_OI_2 <- paste(Group_name_2,"fraction",sep = "_")
# Some editing to adjust custom subtitles:  
if(is.null(xlab)){
  xlab = Group_name_2
}
if(is.null(ylab)){
  ylab = Group_name_1
}
Plot_enrichment <- ggplot(Mutation_stats,
             aes_string(column_OI_2,column_OI_1)) +
        geom_point(shape = 21,
                   size = 10,
                   aes_string(fill = "Threshold",
                              # size = "p_val",
                              alpha = "p_val",
                              )) + # Larger dots for lower p-values
        coord_cartesian(xlim = c(0,max_lim_x),
                        ylim = c(0,max_lim_y)) +
    theme_cowplot() +
        # scale_size_continuous(range = c(.5,15),trans = "reverse") +
        scale_alpha_continuous(range = c(0.01,.85),trans = "reverse") +
        geom_text_repel(aes_string(label = "ifelse(Mutation_stats$Feature %in% Genes_OI|
        p_val < p_value_label|
        Mutation_stats[[column_OI_1]]>cut_off_y_label|
        Mutation_stats[[column_OI_2]]>cut_off_x_label,
                                    as.character(Mutation_stats$Feature),'')"),
                                   segment.colour = "grey",
                                   segment.size = 0.2,
                                   size = 4,
                                   color = "navy",
                                   fontface = "italic",
                                   force = 25,
                        max.overlaps = 60) +
      guides(
            alpha = "none",
             # size = "none",
             fill = guide_legend(title = "p-value",
                                 override.aes = list(size = 7)))  +
                          scale_color_hue() +
            labs(title = title,
                               subtitle = subtitle,
                               x = xlab,
                               y = ylab) +
                           geom_abline(intercept = 0,linetype = 2, color = "grey") +
        scale_fill_manual(values=c("#FF0033","#666666"),na.translate=FALSE) 

print(Plot_enrichment)
name_pdf <- paste(name_pdf,".pdf",sep ="")
name_pdf <- str_replace_all(name_pdf,"-| |:","_")
pdf(name_pdf,width=width,height=height,onefile = F)
print(Plot_enrichment)
dev.off()
}
#### 



