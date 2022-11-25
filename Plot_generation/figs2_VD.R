
library(ggvenn)

removed_taxa <- read.csv('../results/data/Fig_1/No_contamination_removed_taxa.csv')


ggvenn(
    list(
      'microDecon' = ( removed_taxa$microDecon==T ) %>% which(),
      "SCRuB" = ( removed_taxa$SCRuB==T ) %>% which() ,
      "Decontam"= ( removed_taxa$Decontam==T ) %>% which(),
      "Decontam (LB)"= ( removed_taxa$Decontam..LB.==T ) %>% which() 
      ), 
    columns = c("microDecon", 
                "SCRuB" , 
                "Decontam", 
                "Decontam (LB)"),
    fill_color=c('#9467bd', '#c44e52', '#55a868', 'darkgreen'),
    show_percentage = F, 
    fill_alpha = .9
    )

ggsave(  '../results/Supplementary_figures/Fig_S2_a_with_numbers.pdf', device='pdf', dpi=900)



ggvenn(
  list(
    'microDecon' = ( removed_taxa$microDecon==T ) %>% which(),
    "SCRuB" = ( removed_taxa$SCRuB==T ) %>% which() ,
    "Decontam"= ( removed_taxa$Decontam==T ) %>% which(),
    "Decontam (LB)"= ( removed_taxa$Decontam..LB.==T ) %>% which() 
  ), 
  columns = c("microDecon", 
              "SCRuB" , 
              "Decontam", 
              "Decontam (LB)"),
  fill_color=c('#9467bd', '#c44e52', '#55a868', 'darkgreen'),
  show_percentage = F, 
  fill_alpha = .8, 
  show_elements=F, 
  text_size=0, 
  set_name_size=0 ) 

ggsave('../results/Supplementary_figures/Fig_S2_a_without_numbers.pdf', device='pdf', dpi=900)





