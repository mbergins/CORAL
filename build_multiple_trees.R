library(DarkKinaseTools)
library(here)
library(magick)
library(tidyverse)

source('R/writekinasetree.R')
source('R/colorby.R')

#This data structure contains all the info CORAL uses to build the individual
#kinase trees, to change the appearance of the individual figures, we will
#modify the dataframe variable
base_svg_info = readRDS('Data/kintree.RDS')

base_svg_info$header = "<svg viewBox=\"200 20 640 600\"  preserveAspectRatio=\"xMidYMid meet\"\n\nxmlns=\"http://www.w3.org/2000/svg\"\n\nxmlns:xlink=\"http://www.w3.org/1999/xlink\" >\n"

base_svg_info$dataframe$text.label = base_svg_info$dataframe$id.HGNC
base_svg_info$dataframe$node.radius = 0
base_svg_info$dataframe$text.size = 0;

kinase_svg_info = base_svg_info;
dark_rows = kinase_svg_info$dataframe$id.HGNC %in% dark_kinases$hgnc_symbol
kinase_svg_info$dataframe$node.radius[dark_rows] = 0;
kinase_svg_info$dataframe$text.size[dark_rows] = 6;
kinase_svg_info$dataframe$node.selected[dark_rows] = 1;
kinase_svg_info$dataframe$node.strokewidth[dark_rows] = 0.5;
kinase_svg_info$dataframe$node.opacity[dark_rows] = 0.5;
kinase_svg_info$dataframe$text.col = "#000000"

kinase_svg_info$dataframe$branch.col[dark_rows] = "#008080";
kinase_svg_info$dataframe$branch.group[dark_rows] = this_hgnc;

kinase_svg_info$dataframe = kinase_svg_info$dataframe %>%
 arrange(!is.na(id.HGNC), id.HGNC %in% dark_kinases$hgnc_symbol)

writekinasetree(kinase_svg_info,'dark_tree.svg',font = "Helvetica", labelselect = "HGNC",groupcolor = "#000000") 

dir.create(here('single_kinase_trees'),showWarnings = F)
dir.create(here('single_kinase_trees_png'),showWarnings = F)
svg_data = list()
for (this_hgnc in all_kinases$symbol) {
 
 kinase_row = which(base_svg_info$dataframe$id.HGNC == this_hgnc);
 kinase_svg_info = base_svg_info;
 
 #length of zero means no hits for that kinase name in my list, still ouput a blank tree
 if (length(kinase_row) != 0) {
  
  kinase_svg_info$dataframe$node.radius[kinase_row] = 5;
  kinase_svg_info$dataframe$text.size[kinase_row] = 20;
  kinase_svg_info$dataframe$node.selected[kinase_row] = 1;
  kinase_svg_info$dataframe$node.strokewidth[kinase_row] = 0.5;
  kinase_svg_info$dataframe$node.opacity[kinase_row] = 0.5;
  kinase_svg_info$dataframe$text.col = "#000000"
  
  kinase_svg_info$dataframe$branch.col[kinase_row] = "#008080";
  kinase_svg_info$dataframe$branch.group[kinase_row] = this_hgnc;
  
  kinase_svg_info$dataframe = kinase_svg_info$dataframe %>%
   arrange(!is.na(id.HGNC), id.HGNC == this_hgnc)
 }
 
 writekinasetree(kinase_svg_info,
                 here('single_kinase_trees',paste0(this_hgnc,'.svg')),
                 font = "Helvetica", 
                 labelselect = "HGNC",
                 groupcolor = "#000000")
 svg_data[[this_hgnc]] = image_read_svg(here('single_kinase_trees',paste0(this_hgnc,'.svg'))) %>% 
  image_trim() %>% image_background("white") %>% image_scale(800)
 
 image_write(svg_data[[this_hgnc]], 
             path = here('single_kinase_trees_png/',paste0(this_hgnc,'.png')), 
             format = "png")
}


