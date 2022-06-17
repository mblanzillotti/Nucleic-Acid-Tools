#Details========================================================================

#Nucleic Acid Tools for Mass Spectrometry
#Version 15
#MBL 08032021

#This Shiny application is intended for automated fragment identification of
#nucleic acid MS2 data.  Currently, all fragment ion types are considered by
#default. Sequence input uses IDT nomenclature for ease of use, and custom 
#modifications can be indicated immediately after sugar, nucleotide, or base in
#parentheses with the corresponding denotation (s, n, b) e.g.:

    #For a methylated cytosine with a thiophosphate backbone and ribose, the
    #input sequence would be rC(n+14.014)*

#Outputs from this program include IDT sequence parsing, theoretical fragments
#(uncharged and monoisotopic), fragments identified from a deconvoluted mass
#list, identified fragment maps, and residual error plots.

#Please refer any questions to: m.b.lanzillotti@utexas.edu

#Requirements===================================================================

require(tidyverse)
require(plotly)
require(readxl)
require(writexl)
require(shiny)

#Changelog======================================================================

#MBL 04092021
#
#This version is designed for distribution alongside I. Santos' publication. 
#New features added include implementation of IDT-style nomenclature and custom
#modification input, incorporation of tabsets delineating app functions, and 
#plotting of an interactive fragment map and mass errors


#MBL 04192021
#
#This new version added annotated spectrum plotting and more consistent y-bounds
#for fragment maps.  Small changes to map plotting parameters allowed for this, =
#as well as changed to the fragment map and annotated spectrum map plotting functions.
#Furthermore, for better usability for changing between experimental mass lists,
#the "sheet_no" numerical input was moved to the side panel rather than being present
#in each tab.

#MBL 04202021
#
#Added crossring cleavage search checkbox as well as plotting support for these fragments
#in the error plot.  Removed vestigial functions calc_seq_cov() and #comparitor().

#MBL 04222021
#
#Rearranged server side functionality to enable more efficient download handling
#and batch processing.  Now, all deconvoluted mass lists supplied are processed
#at once and accessed one at a time for plots, etc. which enables easy 
#download processing.

#MBL 08032021
#
#Application supplied to C. Crittenden via publication as webpage. Additional features
#added per request include inclusion of sequence coverage in fragment identification
#output, addition of legend to fragment ion map plots, and ability to download .png
#format plots directly from the application.

#MBL 08192021
#
#Per request by C. Crittenden, support for Freestyle or QualBrowser -type headers
#included. FreeStyle headers include less metadata, resulting in fewer rows
#before data. This conflicted with readin_expr_mass_list's skip = 6 argument.
#Adjusted to determine where "Mass" is in first column, and read data from there.

#MBL 09162021
#
#After meeting with C. Crittenden and B. Chen, some additions and bugs were discussed.
#Based on errors appearing to be from deconvolution, +1 & +2 series were added as
#optional searches for any data set, similar to the toggle for cross-ring cleavage.
#Furthermore, a bug was identified when increasing ppm threshold to 100 ppm then
#changing which sheet was being examined.  This appears to be due to a mismatch in
#the scan number selection box and the actual scan being displayed.

#MBL 09282021
#
#Added support for internal ion identification only for unambiguous fragments where
#a terminal ion or its complement was detected. 

#MBL 09292021
#
#Added support for mass lists with no matching fragments, a tibble of NA's is 
#provided from fragment matching functions (terminal and internal) that mitigates
#"subscript out of bounds [[ errors when the rendering functions try to access
#data sets where no fragments are found, which originally resulted in empty tibbles.
#
#Re-worked mass list accession to display experiment name in a selectInput UI
#element rather than a numericInput for more explicit control over which data 
#set to observe. Enabling selection of a specific experiment circumvents the need
#display experiment name above all plots, as it is populated and unambiguous in
#the side bar.

#MBL 09302021
#
#Fixed internal ion filtering batching such that each mass list is properly searched
#without needing to re-search all spectra (it was time consuming). Now works with
#selectInput accession.

#MBL 05252022
#
#Added position arguments to map_params to provide plotting for -1 and -2 ion series
#in fragment map.

#Prerequisites==================================================================

#Required default values for scaffold mass--------------------------------------

sugar_mass <- c(99.0446, 99.0446 + 12.000 + 2*1.007825035, 99.0446 + 12.000)
names(sugar_mass) <- c("r", "m", "+")

nitbase_mass <- c(134.04667, 110.03544, 125.03511, 150.04159, 111.01946)
names(nitbase_mass) <- c("A", "C", "T", "G", "U")

backbone_mass <- c(111.9384)
names(backbone_mass) <- c("*")

#Required parameters for fragment map plot--------------------------------------

map_params <- tribble(
  ~ion_type,	          ~x_adj,	  ~y_adj, 	~shape,	  ~size,   ~color,
  "a",                  0.35,	    0.035,	  15,	      1,       "#228833",
  "a-2'",	              NA,	      NA,	      15,	      1,     	 "#228833",
  "a-2'-3'",      	    NA,	      NA,	      15,	      1,	   	 "#228833",
  "a-3'",	              NA,	      NA,	      15,	      1,	   	 "#228833",
  "a-B",           	    0.35,	    0.08,	    18,	      1,	   	 "#228833",
  "a-B-1'",       	    NA,	      NA,	      15,	      1,	   	 "#228833",
  "a-B-1'-2'-3'-O",	    NA,	      NA,	      15,	      1,	   	 "#228833",
  "a-B-1'-O",	          NA,	      NA,	      15,	      1,	   	 "#228833",
  "a-B+2H",	            0.35,	    0.08,     18,	      1,	   	 "#228833",
  "a-B+H",        	    0.35,	    0.08,     18,	      1,   	 	 "#228833",
  "a+2H",	              0.35,	    0.035,    15,	      1,	   	 "#228833",
  "a+H",          	    0.35,	    0.035,    15,	      1,	   	 "#228833",
  "b",	                0.45,	    0.035,	  15,	      1,	   	 "#4477AA",
  "b+2H",	              0.45,	    0.035,    15,	      1,	   	 "#4477AA",
  "b+H",	              0.45,	    0.035,    15,	      1,     	 "#4477AA",
  "c",	                0.55,   	0.035,	  15,	      1,	   	 "#BB5566",
  "c+2H",          	    0.55,   	0.035,	  15,	      1,	   	 "#BB5566",
  "c+H",           	    0.55,   	0.035,    15,	      1,	   	 "#BB5566",
  "d",	                0.65,	    0.035,	  15,	      1,	   	 "#CCBB44",
  "d+2H",	              0.65,	    0.035,    15,	      1,	   	 "#CCBB44",
  "d+H",          	    0.65,	    0.035,    15,	      1,    	 "#CCBB44",
  "w",	                0.35,	    -0.035,	  15,	      1,	   	 "#228833",
  "w+2H",	              0.35,	    -0.035,   15,	      1,	   	 "#228833",
  "w+H",	              0.35,	    -0.035,   15,	      1,	   	 "#228833",
  "x",	                0.45,	    -0.035,  	15,	      1,	   	 "#4477AA",
  "x-B",           	    0.45,	    -0.08,  	18,	      1,	   	 "#4477AA",
  "x-B+2H",	            0.45,	    -0.08,    18,	      1,    	 "#4477AA",
  "x-B+H",	            0.45,	    -0.08, 	  18,	      1,	   	 "#4477AA",
  "x+2H",	              0.45,	    -0.035,  	15,	      1,	   	 "#4477AA",
  "x+H",          	    0.45,	    -0.035,   15,	      1,	   	 "#4477AA",
  "y",	                0.55,	    -0.035,	  15,	      1,	   	 "#BB5566",
  "y+2H",	              0.55,	    -0.035,	  15,	      1,	   	 "#BB5566",
  "y+H",	              0.55,	    -0.035,	  15,	      1,	   	 "#BB5566",
  "z",	                0.65,	    -0.035,  	15,	      1,   	   "#CCBB44",
  "z-5'",	              NA,	      NA,	      15,	      1,    	 "#CCBB44",
  "z-5'-4'",	          NA,	      NA,	      15,	      1,   	 	 "#CCBB44",
  "z-5'-4'-O",	        NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z-5'-4'-O-B-1'",	    NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z-5'-4'-O-B-1'-2'",	NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z-B",	              NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z-B-1'",	            NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z-B-1'-O",	          NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z-B-1'-O-2'",	      NA,	      NA,	      15,	      1,	   	 "#CCBB44",
  "z+2H",	              0.65,	    -0.035,	  15,	      1,	   	 "#CCBB44",
  "z+H",          	    0.65,	    -0.035,	  15,	      1,	   	 "#CCBB44"
)

#Functions for interpreting IDT nomenclature and mass calc----------------------

parse_IDT_sequence <- function(IDT_seq){
  parsed <- str_extract_all(IDT_seq, "[rm\\+]?([(]s[+\\-]{1}\\d+(\\.\\d+)?([(][12345]{1}'[+\\-]{1}\\d+(\\.\\d+)?[)]){0,5}[)])?[ACGTURYMKSWHBVDN]{1}([(]n[+\\-]{1}\\d+(\\.\\d+)?[)])?[*]{0,2}([(]b[+\\-]{1}\\d+(\\.\\d+)?[)])?") %>% 
    as_vector()
  
  parsed <- tibble(
    IDT_parsed = parsed,
    sugar = str_extract(parsed, "^[rm\\+]"),
    sugar_mod = str_extract(parsed, "[(]s[+\\-]{1}\\d+\\.\\d+[)]"),
    `1'` = str_extract(parsed, "[(]1'[+\\-]{1}\\d+(\\.\\d+)?[)]"),
    `2'` = str_extract(parsed, "[(]2'[+\\-]{1}\\d+(\\.\\d+)?[)]"),
    `3'` = str_extract(parsed, "[(]3'[+\\-]{1}\\d+(\\.\\d+)?[)]"),
    `4'` = str_extract(parsed, "[(]4'[+\\-]{1}\\d+(\\.\\d+)?[)]"),
    `5'` = str_extract(parsed, "[(]5'[+\\-]{1}\\d+(\\.\\d+)?[)]"),
    nitbase = str_extract(parsed, "[ACGTURYMKSWHBVDN]{1}"),
    nitbase_mod = str_extract(parsed, "([(]n[+\\-]{1}\\d+\\.\\d+[)])"), 
    backbone = str_extract(parsed, "[*]"),
    backbone_mod = str_extract(parsed, "([(]b[+\\-]{1}\\d+\\.\\d+[)])")
  )
  
  mass_table <- tibble(
    IDT_parsed = parsed$IDT_parsed,
    sugar = ifelse(!is.na(parsed$sugar), sugar_mass[parsed$sugar], 83.04969) + ifelse(!is.na(parsed$sugar_mod), as.numeric(str_extract(parsed$sugar_mod, "[+\\-]{1}\\d+\\.\\d+")), 0),
    `1'` = ifelse(!is.na(parsed$`1'`), as.numeric(str_extract(parsed$`1'`, "[+\\-]{1}\\d+\\.\\d+")), 0) + 13.00274,
    `2'` = ifelse(!is.na(parsed$`2'`), as.numeric(str_extract(parsed$`2'`, "[+\\-]{1}\\d+\\.\\d+")), 0) + ifelse(!is.na(parsed$sugar), sugar_mass[parsed$sugar]-83.04969, 1.00274) + 13.00274,
    `3'` = ifelse(!is.na(parsed$`3'`), as.numeric(str_extract(parsed$`3'`, "[+\\-]{1}\\d+\\.\\d+")), 0) + 13.00274,
    `4'` = ifelse(!is.na(parsed$`4'`), as.numeric(str_extract(parsed$`4'`, "[+\\-]{1}\\d+\\.\\d+")), 0) + 13.00274,
    `5'` = ifelse(!is.na(parsed$`5'`), as.numeric(str_extract(parsed$`5'`, "[+\\-]{1}\\d+\\.\\d+")), 0) + 14.00548,
    nitbase = nitbase_mass[parsed$nitbase] + ifelse(!is.na(parsed$nitbase_mod), as.numeric(str_extract(parsed$nitbase_mod, "[+\\-]{1}\\d+\\.\\d+")), 0),
    backbone = ifelse(!is.na(parsed$backbone), backbone_mass[parsed$backbone], 95.9612) + ifelse(!is.na(parsed$backbone_mod), as.numeric(str_extract(parsed$backbone_mod, "[+\\-]{1}\\d+\\.\\d+")), 0)
  )
  
  return(mass_table)
}
#prepared input sequence for manual proofreading and downstream calculations

calc_IDT_intactmass <- function(IDT_seq, fiveprime = 17.00274, threeprime = 17.00274){
  
  mass_table <- parse_IDT_sequence(IDT_seq)
  
  intact_mass <- sum(
    mass_table$sugar,
    mass_table$nitbase,
    mass_table$backbone[1:(nrow(mass_table)-1)],
    fiveprime,
    threeprime
  )
  return(intact_mass)
}
#calculates intact mass based on parsed IDT sequence using parse_IDT_sequence() above

calc_IDT_fragments <- function(mass_table, fiveprime = 17.00274, threeprime = 17.00274, add_plus_series = F, add_minus_series, xring_series = F){
  
    fragments <- tibble(
      fwd_cumsum = cumsum(mass_table$sugar + mass_table$nitbase + mass_table$backbone),
      back_cumsum = cumsum(rev(mass_table$sugar + mass_table$nitbase + mass_table$backbone))
    ) %>% 
      mutate(
        `a` = fiveprime + fwd_cumsum - mass_table$backbone - 1.007825035,
        `a-B` = fiveprime + fwd_cumsum - mass_table$backbone - mass_table$nitbase - 2*1.007825035,
        `b` = fiveprime + fwd_cumsum - mass_table$backbone + 15.99491463 + 1.007825035, #increased by 2 H
        `c` = fiveprime + fwd_cumsum - 15.99491463 - 1.007825035,
        `d` = fiveprime + fwd_cumsum + 1.007825035,
        `w` = threeprime + back_cumsum + 1.007825035,
        `x` = threeprime + back_cumsum - 15.99491463 - 1.007825035,
        `x-B` = threeprime + back_cumsum - 15.99491463 - rev(mass_table$nitbase) - 1.007825035,
        `y` = threeprime + back_cumsum - rev(mass_table$backbone) + 15.99491463 + 1.007825035,
        `z` = threeprime + back_cumsum - rev(mass_table$backbone) - 1.007825035#one high according to LNA
      ) %>% 
      select(-fwd_cumsum, -back_cumsum) %>% 
      mutate(
        position = 1:nrow(mass_table)
      ) %>% 
      .[-nrow(.),]
    
  if(add_plus_series == T){
    fragments <- fragments %>% 
      mutate(
        `a+1` = `a` + 1.008665,
        `a+2` = `a` + 2*1.008665,
        `a-B+1` = `a-B` + 1.008665,
        `a-B+2` = `a-B` + 2*1.008665,
        `b+1` = `b` + 1.008665,
        `b+2` = `b` + 2*1.008665,
        `c+1` = `c` + 1.008665,
        `c+2` = `c` + 2*1.008665,
        `d+1` = `d` + 1.008665,
        `d+2` = `d` + 2*1.008665,
        `w+1` = `w` + 1.008665,
        `w+2` = `w` + 2*1.008665,
        `x+1` = `x` + 1.008665,
        `x+2` = `x` + 2*1.008665,
        `x-B+1` = `x-B` + 1.008665,
        `x-B+2` = `x-B` + 2*1.008665,
        `y+1` = `y` + 1.008665,
        `y+2` = `y` + 2*1.008665,
        `z+1` = `z` + 1.008665,
        `z+2` = `z` + 2*1.008665
      ) 
  }
  if(add_minus_series == T){
    fragments <- fragments %>% 
      mutate(
        `a-1` = `a` - 1.008665,
        `a-2` = `a` - 2*1.008665,
        `a-B-1` = `a-B` - 1.008665,
        `a-B-2` = `a-B` - 2*1.008665,
        `b-1` = `b` - 1.008665,
        `b-2` = `b` - 2*1.008665,
        `c-1` = `c` - 1.008665,
        `c-2` = `c` - 2*1.008665,
        `d-1` = `d` - 1.008665,
        `d-2` = `d` - 2*1.008665,
        `w-1` = `w` - 1.008665,
        `w-2` = `w` - 2*1.008665,
        `x-1` = `x` - 1.008665,
        `x-2` = `x` - 2*1.008665,
        `x-B-1` = `x-B` - 1.008665,
        `x-B-2` = `x-B` - 2*1.008665,
        `y-1` = `y` - 1.008665,
        `y-2` = `y` - 2*1.008665,
        `z-1` = `z` - 1.008665,
        `z-2` = `z` - 2*1.008665
      ) 
  } 
    
  if(xring_series == T){
    fragments <- fragments %>% 
      mutate(
        `a-2'` = a - mass_table$`2'`[-nrow(.)],
        `a-3'` = a - mass_table$`3'`[-nrow(.)],
        `a-2'-3'` = a - mass_table$`2'`[-nrow(.)] - mass_table$`3'`[-nrow(.)],
        `a-B-1'` = `a-B` - mass_table$`1'`[-nrow(.)],
        `a-B-1'-O` = `a-B` - mass_table$`1'`[-nrow(.)] - 15.99491463,
        `a-B-1'-2'-3'-O` = `a-B` - mass_table$`1'`[-nrow(.)] - mass_table$`2'`[-nrow(.)] - mass_table$`3'`[-nrow(.)] - 15.99491463,
        `z-5'` = `z` - rev(mass_table$`5'`)[-nrow(.)],
        `z-5'-4'` = `z` - rev(mass_table$`5'`)[-nrow(.)] - rev(mass_table$`4'`)[-nrow(.)],
        `z-5'-4'-O` = `z` - (mass_table$`5'`)[-nrow(.)] - rev(mass_table$`4'`)[-nrow(.)] - 15.99491463,
        `z-B` = `z` - rev(mass_table$nitbase)[-nrow(.)],
        `z-5'-4'-O-B-1'` = `z-B` - rev(mass_table$`1'`)[-nrow(.)] - 15.99491463 - rev(mass_table$`5'`)[-nrow(.)] - rev(mass_table$`4'`)[-nrow(.)],
        `z-5'-4'-O-B-1'-2'` = `z-B` - rev(mass_table$`1'`)[-nrow(.)] - 15.99491463 - rev(mass_table$`2'`)[-nrow(.)] - rev(mass_table$`5'`)[-nrow(.)] - rev(mass_table$`4'`)[-nrow(.)],
        `z-B-1'` = `z-B` - rev(mass_table$`1'`)[-nrow(.)],
        `z-B-1'-O` = `z-B` - rev(mass_table$`1'`)[-nrow(.)] - 15.99491463,
        `z-B-1'-O-2'` = `z-B` - rev(mass_table$`1'`)[-nrow(.)] - 15.99491463 - rev(mass_table$`2'`)[-nrow(.)]
      )
  }
  
  return(fragments)
  
}
#calculates fragment masses for all fragment ion types based on parsed sequence

calculate_internal_fragments <- function(theor_fragments, intact_mass){
  
  threeprime_fragments <- theor_fragments %>% 
    pivot_longer(-position, names_to = "ion_type", values_to = "mass") %>% 
    filter(ion_type != "a-B") %>% 
    filter(ion_type != "x-B") %>% 
    filter(ion_type != "d-H2O") %>% 
    filter(ion_type != "w") %>% 
    filter(ion_type != "x") %>% 
    filter(ion_type != "y") %>% 
    filter(ion_type != "z")
  
  fiveprime_fragments <- theor_fragments %>% 
    pivot_longer(-position, names_to = "ion_type", values_to = "mass") %>% 
    filter(ion_type != "a-B") %>% 
    filter(ion_type != "x-B") %>% 
    filter(ion_type != "d-H2O") %>% 
    filter(ion_type != "a") %>% 
    filter(ion_type != "b") %>% 
    filter(ion_type != "c") %>% 
    filter(ion_type != "d") %>% 
    
  
  internal_fragments <- tibble(
    left_type = vector("character"),
    left_position = vector("integer"),
    right_type = vector("character"),
    right_position = vector("integer"),
    theoretical_mass = vector("numeric")
  )
  
  for(i in 1:nrow(threeprime_fragments)){
    for(j in 1:nrow(fiveprime_fragments)){
      
      if(((nrow(theor_fragments) + 1 - threeprime_fragments$position[i]) - fiveprime_fragments$position[j]) <= 2){
        next
      }
      
      internal_fragments <- internal_fragments %>% add_row(
        left_type = convert_iontype(threeprime_fragments$ion_type[i]),
        left_position = nrow(theor_fragments) + 1 - threeprime_fragments$position[i],
        right_type = convert_iontype(fiveprime_fragments$ion_type[j]),
        right_position = nrow(theor_fragments) + 1 - fiveprime_fragments$position[j],
        theoretical_mass = intact_mass - threeprime_fragments$mass[i] - fiveprime_fragments$mass[j]
      )
    }
  }
  
  return(internal_fragments)
}

ID_IDT_fragments <- function(theor_frags, mass_list, tolerance = 0.00001){
  
  theor_frags <- theor_frags %>% 
    pivot_longer(-position, names_to = "ion_type", values_to = "mass")
  
  DRNA_ID <- tibble(
    expr_mass = vector("numeric"),
    theor_mass = vector("numeric"),
    intensity = vector("numeric"),
    ion_type = vector("character"),
    `position` = vector("integer")
  )
  
  for(i in 1:nrow(theor_frags)){
    carrier <- mass_list %>% 
      filter(near(expr_mass, theor_frags$mass[i], tol = tolerance*theor_frags$mass[i]))
    if(nrow(carrier) > 0){
      DRNA_ID <- add_row(DRNA_ID,
                         `expr_mass` = carrier$expr_mass,
                         `theor_mass` = theor_frags$mass[i],
                         `intensity` = carrier$intensity,
                         `ion_type` = theor_frags$ion_type[i],
                         `position` = theor_frags$`position`[i])
    }
  }
  
  DRNA_ID <- DRNA_ID %>% 
    mutate(
      `ppm_error` = (theor_mass - expr_mass)/theor_mass * 1000000,
      `mass_diff` = theor_mass - expr_mass,
      `position` = `position`,
      `rounded_theor_mass` = round(theor_mass, 5)
    ) %>% 
    filter((position != 1) & (!is.na(position))) %>% 
    group_by(rounded_theor_mass) %>% 
    filter(intensity == pmax(intensity)) %>% 
    ungroup() %>% 
    select(-rounded_theor_mass) %>% 
    arrange(!desc(theor_mass))
  
  if(nrow(DRNA_ID) == 0){
    DRNA_ID <- tibble(
      expr_mass = NA,
      theor_mass = NA,
      intensity = NA,
      ion_type = NA,
      `position` = NA,
      ppm_error = NA,
      mass_diff = NA
    )
  }
  
  return(DRNA_ID)
}
#identifies fragments from a deconvoluted mass list against supplied theoretical fragments

filter_internal_fragments <- function(internal_fragments, theoretical_fragments, identified_fragments, keep_unambiguous = T, keep_terminal_IDs = T){
  
  filtered_internal_fragments <- internal_fragments
  
  if(keep_terminal_IDs == T){
    
    if(any(!is.na(identified_fragments))){
    
      complement_identified_fragments <- mutate(identified_fragments,
                                                `position` = nrow(theoretical_fragments) + 1 - position,
                                                ion_type = convert_iontype(ion_type)
      )
      
      filtered_internal_fragments <- 
        add_row(
          semi_join(internal_fragments, identified_fragments, by = c("left_type" = "ion_type", "left_position" = "position")),
          semi_join(internal_fragments, identified_fragments, by = c("right_type" = "ion_type", "right_position" = "position"))
        ) %>% 
        add_row(semi_join(internal_fragments, complement_identified_fragments, by = c("left_type" = "ion_type", "left_position" = "position"))) %>% 
        add_row(semi_join(internal_fragments, complement_identified_fragments, by = c("right_type" = "ion_type", "right_position" = "position"))) %>% 
        distinct()
    }
  }
  
  if(keep_unambiguous == T){
    
    filtered_internal_fragments <- filtered_internal_fragments %>% 
      mutate(
        theoretical_mass = round(theoretical_mass, 5)
      ) %>% 
      group_by(theoretical_mass) %>% 
      add_tally() %>% 
      filter(n == 1) %>% 
      select(-n) %>% 
      ungroup() %>% 
      arrange(desc(theoretical_mass))
    
  }
  
  return(filtered_internal_fragments)
}

batch_filter_internal_fragments <- function(internal_fragments, theoretical_fragments, batch_identified_fragments, keep_unambiguous = T, keep_terminal_IDs = T){
  
  set_filtered_internals <- tibble(
      left_type = vector("character"),
      left_position = vector("integer"),
      right_type = vector("character"),
      right_position = vector("integer"),
      theoretical_mass = vector("numeric")
    )
  
  for(i in 1:nrow(batch_identified_fragments)){
    temp_filtered <- filter_internal_fragments(internal_fragments, theoretical_fragments, batch_identified_fragments[[2]][[i]], keep_unambiguous, keep_terminal_IDs) %>% 
      mutate(
        expr_name = batch_identified_fragments$expr_name[i]
      )
    
    set_filtered_internals <- bind_rows(set_filtered_internals, temp_filtered)
  }
  
  set_filtered_internals <- set_filtered_internals %>% 
    group_by(expr_name) %>% 
    nest()
  
  return(set_filtered_internals)
}

ID_internal_fragments <- function(filtered_internal_fragments, identified_fragments, mass_list, tolerance = 0.00001){
  
  mass_list <- mass_list %>% 
    anti_join(identified_fragments, by = "expr_mass")
  
  DRNA_internal_IDs <- tibble(
    expr_mass = vector("numeric"),
    theor_mass = vector("numeric"),
    intensity = vector("numeric"),
    left_type = vector("character"),
    left_position = vector("integer"),
    right_type = vector("character"),
    right_position = vector("integer")
  )
  
  for(i in 1:nrow(filtered_internal_fragments)){
    carrier <- mass_list %>% 
      filter(near(expr_mass, filtered_internal_fragments$theoretical_mass[i], tol = tolerance*filtered_internal_fragments$theoretical_mass[i]))
    if(nrow(carrier) > 0){
      DRNA_internal_IDs <- add_row(DRNA_internal_IDs,
                                   expr_mass = carrier$expr_mass,
                                   theor_mass = filtered_internal_fragments$theoretical_mass[i],
                                   intensity = carrier$intensity,
                                   left_type = filtered_internal_fragments$left_type[i],
                                   left_position = filtered_internal_fragments$left_position[i],
                                   right_type = filtered_internal_fragments$right_type[i],
                                   right_position = filtered_internal_fragments$right_position[i]
      )
    }
  }
  
  DRNA_internal_IDs <- DRNA_internal_IDs %>% 
    mutate(
      `ppm_error` = (theor_mass - expr_mass)/theor_mass * 1000000,
      `mass_diff` = theor_mass - expr_mass,
    ) %>% 
    arrange(!desc(theor_mass))
  
  if(nrow(DRNA_internal_IDs) == 0){
    DRNA_internal_IDs <- tibble(
      expr_mass = NA,
      theor_mass = NA,
      intensity = NA,
      left_type = NA,
      left_position = NA,
      right_type = NA,
      right_position = NA
    )
  }
  
  return(DRNA_internal_IDs)
}

calc_seq_cov <- function(mass_table, IDT_IDs){
  length(
    unique(
      c(
        filter(IDT_IDs, str_detect(ion_type, "[abcd]"))$position,
        nrow(mass_table) - filter(IDT_IDs, str_detect(ion_type, "[wxyz]"))$position
      )
    )
  ) / (nrow(mass_table)-1) * 100
}
#calculates sequence coverage based on identified fragments

batch_seq_cov <- function(IDT_IDs, theor_frags){
  length(
    unique(
      c(
        filter(IDT_IDs, str_detect(ion_type, "[abcd]"))$position,
        (nrow(theor_frags) + 1) - filter(IDT_IDs, str_detect(ion_type, "[wxyz]"))$position
      )
    )
  ) / (nrow(theor_frags)) * 100
}

batch_IDT_frags <- function(theor_frags, MSexperiments, tolerance){
  IDs_carrier <- tibble(
    expr_name = vector("character"),
    expr_mass = vector("numeric"),
    theor_mass = vector("numeric"),
    intensity = vector("numeric"),
    ion_type = vector("character"),
    `position` = vector("integer")
  )
  
  for(i in 1:nrow(MSexperiments)){
    temp_mass_list <- MSexperiments[[2]][[i]]
    temp_frag_IDs <- ID_IDT_fragments(theor_frags, temp_mass_list, tolerance) %>% 
      mutate(
        expr_name = MSexperiments$expr_name[i]
      )
    IDs_carrier <- bind_rows(IDs_carrier, temp_frag_IDs)
  }
  
  IDs_carrier <- group_by(IDs_carrier, expr_name) %>% 
    nest() %>% 
    mutate(
      seq_cov = lapply(data, batch_seq_cov, theor_frags = theor_frags)
    ) %>% 
    mutate(
      seq_cov = lapply(seq_cov, signif, digits = 3)
    )
  
  return(IDs_carrier)
}#feed in a single sequence or a set of sequences, give back a tibble: cols(name, matched_fragments)

batch_internal_frags <- function(batch_internal_fragments, batch_IDs, MSexperiments, tolerance){
  
  internal_IDs_carrier <- tibble(
    expr_name = vector("character"),
    expr_mass = vector("numeric"),
    theor_mass = vector("numeric"),
    intensity = vector("numeric"),
    left_type = vector("character"),
    left_position = vector("integer"),
    right_type = vector("character"),
    right_position = vector("integer")
  )
  
  for(i in 1:nrow(MSexperiments)){
    
    temp_internal_IDs <- ID_internal_fragments(batch_internal_fragments[[2]][[i]], batch_IDs[[2]][[i]], MSexperiments[[2]][[i]], tolerance) %>% 
      mutate(
        expr_name = MSexperiments$expr_name[i]
      ) #%>% 
      # select(-ion_type, -position)
    internal_IDs_carrier <- bind_rows(internal_IDs_carrier, temp_internal_IDs)
  }
  
  internal_IDs_carrier <- group_by(internal_IDs_carrier, expr_name) %>% 
    nest()
  
  return(internal_IDs_carrier)
}


#Read-in support functions------------------------------------------------------
readin_expr_mass_list <- function(path, sheet){
  mass_list <- read_xlsx(path, sheet = sheet) %>% 
    filter(1:nrow(.) %in% (which(.[1] == "Mass"):nrow(.)+1))
  
  names(mass_list) <- c("Mass", "Intensity")
  
  mass_list <- mutate(mass_list, 
                      Mass = as.double(Mass),
                      Intensity = as.double(Intensity)
  )
  mass_list <- distinct(mass_list, Mass, .keep_all = T) %>% 
    filter(Intensity > 10) %>% 
    group_by(Mass) %>% #removes duplicates, takes highest intensity
    filter(Intensity == pmax(Intensity)) %>% 
    ungroup() %>% 
    as_tibble()
  
  return(mass_list)
}#return a single mass list


readin_MSexperiments <- function(path){
  lists_carrier <- tibble(
    expr_name = vector("character"),
    expr_mass = vector("numeric"),
    intensity = vector("numeric")
  )
  
  for(i in seq_along(excel_sheets(path))){
    
    temp_mass_list <- readin_expr_mass_list(path, i)
    
    lists_carrier <- add_row(lists_carrier,
                             expr_name = excel_sheets(path)[i],
                             expr_mass = temp_mass_list$Mass,
                             intensity = temp_mass_list$Intensity
    )
  }
  
  lists_carrier <- group_by(lists_carrier, expr_name) %>% 
    nest() %>% 
    mutate(expr_name = factor(expr_name))
  
  return(lists_carrier)
}#return a nested tibble of mass lists

convert_iontype <- function(ion_type){
  
  conversions <- vector("character")
  
  for(i in seq_along(ion_type)){
    if(ion_type[i] == "a"){
      conversions[i] <- "w"
    } else if(ion_type[i] == "b"){
      conversions[i] <- "x"
    } else if(ion_type[i] == "c"){
      conversions[i] <- "y"
    } else if(ion_type[i] == "d"){
      conversions[i] <- "z"
    } else if(ion_type[i] == "w"){
      conversions[i] <- "a"
    } else if(ion_type[i] == "x"){
      conversions[i] <- "b"
    } else if(ion_type[i] == "y"){
      conversions[i] <- "c"
    } else if(ion_type[i] == "z"){
      conversions[i] <- "d"
    } else {
      conversions[i] <- "NA" 
    }
  }
  
  return(conversions)
  
}

#Plotting-----------------------------------------------------------------------

gen_fragment_map <- function(parsed, IDT_IDs){
  fragment_map_data <- 
    left_join(IDT_IDs, map_params, by = "ion_type") %>%
    mutate(
      abs_position = ifelse(str_detect(ion_type, "[abcd]"), position, nrow(parsed) - position),
      x_coord = x_adj + (abs_position-1) %% 10,
      y_coord = y_adj - (abs_position-1) %/% 10
    ) %>% 
    filter(!is.na(x_adj))
  
  temp_map_params <- filter(map_params, map(map_params$ion_type, `==`, unique(IDT_IDs$ion_type)) %>% map(., any) %>% as_vector()) %>% filter(!is.na(x_adj)) %>% 
    arrange(str_sort(ion_type))
  
  
  figure <-
    ggplot(data = fragment_map_data)+
      geom_text(data = parsed, aes(x = (0:(nrow(parsed)-1) %% 10), y = -(0:(nrow(parsed)-1) %/% 10)), label = parsed$IDT_parsed)+
      geom_point(aes(x = x_coord, y = y_coord, color = ion_type,
                     text = sprintf(paste("</br>", ion_type, position, "</br>Mass: ", signif(expr_mass, 8), "</br>ppm Error: ", signif(ppm_error, 5), sep = ""))),
                 size = fragment_map_data$size, shape = fragment_map_data$shape
      )+
      scale_color_manual(labels = temp_map_params$ion_type, values = temp_map_params$color)+
      labs(color = "Ion Type", shape = "Ion Type")+
      xlim(c(-0.5, 10))+
      ylim(c(-((nrow(parsed)-1) %/% 10) - 0.5, 0.5))+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())+
      xlab("")+
      ylab("")
  
  return(figure)
} 

gen_error_plot <- function(IDT_IDs, tolerance = 10){
  IDT_IDs <-
    left_join(IDT_IDs, map_params, by = "ion_type") %>%
    select(expr_mass, theor_mass, ion_type, color, position, ppm_error, intensity) %>% 
    mutate(color = ifelse(is.na(color), "#000000", color))
    
  temp_map_params <- filter(map_params, map(map_params$ion_type, `==`, unique(IDT_IDs$ion_type)) %>% map(., any) %>% as_vector()) %>% 
    arrange(str_sort(ion_type))
  
  error_plot <- ggplot(data = IDT_IDs)+
    geom_point(aes(x = expr_mass, y = ppm_error, color = ion_type,
                   text = sprintf(paste("</br>", ion_type, position, "</br>Mass: ", signif(expr_mass, 8), "</br>Intensity: ", intensity, "</br>ppm Error: ", signif(ppm_error, 5), sep = ""))
                   ), size = 1)+
    scale_color_manual(labels = temp_map_params$ion_type, values = temp_map_params$color)+
    labs(color = "Ion Type")+
    ylim(c(-tolerance, tolerance))+
    ylab("Mass Error [ppm]")+
    xlab("Experimental Fragment Mass [Da]")
    
    return(error_plot)
}

gen_annotated_spectrum <- function(expr_mass_list, IDT_IDs, how_wide = 5){
  annotated_spectrum_data <- 
    full_join(IDT_IDs, expr_mass_list, by = c("expr_mass", "intensity")) %>% 
    left_join(map_params, by = "ion_type") %>% 
    select(expr_mass, theor_mass, intensity, ion_type, color, position, ppm_error) %>% 
    mutate(
      color = ifelse(is.na(color), "#000000", color)
    )
  
  temp_map_params <- filter(map_params, map(map_params$ion_type, `==`, unique(annotated_spectrum_data$ion_type)) %>% map(., any) %>% as_vector()) %>% 
    arrange(str_sort(ion_type))
  
  figure <- 
    ggplot(data = annotated_spectrum_data)+
      geom_col(aes(x = expr_mass, y = intensity, fill = ion_type,
                   text = sprintf(paste("</br>", ion_type, position, "</br>Mass: ", signif(expr_mass, 8), "</br>Intensity: ", intensity, "</br>ppm Error: ", signif(ppm_error, 5), sep = ""))
      ), position = "identity", width = how_wide)+
    scale_fill_manual(labels = temp_map_params$ion_type, values = temp_map_params$color)+
    labs(fill = "Ion Type")+
    xlab("Deconvoluted Mass")+
    ylab("Intensity")
  
  return(figure)
}

#NATMS App Body=================================================================

ui <- fluidPage(
  titlePanel("Nucleic Acid MS Deconvoluted Fragment Identification"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("sequence", "Sequence", value = "ACGTU"),
      numericInput("fiveprime", "5' Modification", value = "17.00274"),
      numericInput("threeprime", "3' Modification", value = "17.00274"),
      
      textOutput("sequence_mass"),#output intact mass of the strand
      
      tags$hr(),
      
      fileInput("experiment", "Deconvoluted Mass List"),
      
      selectInput("sheet_selection", "Experiment", "", selected = 1),
      
      numericInput("ppm_tolerance", "ppm Tolerance", value = 10, min = 1, max = 50, step = 1),
      
      checkboxInput("plus_series", "Add +1 and +2 Series", value = F),
      
      checkboxInput("minus_series", "Add -1 and -2 Series", value = F),
      
      checkboxInput("xring_series", "Search cross-ring cleavages", value = F),
      
      tags$hr(),
      
      downloadButton("xlsx_download_terminal", "Download Fragment IDs (.xlsx)")
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Sequence Parsing", 
                 tableOutput("parsed_sequence")
        ),
        tabPanel("Theoretical Fragments",
                 tableOutput("theoretical_fragments")
        ),
        tabPanel("Fragment Identification",
                 textOutput("FI_sequence_coverage"), 
                 textOutput("number_searched"),
                 textOutput("number_identified"),
                 # textOutput("FI_experiment_name"),
                 plotlyOutput("fragment_map", width = "900px"),
                 tableOutput("identified_fragments")
        ),
        tabPanel("Annotated Spectrum",
                 textOutput("AS_sequence_coverage"), 
                 # textOutput("AS_experiment_name"),
                 numericInput("how_wide", "Column Width", value = 5, min = 0.1, max = 20, step = 0.1),
                 plotlyOutput("annotated_spectrum", width = "800px")
        ),
        tabPanel("Matched Fragment Error",
                 textOutput("MFE_sequence_coverage"), 
                 # textOutput("MFE_experiment_name"),
                 plotlyOutput("ppm_error_plot", width = "700px")
        ),
        tabPanel("Internal Fragments",
                 checkboxInput("do_internal_fragments", "Search Internal Ions*", value = F),
                 downloadButton("xlsx_download_internal", "Download Internal Fragment IDs"),
                 tags$hr(),
                 textOutput("number_internals_searched"),
                 textOutput("number_internals_identified"),
                 tableOutput("identified_internals")
                 
        )
      )
    )
  )
)
server <- function(input, output){
  
  #Outputs----------------------------------------------------------------------
  output$sequence_mass <- renderText({
    req(input$sequence)
    paste("Monoisotopic Mass: ", signif(intact_mass(), 10), sep = "")
  })
  output$parsed_sequence <- renderTable({
    parsed_sequence()
  })
  output$theoretical_fragments <- renderTable({
    req(input$sequence)
    calculated_fragments()
  })
  output$identified_fragments <- renderTable({
    pulled_IDs()
  })
  # output$FI_experiment_name <- renderText({
  #   paste("Experiment: ", experiment_name())
  # })
  output$FI_sequence_coverage <- renderText({
    paste("Sequence Coverage: ", sequence_coverage(), "%")
  })
  # output$AS_experiment_name <- renderText({
  #   paste("Experiment: ", experiment_name())
  # })
  output$AS_sequence_coverage <- renderText({
    paste("Sequence Coverage: ", sequence_coverage(), "%")
  })
  # output$MFE_experiment_name <- renderText({
  #   paste("Experiment: ", experiment_name())
  # })
  output$MFE_sequence_coverage <- renderText({
    paste("Sequence Coverage: ", sequence_coverage(), "%")
  })
  output$identified_internals <- renderTable({
    pulled_internal_IDs()
  })
  output$number_searched <- renderText({
    paste("Masses Searched:", n_frags_searched())
  })
  output$number_identified <- renderText({
    n_frags_identified()
    paste("Fragments Identified:", n_frags_identified())
  })
  output$number_internals_searched <- renderText({
    paste("Internal Fragments Searched:", n_internals_searched())
  })
  output$number_internals_identified <- renderText({
    paste("Internal Fragments Identified", n_internals_identified())
  })
  
  output$fragment_map <- renderPlotly({
    ggplotly(
      gen_fragment_map(parsed_sequence(), pulled_IDs()),
      tooltip = "text"
    )
  })
  output$annotated_spectrum <- renderPlotly({
    ggplotly(
      gen_annotated_spectrum(pulled_experiment(),
                             pulled_IDs(),
                             input$how_wide),
      tooltip = "text"
    )
  })
  output$ppm_error_plot <- renderPlotly({
    ggplotly(
      gen_error_plot(pulled_IDs(), input$ppm_tolerance),
      tooltip = "text"
    )
  })
  
  output$xlsx_download_terminal <- downloadHandler(
    filename = function() {paste0("Fragment Identifications ", input$experiment$name, ".xlsx")},
    content = function(file) {
      terminal_prep <- batched_fragment_IDs()$data
      names(terminal_prep) <- paste(batched_fragment_IDs()$expr_name, " (", batched_fragment_IDs()$seq_cov, "%)", sep = "")
      write_xlsx(terminal_prep, path = file)
    }
  )
  
  output$xlsx_download_internal <- downloadHandler(
    filename = function() {paste0("Internal Fragment Identificaitons ", input$experiment$name, ".xlsx")},
    content = function(file) {
      internal_prep <- batched_internal_IDs()$data
      names(internal_prep) <- paste(batched_internal_IDs()$expr_name)
      write_xlsx(internal_prep, path = file)
    }
  ) 
  
  #Reactives--------------------------------------------------------------------
  intact_mass <- reactive({
    req(input$sequence)
    calc_IDT_intactmass(input$sequence, input$fiveprime, input$threeprime)
  })
  parsed_sequence <- reactive({
    req(input$sequence)
    parse_IDT_sequence(input$sequence)
  })
  calculated_fragments <- reactive({
    req(input$sequence)
    calc_IDT_fragments(parsed_sequence(), input$fiveprime, input$threeprime, input$plus_series, input$minus_series, input$xring_series)
  })
  calculated_internal_fragments <- reactive({
    req(input$sequence)
    req(input$do_internal_fragments)
    calculate_internal_fragments(calculated_fragments(), intact_mass())
  })
   batch_filtered_internal_fragments <- reactive({
    req(input$sequence)
    req(input$do_internal_fragments)
    batch_filter_internal_fragments(calculated_internal_fragments(), calculated_fragments(), batched_fragment_IDs())
  })
  readin_experiments <- reactive({
    req(input$experiment)
    readin_MSexperiments(input$experiment$datapath)
  })
  outVar <- reactive({
    req(input$experiment)
    experiments <- readin_experiments()$expr_name
  })
  observe({
    updateSelectInput(inputId = "sheet_selection",
                      choices = outVar()
                      )
  })
  experiment_to_numeric <- reactive({
    which(readin_experiments()$expr_name == input$sheet_selection)
  })
  pulled_experiment <- reactive({
    req(input$ppm_tolerance)
    req(input$sheet_selection)
    readin_experiments()[[2]][[experiment_to_numeric()]]
  })
  sequence_coverage <- reactive({
    signif(calc_seq_cov(parsed_sequence(), pulled_IDs()), 3)
  })
  # experiment_name <- reactive({
  #   req(input$ppm_tolerance)
  #   readin_experiments()$expr_name[experiment_to_numeric()]
  # })
  batched_fragment_IDs <- reactive({
    req(input$experiment)
    req(input$sequence)
    req(input$ppm_tolerance)
    batch_IDT_frags(calculated_fragments(), readin_experiments(), input$ppm_tolerance/1e6)
  })
  batched_internal_IDs <- reactive({
    req(input$experiment)
    req(input$sequence)
    req(input$ppm_tolerance)
    batch_internal_frags(batch_filtered_internal_fragments(), batched_fragment_IDs(), readin_experiments(), input$ppm_tolerance/1e6)
  }) 
  pulled_IDs <- reactive({
    req(input$ppm_tolerance)
    req(input$sheet_selection)
    batched_fragment_IDs()[[2]][[experiment_to_numeric()]]
  })
  pulled_internal_IDs <- reactive({
    req(input$ppm_tolerance)
    req(input$sheet_selection)
    batched_internal_IDs()[[2]][[experiment_to_numeric()]]
  })
  n_frags_searched <- reactive({
    req(pulled_IDs())
    nrow(readin_experiments()[[2]][[experiment_to_numeric()]])
  })
  n_frags_identified <- reactive({
    req(pulled_IDs())
    if(any(is.na(pulled_IDs()))){
      0
    } else {
      nrow(pulled_IDs()) 
    }
  })
  n_internals_searched <- reactive({
    req(pulled_internal_IDs)
    nrow(batch_filtered_internal_fragments()[[2]][[experiment_to_numeric()]])
  })
  n_internals_identified <- reactive({
    req(pulled_internal_IDs())
    if(any(is.na(pulled_internal_IDs()))){
      0
    } else {
      nrow(pulled_internal_IDs()) 
    }
  })
}
#Spinraza PT: rU*rC*rA*rC*rU*rU*rU*rC*rA*rU*rA*rA*rU*rG*rC*rU*rG*rG*
#Spinraza LNA: rU+C(n+14.01565)+A+C(n+14.01565)+T+T+T+C(n+14.01565)+A+T+A+A+T+G+C(n+14.01565)+T+G+G


#RUN APP========================================================================

shinyApp(ui, server)

