require(tidyverse) 
require(scales)
require(grid)
require(ggthemes)
require(IsoSpecR)
require(rawrr)
require(plotly)

#Changelog----------------------------------------------------------------------

# 06282022
#Increased ribobse and dioxyribose oxygen count to reflect 3' and 5' oxygens.
#Removed 2 oxygens from phosphate and thiophosphate compositions.
#Enables more specific modification with 3'+/- or 5'+/- nomenclature rather
#than having these oxyns be immutable as apart of the backbone.

# 06302022
#Updated fragment_compositions formula to accurately reflect a/b/c/d/w/x/y/z
#formulae.  Incorporated isotope generation, charge convolution, and fragment
#matching

#===============================================================================
#Presets and Global Variables
#===============================================================================

element_masses <- 
tribble(
  ~element, ~mono_mass,
  "C",      12.00000,
  "H",      1.00783,
  "N",      14.00307,
  "O",      15.99492,
  "P",      30.97376,
  "S",      31.97207,
  "Na",     22.98977
)

symbol_composition <- 
tribble(
  ~symbol,  ~C, ~H, ~N, ~O, ~P, ~S, ~Na,
  "r",      5,  7,  0,  2,  0,  0,  0,
  "d",      5,  7,  0,  1,  0,  0,  0,
  "m",      6,  9,  0,  2,  0,  0,  0,
  "+",      6,  7,  0,  2,  0,  0,  0,
  'A',      5,  4,  5,  0,  0,  0,  0,
  "C",      4,  4,  3,  1,  0,  0,  0,
  "G",      5,  4,  5,  1,  0,  0,  0,
  "T",      5,  5,  2,  2,  0,  0,  0,
  "U",      4,  3,  2,  2,  0,  0,  0,
  "*",      0,  1,  0,  1,  1,  1,  0,
  "p",      0,  1,  0,  2,  1,  0,  0, 
  "Acetyl", 2,  2,  0,  1,  0,  0,  0,
  "Methyl", 1,  2,  0,  0,  0,  0,  0#,
) %>% 
  pivot_longer(-symbol, names_to = "element", values_to = "amount")


#Functions======================================================================

##Helpers-----------------------------------------------------------------------

interpret_formula <- function(encoding){
  elemental_composition <- tibble(
    element = str_extract_all(encoding, "[:alpha:]{1}[:lower:]?\\d+") %>% 
      as_vector() %>% 
      str_extract("[:alpha:]{1}[:lower:]?"),
    amount = str_extract_all(encoding, "[:upper:]{1}[:lower:]?\\d+") %>% 
      as_vector() %>% 
      str_extract("\\d+") %>% 
      as.numeric()
  ) %>% 
    group_by(element) %>% 
    summarise(amount = sum(amount))
  
  return(elemental_composition)
}

interpret_digits <- function(encoding){
  return(
    str_remove_all(encoding, "[\\(\\)]") %>% 
      str_extract("[\\-\\+]?\\d*\\.?\\d*") %>% 
      as.numeric()
  )
}

interpret_symbol <- function(encoding){
  return(filter(symbol_composition,
                str_detect(symbol, paste("^", encoding, "$", sep = "")) 
    )%>% 
      select(-symbol)
  )
}

composition_to_mass <- function(composition){
  return(
    composition %>%
      left_join(element_masses, by = "element") %>%
      summarise(sum(amount*mono_mass)) %>%
      as_vector() %>%
      unname()
  )
}

composition_to_nominal <- function(composition){
  return(
    composition %>%
      left_join(element_masses, by = "element") %>%
      summarise(sum(amount*mono_mass)) %>%
      as_vector() %>%
      unname() %>% 
      round(0)
  )
}

determine_composition <- function(encoding){
  
  if(str_detect(encoding, "[:alpha:]") && str_detect(encoding, "[:digit:]")){
    #Letters and numbers means formula
    composition <- interpret_formula(encoding)
    
  } else if(str_detect(encoding, "[:graph:]") && !str_detect(encoding, "[:digit:]")){
      
    if(str_detect(encoding, "[\\*\\+]")){
      
      composition <- interpret_symbol(paste("\\", encoding, sep = ""))
      
    } else {
      
      composition <- interpret_symbol(encoding)
      
    }
    
  } else if(!str_detect(encoding, "[:alpha:]") && str_detect(encoding, "[:digit:]")){
    
    composition <- interpret_digits(encoding)
    
  } else {
    
    composition <- "Invalid input"
    
  }
  
  return(composition)
  
}

composition_to_IsoSpec <- function(composition){
  IsoSpec_comp <- composition$amount %>% as_vector()
  names(IsoSpec_comp) <- composition$element %>% as_vector()
  
  return(IsoSpec_comp)
} 

##Sequence Handling-------------------------------------------------------------

parse_nucleic_input <- function(input_seq){
  
  parsed <- str_extract_all(input_seq, "[rmd\\+]?(\\([12345]?'?[\\+\\-]?[:alnum:]+\\.?[:alnum:]*\\)){0,5}[ACGTU]{1}(\\([\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\))?[p*]?(\\([\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\))?") %>% 
    as_vector()
  
  parsed <- tibble(
    IDT_parsed = parsed,
    position = 1:length(parsed),
    sugar = str_extract(parsed, "^[rmd\\+]"),
    sugar_mod = str_extract(parsed, "[rmd\\+](\\([12345]?'?[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)){1,5}") %>% str_remove("[rmd\\+]"),
    nitbase = str_extract(parsed, "[ACGTU]{1}"),
    nitbase_mod = str_extract(parsed, "[ACGTU]{1}(\\([\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\))") %>% str_remove("[ACGTU]"), 
    backbone = str_extract(parsed, "[p*]"),
    backbone_mod = str_extract(parsed, "[p*](\\([\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\))") %>% str_remove("[p*]")
  ) %>% 
    mutate(
      sugar = ifelse(is.na(sugar), "d", sugar),
      backbone = ifelse(is.na(backbone), "p", backbone) %>% ifelse(position == last(position), NA, .),
      `1'` = str_extract(sugar_mod, "\\(1'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)"),
      `2'` = str_extract(sugar_mod, "\\(2'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)"),
      `3'` = ifelse(is.na(str_extract(sugar_mod, "\\(3'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)")),
               "O1", str_extract(sugar_mod, "\\(3'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)")
             ) %>% ifelse(position == last(position), paste(., "H1", sep = ""), .),
      `4'` = str_extract(sugar_mod, "\\(4'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)"),
      `5'` = ifelse(is.na(str_extract(sugar_mod, "\\(5'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)")),
                    "O1", str_extract(sugar_mod, "\\(5'[\\+\\-]?[:alnum:]+\\.?[:alnum:]+\\)")
      ) %>% ifelse(position == first(position), paste(., "H1", sep = ""), .),
    ) %>% 
    select(!IDT_parsed) %>% 
    pivot_longer(-position, values_to = "symbol", names_to = "item") %>% 
    filter(!is.na(symbol))
  
  return(parsed)
}

get_nucleic_composition <- function(parsed_seq){
  
  parsed_seq %>% mutate(
    composition = map_chr(symbol, as_vector) %>% map(determine_composition)
  )
  
}

fragment_composition <- function(parsed_input, pos, ion_type){
  
  if(str_detect(ion_type, "^a$")){
    
    carrier <- parsed_input %>% 
      filter(
        (position < pos) | (position == pos & !str_detect(item, "(3')|(backbone)"))
      ) 
    
  } else if(str_detect(ion_type, "^a\\-B$")){
    
    carrier <- parsed_input %>% 
      filter(
        (position < pos) | (position == pos & !str_detect(item, "(3')|(backbone)|(nitbase)"))
      )
    
  } else if(str_detect(ion_type, "^b$")){
    
    carrier <- parsed_input %>% 
      filter(
        (position < pos) | (position == pos & !str_detect(item, "(backbone)"))
      )
      
  } else if(str_detect(ion_type, "^c$")){
    
    carrier <- parsed_input %>% 
      filter(
        (position <= pos)
      )
      
  } else if(str_detect(ion_type, "^d$")){
    
    carrier <- parsed_input %>% 
      filter(
        (position <= pos) | (position == (pos + 1) & str_detect(item, "(5')"))
      )
    
  } else if(str_detect(ion_type, "^w$")){
    
    carrier <- parsed_input %>%
       filter(
         (position > pos) | ((position == pos) & str_detect(item, "(3')|(backbone)"))
       )
    
  } else if(str_detect(ion_type, "^x$")){
    
    carrier <- parsed_input %>%
      filter(
        (position > pos) | ((position == pos) & str_detect(item, "(backbone)"))
      )
      
  } else if(str_detect(ion_type, "^x\\-B$")){
    
    carrier <- parsed_input %>%
      filter(
        (position > pos) | ((position == pos) & str_detect(item, "(backbone)|(nitbase)"))
      )
      
  } else if(str_detect(ion_type, "^y$")){
    
    carrier <- parsed_input %>%
      filter(
        (position > pos)
      )

  } else if(str_detect(ion_type, "^z$")){
    
    carrier <- parsed_input %>%
      filter(
        (position > pos) | ((position == pos) & !str_detect(item, "(5')"))
      )
    
  } else {
    
    return("Please input a valid ion type")
    
  }
  
  H_adjust <- c(-1, -1, +1, -1, +1, +1, -1, -3, +1, -1)
  names(H_adjust) <- c("a", "a-B", "b", "c", "d", "w", "x", "x-B", "y", "z")
  
  composition_carrier <- carrier %>%
    filter(!map_lgl(composition, is.numeric)) %>% 
    unnest(composition) %>%
    select(element, amount) %>%
    group_by(element) %>%
    summarise(amount = sum(amount)) %>% 
    mutate(
      amount = ifelse(str_detect(element, "H"), amount+H_adjust[ion_type], amount)
    )

  digits_carrier <- carrier %>%
    filter(map_lgl(composition, is.numeric)) %>%
    mutate(composition = map_dbl(composition, unlist)) %>% 
    .$composition %>%
    as.vector() %>%
    sum()

  frag_comp <- list(
    composition = composition_carrier,
    digits = digits_carrier
  )
  
  return(
    frag_comp
  )
  
}

generate_frag_formulae <- function(parsed_input, ion_types = c("a", "a-B", "b", "c", "d", "w", "x", "x-B", "y", "z")){
  
  fragment_ion_compositions <- tibble(
    ion_type = vector("character"),
    position = vector("integer"),
    digit_mass = vector("numeric"),
    element = vector("character"),
    amount = vector("numeric")
  )
  
  for(i in 1:(max(parsed_input$position)-1)){
    
    for(j in seq_along(ion_types)){
      
      carrier <- fragment_composition(
        parsed_input,
        i,
        ion_types[j]
      )
      
      if(typeof(carrier) == "character"){
        
        next
        
      } else {
        
        fragment_ion_compositions <- fragment_ion_compositions %>% 
          add_row(
            ion_type = ion_types[j],
            position = i,
            digit_mass = carrier[[2]],
            element = carrier[[1]]$element,
            amount = carrier[[1]]$amount
          )
        
      }
      
    }
  }
  
  fragment_ion_compositions <- fragment_ion_compositions %>% 
    group_by(ion_type, position) %>%
    nest(composition = c("element", "amount")) %>% 
    mutate(
      formula_mass = map(composition, composition_to_mass) %>% as.numeric()
    )
  
  return(fragment_ion_compositions)
  
}

##Isotope Generation------------------------------------------------------------

IsoSpec_nominal_dist <- function(IsoSpec_comp){
  
  nominal_dist <- IsoSpecR::IsoSpecify(IsoSpec_comp, 0.95) %>% 
    as_tibble() %>% 
    mutate(
      nominal_mass = round(mass, 0)
    ) %>% 
    group_by(nominal_mass) %>% 
    summarise(
      nominal_prob = sum(prob)
    ) %>% 
    select(nominal_mass, nominal_prob) %>% 
    slice_max(order_by = nominal_mass, n = 15) %>% 
    ungroup()
  
  return(nominal_dist)
  
}

generate_iso_dist <- function(fragment_ion_compositions){
  
  fragment_ion_distributions <- fragment_ion_compositions %>% 
    mutate(
      nominal_formula = round(formula_mass, 0),
      nominal_isos = map(composition, composition_to_IsoSpec) %>%
        map(IsoSpec_nominal_dist)
    ) %>% 
    unnest(cols = c(nominal_isos)) %>%
    mutate(
      exact_mass = (nominal_mass - nominal_formula)*1.008665 + formula_mass + digit_mass,
      mono_mass = formula_mass + digit_mass
    ) %>% 
    select(ion_type, position, nominal_prob, exact_mass, mono_mass) %>% 
    nest(exact_isotopes = c(nominal_prob, exact_mass))
  
  return(fragment_ion_distributions)
}

charge_convolve <- function(fragment_distributions, min_charge = 1, max_charge, low_mz = 200, high_mz = 3000){
  
  fragment_mz <- fragment_distributions %>% 
    full_join(tibble(z = min_charge:max_charge), by = character(0)) %>% 
    unnest(cols = exact_isotopes) %>% 
    mutate(
      mz = (exact_mass + z*1.00727647)/abs(z)
    ) %>% 
    ungroup() %>% 
    group_by(ion_type, position, z) %>% 
    filter(mz > low_mz) %>% 
    filter(mz < high_mz) %>% 
    nest(iso_mz = c(exact_mass, mz, nominal_prob)) %>% 
    ungroup()
  
  return(fragment_mz)
}

ID_zcon_frags <- function(fragments, centroids, tolerance = 0.00001){
  
  matched_fragments <- tibble(
    ion_type = vector("character"),
    position = vector("integer"),
    charge = vector("integer"),
    exact_mass = vector("numeric"),
    mz = vector("numeric"),
    nominal_prob = vector("numeric"),
    expr_mz = vector("numeric"),
    intensity = vector("numeric"),
    mz_error = vector("numeric"),
    ppm_error = vector("numeric"),
    predicted_i = vector("numeric")
  )
  
  for(i in 1:nrow(fragments)){
    
    current_fragment <- fragments %>% 
      slice(i) %>% 
      unnest(iso_mz)
    
    isotope_matches <- tibble(
      ion_type = vector("character"),
      position = vector("integer"),
      charge = vector("integer"),
      exact_mass = vector("numeric"),
      mz = vector("numeric"),
      nominal_prob = vector("numeric"),
      expr_mz = vector("numeric"),
      intensity = vector("numeric"),
      series_no = vector("integer")
    )
    
    for(j in 1:nrow(current_fragment)){
      
      ID <- filter(centroids, near(centroids$expr_mz, current_fragment$mz[j], tol = current_fragment$mz[j]*(tolerance+5e-6)))
      
      if(nrow(ID) > 0 && min(ID$intensity) > 10){
        
        ID <- mutate(ID,
                     mz_error = current_fragment$mz[j] - expr_mz,
        ) %>% 
          #filter(abs(mz_error) == min(abs(mz_error)))
          filter(intensity == max(intensity))
        
        isotope_matches <- add_row(isotope_matches,
                                   ion_type = current_fragment$ion_type[j],
                                   position = current_fragment$position[j],
                                   charge = current_fragment$z[j],
                                   exact_mass = current_fragment$exact_mass[j],
                                   mz = current_fragment$mz[j],
                                   nominal_prob = current_fragment$nominal_prob[j],
                                   expr_mz = ID$expr_mz,
                                   intensity = ID$intensity,
                                   series_no = j
        )
        
      }
      
    }
    
    # #Diagnostic to evaluate ppm matching
    # if(sum(isotope_matches$nominal_prob) > 0.6){
    # 
    #   matched_fragments <- add_row(matched_fragments,
    #                                  select(isotope_matches, -series_no),
    #                                  predicted_i = sum(isotope_matches$intensity)*nominal_prob,
    #                                  mz_error = mz - expr_mz,
    #                                  ppm_error = mz_error / mz * 1e6
    #                        )
    # 
    # 
    # }
    
    if(nrow(isotope_matches) > 0){

      if(sum(isotope_matches$nominal_prob) > 0.70){

        isotope_fit_check <- mutate(isotope_matches,
                                    normalized_intensity = intensity / max(intensity),
                                    normalized_prob = nominal_prob / max(nominal_prob),
                                    adjusted_prob = nominal_prob / sum(nominal_prob),
                                    frac_intensity = intensity / sum(intensity),
                                    prob_diff = abs(normalized_prob - normalized_intensity),
                                    #frac_diff = nominal_prob - frac_intensity,
                                    frac_diff = adjusted_prob - frac_intensity,
                                    mz_error = mz - expr_mz,
                                    ppm_error = mz_error / mz * 1e6,
                                    series_no = series_no,
                                    #predicted_i = max(intensity)*normalized_prob,
                                    predicted_i = adjusted_prob*sum(intensity),
                                    predicted_i_diff_pct = abs((predicted_i - intensity) / ((predicted_i + intensity)/2))
        )

        # distribution_check <- tibble(
        #   delta_intensity = diff(isotope_fit_check$intensity),
        #   delta_predicted_i = diff(isotope_fit_check$predicted_i),
        #   raw_delta = delta_intensity - delta_predicted_i,
        #   model_agreement = abs(delta_intensity/delta_predicted_i)
        # )

        if(all(abs(isotope_fit_check$frac_diff) < 0.1) &&
           all(abs(diff(isotope_fit_check$series_no)) == 1) &&
           #(sum(abs(isotope_fit_check$predicted_i_diff_pct) > 1) <= 1) &&
           all((abs(isotope_fit_check$predicted_i_diff_pct) <= 1.5)) &&
           (abs(mean(isotope_fit_check$ppm_error)) < tolerance*1e6) #&&
           
        ){

          matched_fragments <- add_row(matched_fragments,
                                       select(isotope_matches, -series_no),
                                       select(isotope_fit_check, predicted_i, mz_error, ppm_error)
          )

        }

      }

    }
    
  }
  
  
  return(matched_fragments)
}

#Usage==========================================================================

parsed <- parse_nucleic_input("mG*mC*mG*rUrArGrArCrArCrGrGrArArGrArGrCrGrArGrUrUrUrUrArGrArGrCrUrArGrArArArUrArGrCrArArGrUrUrArArArArUrArArGrGrCrUrArGrUrCrCrGrUrUrArUrCrArArCrUrUrGrArArArArArGrUrGrGrCrArCrCrGrArGrUrCrGrGrUrGrCmU*mU*mU*rU")

parsed %>% 
  get_nucleic_composition %>% 
  unnest(composition) %>% 
  group_by(element) %>% 
  summarise(amount = sum(amount)) %>% 
  composition_to_mass()

theor_frags <- parsed %>% 
  get_nucleic_composition() %>% 
  generate_frag_formulae() %>% 
  generate_iso_dist() %>% 
  charge_convolve(min_charge = -1, max_charge = -30, high_mz = 2000)

setwd("D:/UT Austin/Research/DNA/Genentech/Representative sgRNA Top Down/ALK")

spectrum <- readSpectrum("20211105_chenb_sgRNA_ALK_CID_30-_16NCE_p25Q Scan[15].raw", 1)
  
profile <- tibble(
  expr_mz = spectrum[[1]]$mZ,
  intensity = spectrum[[1]]$intensity,
)

centroids <- tibble(
    expr_mz = spectrum[[1]]$centroid.mZ,
    intensity = spectrum[[1]]$centroid.intensity,
  )


frag_IDs <- ID_zcon_frags(theor_frags, centroids, tolerance = 1.5e-5)

cooloors <- c("#4F81BD", "#4BACC6", "#C0504D", "#9BBB59", "#8064A2", "#2C4D75", "#772C2A", "#F79646", "#5F7530", "#4D3B62")
names(cooloors) <- c("a", "a-B", "b", "c", "d", "w", "x", "x-B", "y", "z")

ggplotly( 
  ggplot()+
    geom_line(data = profile, aes(x = expr_mz, y = intensity#, 
                                  #text = sprintf(paste("</br>/z: ", expr_mz, "</br>Intensity: ", intensity, sep = ""))
                                  ))+
    geom_col(data = centroids, aes(x = expr_mz, y = intensity), width = 0.001)+
    geom_point(data = frag_IDs, aes(x = mz, y = predicted_i, color = ion_type#,
                                    #text = sprintf(paste("</br>", ion_type, position, ", ", charge, "</br>m/z: ", signif(expr_mz, 8), "</br>Intensity: ", intensity, "</br>ppm Error: ", signif(ppm_error, 5), sep = ""))
                                    )
               )+
    scale_color_manual(values = cooloors)+
    theme_publication()+
    labs(
      #x = expression(paste(bolditalic("m/z"))),
      x = "m/z",
      y = "Intensity",
      color = "Ion Type"
    )#,
  #tooltip = "text"
)


ion_type_graph <- c(1, 1, 2, 3, 4, -1, -2, -2, -3, -4)
names(ion_type_graph) <- c("a", "a-B", "b", "c", "d", "w", "x", "x-B", "y", "z")


  ggplot()+
    geom_text(
      data = parsed %>% 
        pivot_wider(names_from = item, 
                    values_from = symbol, 
                    values_fill = ""
        ) %>% 
        mutate(
          n_side = paste(sep = "", sugar, nitbase, backbone)
        ) %>% 
        select(position, n_side),
      aes(
        x = (position - 1) %% 10, y = -((position - 1) %/% 10), label = n_side
      ),
      size = 5,
      hjust = 0.5,
      vjust = 0.5
    )+
    geom_text(
      data = parsed %>% 
        pivot_wider(names_from = item, 
                    values_from = symbol, 
                    values_fill = ""
        ) %>% 
        filter(position == 1 | (position %% 10 == 1)),
      aes(
        x = -1.15, y = -((position - 1) %/% 10), label = position
      ),
      color = "lightgrey",
      size = 4,
      hjust = 0.5,
      vjust = 0.5
    )+
    geom_text(
      data = frag_IDs %>%
        select(ion_type, position) %>%
        distinct() %>%
        mutate(
          pos_x = ((position - 1) %% 10) + abs(ion_type_graph[ion_type])*0.125 + 0.24 %>% 
            ifelse(
              (position %% 10 == 0) & str_detect(ion_type, "[wxyz]"),
              . - 10,
              .
            ),
          pos_y = -((position - 1) %/% 10) + ifelse(str_detect(ion_type, "[abcd]"), 0.15, -0.25) %>% 
            ifelse(
              (position %% 10 == 0) & str_detect(ion_type, "[wxyz]"),
              . - 1,
              .
            )
        ),
      aes(
        angle = ifelse(str_detect(ion_type, "[abcd]"), 180, 0),
        x = ifelse(str_detect(ion_type, "[abcd]"), pos_x-0.09, pos_x),
        y = pos_y,
        color = ion_type
      ),
      label = "L", size = 4.5, key_glyph = "point"
    )+
    # geom_point(
    #   data = frag_IDs %>% 
    #     select(ion_type, position) %>% 
    #     distinct() %>% 
    #     mutate(
    #       pos_x = ((position - 1) %% 10) + abs(ion_type_graph[ion_type])*0.05 + 0.35,
    #       pos_y = -((position - 1) %/% 10) + ifelse(str_detect(ion_type, "[abcd]"), 0.30, -0.10)
    #     ),
    #   aes(
    #     x = pos_x,
    #     y = pos_y,
    #     color = ion_type,
    #     #shape = ifelse(str_detect(ion_type, "[abcd]"), 209, 211),
    #     text = sprintf(paste("</br>Ion Type: ", ion_type, position, sep = ""))
    #   ), size = 5, shape = 76  #3, 142 plotly | 5, 124 ggplot
    # )+
    ylim(-10, 0.25)+
    labs(
      color = "Ion Type"
    )+
    scale_shape_identity()+
    scale_color_manual(values = cooloors)+
    theme_void(base_size = 14)+
    theme(
      panel.grid.major = element_blank(),
      axis.line = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      legend.spacing = unit(0, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin = unit(c(2,2,2,2),"mm")
    )+
    guides(color = guide_legend(nrow = 1))

# ggsave(
#   "F:/UT Austin/Research/DNA/Coding/Spinraza Thiophosphate Data/05122021_SpinrazaPT_FragmentMap",
#   width = 7,
#   height = 4,
#   dpi = 220,
#   units = "in",
#   device = "tiff"
# )  
  
# dat <- tibble(p = c(0:200),
#               x = p %% 16,
#               y = p %/% 16)
# 
# 
# ggplotly(
#   ggplot(dat, aes(x, y)) +
#     geom_text(aes(label = p), size = 3, nudge_y = -.25) +
#     geom_point(aes(shape = p), size = 5, fill = "red") +
#     scale_shape_identity() +
#     theme_void()
# )

