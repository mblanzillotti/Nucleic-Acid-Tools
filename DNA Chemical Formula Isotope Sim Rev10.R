require(tidyverse)
 #Imports "tidy" data functions, the basis of database manipulations in this program
require(scales)
 #Used for custom ggplot theme, not necessary
require(grid)
 #Used for custom ggplot theme, not necessary
require(ggthemes)
 #Used for custom ggplot theme, not necessary
require(IsoSpecR)
 #Generates fine isotope structures for each chemical formula, basis of isotope fitting
require(rawrr)
 #Enables read-in of .raw Thermo files. As of 13 June 2023, validate function must be overwritten
require(plotly)
 #Used for interactive plots, not necessary, but very useful.

#Changelog----------------------------------------------------------------------

# 06282022
#Increased ribobse and dioxyribose oxygen count to reflect 3' and 5' oxygens.
#Removed 2 oxygens from phosphate and thiophosphate compositions.
#Enables more specific modification with 3'+/- or 5'+/- nomenclature rather
#than having these oxygens be immutable as apart of the backbone.

# 06302022
#Updated fragment_compositions formula to accurately reflect a/b/c/d/w/x/y/z
#formulae.  Incorporated isotope generation, charge convolution, and fragment
#matching

# 06132023
#Lots of updates.
#Added hydrogen shift fitting support
#Altered isotope generation to reflect accurate mass deltas between peaks (not a neutron)
#Added more thorough annotations throughout
#Implemented separate usage script from functions
#Implemented fix for rawrr validation function in separate script
#In progress fix for combinations of losses and gains in chemical modification, for now use digit modification if applicable

#===============================================================================
#Presets and Global Variables
#===============================================================================

element_masses <- 
tribble(
  ~element, ~mono_mass,
  "C",      12.00000,
  "H",       1.00783,
  "N",      14.00307,
  "O",      15.99492,
  "P",      30.97376,
  "S",      31.97207,
  "Na",     22.98977
)
#Elements and associated monoisotopic masses for consideration. If a modification
#contains a different element than those listed, it must be included here and in
#the followin "symbol_composition" tribble below.

symbol_composition <- 
tribble(
  ~symbol,  ~C, ~H, ~N, ~O, ~P, ~S, ~Na,
  "r",      5,  7,  0,  2,  0,  0,  0,#ribose
  "d",      5,  7,  0,  1,  0,  0,  0,#deoxyribose
  "m",      6,  9,  0,  2,  0,  0,  0,#2'O-methylribose
  "+",      6,  7,  0,  2,  0,  0,  0,#Locked 4'-methyl-2' (I think)
  'A',      5,  4,  5,  0,  0,  0,  0,
  "C",      4,  4,  3,  1,  0,  0,  0,
  "G",      5,  4,  5,  1,  0,  0,  0,
  "T",      5,  5,  2,  2,  0,  0,  0,
  "U",      4,  3,  2,  2,  0,  0,  0,
  "*",      0,  1,  0,  1,  1,  1,  0,#thiophosphate backbone
  "p",      0,  1,  0,  2,  1,  0,  0,#phosphate backbone
  "Acetyl", 2,  2,  0,  1,  0,  0,  0,
  "Methyl", 1,  2,  0,  0,  0,  0,  0#,
) %>% 
  pivot_longer(-symbol, names_to = "element", values_to = "amount")
#Contains all pre-programmed modifications including sugars, bases, and backbone.
#Please add any commonly used modifications or any that cannot be easily implemented
#via the modification nomenclature s/n/b(+/-). Utilize symbols that play nice with
#stringr regular expressions (check out the stringr cheat-sheet).


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
#Helper which generates a tibble of the format [element, amount] from an input chemical formula

interpret_digits <- function(encoding){
  return(
    str_remove_all(encoding, "[\\(\\)]") %>% 
      str_extract("[\\-\\+]?\\d*\\.?\\d*") %>% 
      as.numeric()
  )
}
#Helper that extracts a net loss or gain in mass from a modification

interpret_symbol <- function(encoding){
  return(filter(symbol_composition,
                str_detect(symbol, paste("^", encoding, "$", sep = "")) 
    )%>% 
      select(-symbol)
  )
}
#Helper that extracts the associated chemical formula from a symbol in "symbol_compositions"

composition_to_mass <- function(composition){
  return(
    composition %>%
      left_join(element_masses, by = "element", multiple = "all") %>%
      summarise(sum(amount*mono_mass)) %>%
      as_vector() %>%
      unname()
  )
}
#Helper that converts a chemical formula tibble of the form [element, amount] to a monoisotopic mass

composition_to_nominal <- function(composition){
  return(
    composition %>%
      left_join(element_masses, by = "element", multiple = "all") %>%
      summarise(sum(amount*mono_mass)) %>%
      as_vector() %>%
      unname() %>% 
      round(0)
  )
}
#Helper that converts a chemical formula tibble of the form [element, amount] to an integer monoisotopic mass

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
#Interpreter function that selects which helper function to implement based on the parsed input,
#either determining the presence of a symbol, formula, or digit based on regular expressions.

composition_to_IsoSpec <- function(composition){
  IsoSpec_comp <- composition$amount %>% as_vector()
  names(IsoSpec_comp) <- composition$element %>% as_vector()
  
  return(IsoSpec_comp)
}
#Helper that converts a chemical formula tibble to a named vector to play nice with IsoSpec

IsoSpec_to_mass <- function(IsoSpec_comp){
  
  mass <- element_masses %>% mutate(
    total_mass = IsoSpec_comp[element] * mono_mass
  ) %>%
    filter(!is.na(total_mass))
  
  return(sum(mass$total_mass))
  
}
#Helper that converts a chemical formula vector to a monoisotopic mass

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
#Parsing function that takes an input of the form S(+/-MOD)N(+/-MOD)B(+/-MOD)...
#where S indicates a sugar (defaults omission to "d" = deoxyribose), N indicates
#a nitrogenous base (no default), and B indicates backbone (defaults omission to
#"p" = phosphate). Modifications to each component can be provided as digits (+/-15.995) representing
#net change in mass, chemical formula (+/-C6H12O6), or symbol (+/-Acetyl). Furthermore
#modifications on sugars can be specified at the carbon level (4'+O1C1H2).
#Utilization of presets of formulae will yield the most accurate isotopic fits.

get_nucleic_composition <- function(parsed_seq){
  
  parsed_seq %>% mutate(
    composition = map_chr(symbol, as_vector) %>% map(determine_composition)
  )
  
}
#Condenses a parsed nucleic acid composition to a single chemical formula

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
  
  H_adjust <- c(-1, -2, +1, -1, +1, +1, -1, -3, +1, -2)
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
#Generates a composition for a single parsed nucleic acid, position, and ion type
#while considering any empirical hydrogen shifts for the canonical backbone fragments
#and two commonly observed base-loss series. Implemented in looped form in
#"generate_frag_formulae"

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
      formula_mass = map(composition, composition_to_mass) %>% as.numeric(),
      total_mass = formula_mass + digit_mass
    )
  
  return(fragment_ion_compositions)
  
}
#Function that generates all possible fragment ions for a given parsed nucleic acid.
#This function iterates along the parsed sequence and through each position and ion type 
#to generate fragment ions. Each of these is put into a tidy tibble containing
#composition, position, type, and monoisotopic mass information.

##Isotope Generation------------------------------------------------------------

IsoSpec_nominal_dist <- function(IsoSpec_comp){
  
  nominal_dist <- IsoSpecR::IsoSpecify(IsoSpec_comp, 0.97) %>% 
    as_tibble()
  
  mass_vector <- IsoSpec_to_mass(IsoSpec_comp) + c(0:ceiling(max(nominal_dist$mass) - min(nominal_dist$mass)))*1.008665
  
  nominal_dist <- nominal_dist %>% 
    mutate(
      nuclide = map(mass, `-`, mass_vector) %>% map(abs) %>% map_int(which.min) - 1
    ) %>% 
    group_by(nuclide) %>% 
    summarise(
      nominal_prob = sum(prob)
    ) %>% 
    arrange(nuclide) %>% 
    mutate(
      nominal_mass = mass_vector[1:nrow(.)]
    ) %>%
    slice_max(order_by = nuclide, n = 15) %>% 
    ungroup()
  
  return(nominal_dist)
  
}
#Wraps around IsoSpecR::IsoSpecify to generate 97% abundance isotopic distributions
#and up to 15 peaks (modifiable in "slice_max")

generate_iso_dist <- function(fragment_compositions){
  
  fragment_distributions <- fragment_compositions %>% 
    mutate(
      nominal_isos = map(composition, composition_to_IsoSpec) %>% 
        map(IsoSpec_nominal_dist)
    ) %>% 
    unnest(cols = c(nominal_isos)) %>%
    mutate(                                        #Rather than the mass of a neutron (1.008...)
      exact_mass = nuclide*1.002185 + (total_mass) #mean(c(1.00422, 1.00003, 0.997035, 1.00336, 1.00628))
    ) %>%                                           #^ average mass difference between C, H, N, O isotopes
    select(ion_type, position, nominal_prob, exact_mass) %>% 
    nest(exact_isotopes = c(nominal_prob, exact_mass))
  
  return(fragment_distributions)
  
}
#Overlays abundance information generated from "IsoSpec_nominal_dist"
#on known monoisotopic mass information spaced by avergae isotopic spacing from
#C, N, H, and O isotopes.

charge_convolve <- function(fragment_distributions, min_charge = 1, max_charge, low_mz = 200, high_mz = 3000){
  
  fragment_mz <- fragment_distributions %>% 
    cross_join(tibble(z = min_charge:max_charge)) %>% 
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
#Based on maximum and minimum charge states and m/z range constraints, this produces a tibble
#of m/z domain isotope values, typically piped directly after "generate_iso_dist"

ID_centroid <- function(theoretical, experimental, tolerance = 1e-5){
  
  ID <- near(experimental, theoretical, tolerance*theoretical) %>% 
    which() %>% 
    experimental$expr_mz[.] %>% 
    .[which.min(
      abs(. - theoretical)
    )
    ] %>% 
    ifelse(identical(., integer(0)), NA, .)
  
  ID <- filter(centroids, expr_mz == ID)
  
  if(nrow(ID) == 0){
    
    return(
      tibble(
        expr_mz = NA,
        intensity = NA
      )
    )
    
  } else {
    
    return(ID)
    
  }
  
}
#To limit copy/pasting, this function simply gives an identification if a theoretical
#peak (one value) is matched within the tolerance of an experimental mass list


identify_hydrogen_shifts <- function(isotope_match, theor_frag, centroids, tolerance = 5e-6){
  #rewrite this incorporating theor_frag for situations where the monoisotopic peak isn't
  #present, and use head(1) and tail(1) notation rather than reordering all your dang tables
  theor_frag <- theor_frag %>% 
    arrange(mz) %>% 
    mutate(
      #nuclide = 0:(nrow(theor_frag)-1),
      h_shifts = 0
    )
  
  isotope_match <- isotope_match %>%
    arrange(mz) %>% 
    mutate(
      #nuclide = 0:(nrow(theor_frag)-1),
      h_shifts = 0
    )
  
  base_ion_type <- theor_frag$ion_type[1]
  charge_state <- theor_frag$z[1]
  
  if(str_detect(base_ion_type, "[awbxcydz]")){
    
    proposed_fragment <- isotope_match
    
    carrier <- isotope_match
    
    expected_shift <- -1
    
    i <- 0
    
    while(T){
      
      i <- i + expected_shift
      shift_mz <- (head(theor_frag, 1) %>% .$exact_mass + i*1.00783 + charge_state*1.00783) / abs(charge_state)
      shift_ID <- ID_centroid(shift_mz, centroids, tolerance)
      
      # if(!is.na(shift_ID$expr_mz) && shift_ID$intensity > carrier$intensity[1]*0.01){
      if(!is.na(shift_ID$expr_mz)){
        
        carrier <- carrier %>% 
          mutate(
            ion_type = paste(base_ion_type, i, sep = ""),
            exact_mass = exact_mass + expected_shift*1.00783,
            mz = (exact_mass + z*1.00783) / abs(z),
            #nuclide = nuclide + expected_shift,
            h_shifts = i,
            nominal_prob = nominal_prob,
            expr_mz = lag(expr_mz, default = shift_ID$expr_mz[1]),
            intensity = lag(intensity, default = shift_ID$intensity[1]),
            #mz_error = mz - expr_mz,
            #ppm_error = mz_error / mz,
          )
        
        proposed_fragment <- proposed_fragment %>% 
          add_row(carrier)
        
      } else {
        
        break
        
      }
      
      if(i == -3){
        
        break
        
      }
      
    }
    
  } else {
    
    return(
      isotope_match %>% 
        filter(!is.na(expr_mz))
    )
    
  }
  
  return(
    proposed_fragment %>% 
      arrange(mz) %>% 
      filter(!is.na(expr_mz))
  )
  
}
#This function takes an algorithmic approach to search hydrogen shifts (losses here)
#only identifying those masses, and passing any potential identifications to the
#following "fit_shift_contributions" function. Custom hydrogen shifts (e.g. +H)
#can be implemented, like for proteins.

fit_shift_contributions <- function(proposed_fragment){
  
  matrix_fragment <- proposed_fragment %>% 
    select(expr_mz, nominal_prob, ion_type, intensity) %>%  
    pivot_wider(names_from = ion_type, values_from = nominal_prob, values_fill = 0)
  
  modelled_intensities <- matrix(
    matrix_fragment %>%
      select(-c(1,2)) %>% 
      as.matrix(),
    byrow = F,
    nrow = nrow(matrix_fragment),
    ncol = ncol(matrix_fragment)-2
  ) %>% 
    qr() %>% 
    qr.solve(b = matrix_fragment$intensity)
  
  while(any(modelled_intensities < 0)){
    
    drops <- which(modelled_intensities < 0)
    
    matrix_fragment <- matrix_fragment %>%
      select(-(drops+2))
    
    modelled_intensities <- matrix(
      matrix_fragment %>%
        select(-c(1,2)) %>% 
        as.matrix(),
      byrow = F,
      nrow = nrow(matrix_fragment),
      ncol = ncol(matrix_fragment)-2
    ) %>%
      qr() %>%
      qr.solve(b = matrix_fragment$intensity)
    
  }
  
  names(modelled_intensities) <- matrix_fragment %>% names() %>% .[c(-1:-2)]
  
  return(
    proposed_fragment %>% 
      mutate(predicted_i = nominal_prob * modelled_intensities[ion_type]) %>% 
      filter(!is.na(predicted_i)) %>% 
      filter(!is.na(expr_mz))
  )
  
}
#Utilizes identified hydrogen shifts and treats resulting identified intensities grouped
#by experimental m/z as a system of linear equations. This function converts these values
#into a matrix and employs qr decomposition to solve for this function. If a negative
#component is observed for any one fragment, it is removed, and the fitting is re-run.

ID_zcon_frags <- function(fragments, centroids, tolerance = 1e-5){
  
  fragments <- fragments %>% 
    mutate(
      ion_type = factor(ion_type, levels = c("w", "a", "a-B", "d", "y", "x", "x-B", "z", "b", "c"))
    ) %>% 
    arrange(ion_type, desc(z)) %>% 
    mutate(ion_type = as.character(ion_type))
  
  matched_fragments <- tibble(
    ion_type = vector("character"),
    position = vector("integer"),
    z = vector("integer"),
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
      unnest(iso_mz) %>%
      arrange(mz) %>% 
      mutate(
        series_no = 1:nrow(.)
      ) %>% 
      arrange(desc(nominal_prob))
    
    isotope_match <- tibble(
      ion_type = vector("character"),
      position = vector("integer"),
      z = vector("integer"),
      exact_mass = vector("double"),
      mz = vector("double"),
      nominal_prob = vector("double"),
      series_no = vector("integer"),
      expr_mz = vector("double"),
      intensity = vector("double")
    )
    
    missed_nominal <- 0
    stop_matching <- F
    
    for(j in 1:nrow(current_fragment)){
      
      isotope_identification <- ID_centroid(current_fragment$mz[j], centroids, tolerance = tolerance)
      
      isotope_match <- add_row(isotope_match,
                               current_fragment %>% slice(j),
                               isotope_identification
      )
      
      if(is.na(isotope_identification$expr_mz[1])){
        
        missed_nominal <- current_fragment$nominal_prob[j] + missed_nominal
        
        if(missed_nominal > 0.6){
          stop_matching <- T
          break
        }
        
      }
      
    }
    
    if(stop_matching){
      
      next
      
    }
    
    isotope_match <- isotope_match %>% 
      arrange(mz) %>% 
      mutate(
        filtrate = is.na(expr_mz) & (lag(expr_mz) %>% is.na() %>% !.)
      ) %>% 
      filter(
        (series_no <= match(T, filtrate, nomatch = nrow(.)))
      ) %>% 
      select(-filtrate)
    
    if(nrow(isotope_match) > 0){
      
      if(isotope_match %>% filter(!is.na(expr_mz)) %>% .$nominal_prob %>% sum() > 0.60){
        
        candidate_fragment <- identify_hydrogen_shifts(isotope_match, current_fragment, centroids, tolerance = tolerance/2) %>% 
          fit_shift_contributions() %>%
          # evaluate_hydrogen_shifts(isotope_match) %>%
          mutate(
            mz_error = mz - expr_mz,
            ppm_error = mz_error / mz * 1e6,
          ) %>%
          filter(!is.na(predicted_i))
        # return(candidate_fragment)
        isotope_fit_check <- candidate_fragment %>%
          # mutate(
          #   mz_error = mz - expr_mz,
          #   ppm_error = mz_error / mz * 1e6,
          # ) %>% 
          group_by(expr_mz) %>% 
          summarise(
            predicted_i = sum(predicted_i),
            intensity = mean(intensity),
            i_delta = predicted_i - intensity,
            ppm_error = mean(ppm_error)
          ) %>% 
          ungroup() %>% 
          mutate(
            #prediction_cutoff = near(predicted_i, intensity, tol = predicted_i*0.2)
            prediction_cutoff = near(predicted_i, intensity, tol = sum(predicted_i)*0.035)#0.05 not as much
          )
        # return(isotope_fit_check)
        if(
          (sum(!isotope_fit_check$prediction_cutoff) <= ceiling(nrow(isotope_fit_check)*0.25)) &&
          all(abs(isotope_fit_check$i_delta) <= 0.1*sum(isotope_fit_check$intensity)) &&
          (weighted.mean(abs(isotope_fit_check$ppm_error), isotope_fit_check$intensity) < tolerance*1e6)
        ){
          
          matched_fragments <- add_row(matched_fragments,
                                       # select(candidate_fragment, -c(series_no, nuclide, h_shifts))
                                       select(candidate_fragment, -c(series_no, h_shifts))
          )
          
          centroids <- filter(centroids,
                              map(expr_mz, `!=`, candidate_fragment$expr_mz) %>% map_lgl(all)
          )
          
        }
        
      }
      
    }
    
  }
  
  
  return(matched_fragments)
}
#Large function incorporating many of the previous components to identify fragment ions
#from a theroetical list in an experimental spectrum. Isotope fitting parameters
#can be adjusted, but currently, 60% of predicted abundance must be accounted for
#peaks must be consecutive, and predicted abundance and actual abundance must be 
#within 3.5% of the total distribution's abundance. This weights towards tolerances
#at the most abundant peaks to be tighter than those at lower abundance.

#Usage==========================================================================

#parsed <- parse_nucleic_input("mG*mC*mG*rUrArGrArCrArCrGrGrArArGrArGrCrGrArGrUrUrUrUrArGrArGrCrUrArGrArArArUrArGrCrArArGrUrUrArArArArUrArArGrGrCrUrArGrUrCrCrGrUrUrArUrCrArArCrUrUrGrArArArArArGrUrGrGrCrArCrCrGrArGrUrCrGrGrUrGrCmU*mU*mU*rU")
parsed <- parse_nucleic_input("ATGCTGCCCGGG")

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
  charge_convolve(min_charge = -1, max_charge = -4, high_mz = 2000)

#setwd("D:/UT Austin/Research/DNA/Genentech/Representative sgRNA Top Down/ALK")
setwd("D:/UT Austin/Research/DNA/Jess")

spectrum <- readSpectrum("a1_-4_UVPD_3mj_1pulse_ms2_120k.raw", 1)
  
profile <- tibble(
  expr_mz = spectrum[[1]]$mZ,
  intensity = spectrum[[1]]$intensity,
)

centroids <- tibble(
    expr_mz = spectrum[[1]]$centroid.mZ,
    intensity = spectrum[[1]]$centroid.intensity,
  )


frag_IDs <- ID_zcon_frags(theor_frags, centroids, tolerance = 1.0e-5)

cooloors <- c("#4F81BD", "#4BACC6", "#C0504D", "#9BBB59", "#8064A2", "#2C4D75", "#772C2A", "#F79646", "#5F7530", "#4D3B62")
names(cooloors) <- c("a", "a-B", "b", "c", "d", "w", "x", "x-B", "y", "z")

ggplotly( 
  ggplot()+
    geom_line(data = profile,
              aes(x = expr_mz, y = intensity), color = "black"        
    )+
    geom_col(data = centroids,
             aes(x = expr_mz, y = intensity), color = "lightgrey", width = 0.0001
    )+
    geom_point(data = frag_IDs %>% 
                 mutate(ion_family = str_extract(ion_type, "[:alpha:]")) %>% 
                 group_by(expr_mz, ion_family, position, z) %>% 
                 summarise(predicted_i = sum(predicted_i)),
               aes(x = expr_mz, y = predicted_i, color = ion_family)
    )+
    geom_line(data = frag_IDs %>% 
                mutate(
                  type_pos_charge = paste(ion_type, position, z, sep = "|"),
                  ion_family = str_extract(ion_type, "[:alpha:]")
                ),
              aes(x = mz, y = predicted_i, color = ion_family, group = type_pos_charge)
    )+
    scale_color_manual(values = cooloors)+
    theme_publication()+
    labs(
      x = "m/z",
      y = "Intensity"
    )
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

