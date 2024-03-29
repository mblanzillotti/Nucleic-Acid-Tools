#Nucleic Acid Fragment Fitting Usage============================================

#sgRNA ALK
#parsed <- parse_nucleic_input("mG*mC*mG*rUrArGrArCrArCrGrGrArArGrArGrCrGrArGrUrUrUrUrArGrArGrCrUrArGrArArArUrArGrCrArArGrUrUrArArArArUrArArGrGrCrUrArGrUrCrCrGrUrUrArUrCrArArCrUrUrGrArArArArArGrUrGrGrCrArCrCrGrArGrUrCrGrGrUrGrCmU*mU*mU*rU")

#miR145
#parsed <- parse_nucleic_input("rGrUrCrCrArGrUrUrUrUrCrCrCrArGrGrArArUrCrCrCrU")

#60mer DNA
#parsed <- parse_nucleic_input("GTAACCAAGAGTATTCCATTTAAGAGTATTCCATTTACCAAGAGTCCATTTACCAAGAAG")

#GQP T4(G4T4)4
#parsed <- parse_nucleic_input("TTTTGGGGTTTTGGGGTTTTGGGGTTTTGGGGTTTT")
#parsed <- parse_nucleic_input("TTTTG(K1)G(H-1)GGTTTTG(K1)G(H-1)GGTTTTG(K1)G(H-1)GGTTTTGGGGTTTT")

parsed <- parse_nucleic_input("r(P1O3H2)GrCrGrGrArUrUrUrArG(Methyl)rCrUrCrArGrDrDrGrGrGrArGrArGrCrG(Dimethyl)rCrCrArGrAmCrUmGrArArYrArUrC(Methyl)rUrGrGrArGrG(Methyl)rUrCrC(Methyl)rUrGrUrGrTrUrCrGrA(Methyl)rUrCrCrArCrArGrArArUrUrCrGrCrArCrCrA")
#  Gene sequence, m/z 890
#"parsed" here will output a tibble broken down by backbone position, interpreting
#each chunk of the input sequence


parsed %>% 
  get_nucleic_composition() %>% 
  unnest(composition) %>%
  group_by(element) %>% 
  summarise(amount = sum(amount)) %>% 
  composition_to_mass()
#This set of functions outputs the monoisotopic mass of the parsed sequence

M_isotopes <- parsed %>% 
  get_nucleic_composition() %>% 
  unnest(composition) %>% 
  group_by(element) %>% 
  summarise(amount = sum(amount)) %>% 
  composition_to_IsoSpec() %>% 
  IsoSpec_nominal_dist() %>% view()
#This generates isotopes for the intact, parsed sequence


theor_frags <- parsed %>% 
  get_nucleic_composition() %>% 
  generate_frag_formulae() %>% 
  generate_iso_dist() %>% 
  charge_convolve(min_charge = -1, max_charge = -28, high_mz = 2000) %>% 
  mutate(
    charge_proportion = position/z - 75/28 %>% abs()
  ) %>% 
  group_by(ion_type, position) %>% 
  slice_min(charge_proportion, n = 12) %>% 
  select(!charge_proportion) %>% ungroup()
#This generates a set of theoretical fragments in the m/z domain. Evaluation can
#be stopped at "generate_frag_formulae" or "generate_iso_dist" to retrieve
#monoisotopic and compsition or mass-domain isotopic information, respectively.
 

rawfile_name <- "YOUR AVERAGED, SINGLE-SCAN RAWFILE"
#Here, the program is expecting a Thermo rawfile containing a single, averaged MS/MS scan.
#Alternatively, centroids can be provided as a tibble of two columns named 
#"expr_mz" and "intensity

profile <- fetch_profile(rawfile_name)
#This extracts profile information from the rawfile, only used in plotting

centroids <- fetch_centroid(rawfile_name)
#This extracts centroid information from the rawfile, used in identifications

frag_IDs <- ID_zcon_frags(theor_frags, centroids, tolerance = 1.0e-5)
#This identifies fragments based on the m/z domain theoretical fragments and
#experimental centroid list

writexl::write_xlsx(frag_IDs, "YOUR FILENAME HERE")

nucleic_colors <- c("#4F81BD", "#4BACC6", "#C0504D", "#9BBB59", "#BFDF37", "#8064A2", "#8E81D8", "#2C4D75", "#772C2A", "#F79646", "#5F7530", "#4D3B62")
names(nucleic_colors) <- c("a", "a-B", "b", "c", "c-B", "d", "d-B", "w", "x", "x-B", "y", "z")
#Paired color scales are important.

ggplotly( 
  ggplot()+
    geom_col(data = centroids,
             aes(x = expr_mz, y = intensity), color = "lightgrey", width = 0.0001
    )+
    geom_line(data = profile,
              aes(x = expr_mz, y = intensity), color = "black"        
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
    scale_color_manual(values = nucleic_colors)+
    theme_publication()+
    guides(color = guide_legend(title = "Ion Type", nrow = 1))+
    # theme(legend.position = "none")+
    # xlim(900, 1000)+
    # ylim(0, .5)+
    labs(
      x = "m/z",
      y = "Intensity"
      # x = expression(bolditalic("m/z")),
      # y = expression(bold("Intensity"~"(10"^"5"*")"))
      # x = "",
      # y = ""
    )
)
#This plot generates an interactive spectrum with any annotations

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
ylim(-10, 0.25)+
  labs(
    color = "Ion Type"
  )+
  scale_shape_identity()+
  scale_color_manual(values = nucleic_colors)+
  theme_void(base_size = 14)+
  theme(
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size= unit(0.2, "cm"),
    legend.spacing = unit(0, "cm"),
    legend.title = element_text(face="italic"),
    plot.margin = unit(c(2,2,2,2),"mm")
  )+
  guides(color = guide_legend(nrow = 1))
#This generates a fragment ion map based on the sequence and identified fragments.

