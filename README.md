# Nucleic Acid Tools
---
This repository is for an R Shiny application designed for automated annotation of deconvoluted top-down nucleic acid fragmentation data

---
## Structure

The **Nucleic Acid Tools for Mass Spectrometry** application is designed for use with de-charged and de-isotoped mass spectrometry fragmentation data.  A single sequence can be input using the syntax outlined in the *Inputs* section to follow, which will automatically populate a parsed sequence indicative of the sugar, base, and backbone for each position and a corresponding list of theoretical fragments.  An experimental mass list can be uploaded to the application, where the expected format is a .xlsx file with two columns per sheet, each containing the deconvoluted mass (entitled "Mass") and intensity ("Intensity") from a single spectrum - any metadata present will be ignored, save for the name of the sheet.  Once uploaded, fragments will be idetified for each mass list provided, and sequence coverage maps, annotated spectra, and ppm error distributions will be generated.  Lists of identified fragments can be downloaded.

---
## Inputs

### Sequence

Sequence input follows a similar structure to IDT nomenclature, with some important differences.  The only *required* inputs are cannonical bases (ACGTU), and the interpreter will assume a deoxyribose and phosphate backbone present at these bases.  To specify a particular sugar or custom modification, a sugar-encoding symbol must be used *prior to* the base symbol.  For example, to specify a ribose sugar with a guanine base, the notation is: **rG**

Similarly, to specify a particular backbone, a backbone-encoding symbol must be used *after* the requisite base.  For example, to specify a thiophosphate on the **rG** from the previous example, the notation is:  **rg\***

Custom modifications can be specific for any sugar, base, or backbone by using parentheses directly *after* the appropriate symbol.  Modifications currently support a + or - followed by a number of arbitrary precision indicating the *net* gain or loss in mass at this site.  To specify a methylation of the ribose from the previous example, the notation is: **r(+14.014)G/***

More specifically for sugar modifications, however, a carbon can be indicated for a given modification, which is useful for the cross-ring cleavage fragment identification feature described in *Advanced Features* below. Notation in this case would be: **r(2'+14.014)G/***

Furthermore, methylribose can be easily indicated instead with the following notation: **mG/***

A list of encoded symbols is as follows:

Sugars:
    d = deoxyribose
    r = ribose
    m = methylribose
    + = 2'O - 4'C methyl bridged ribose (LNA)
    
Bases:
    A = adenine
    C = cytidine
    G = guanine
    T = thymine
    U = uracil
    
Backbone:
    p = phosphate
    /* = thiophosphate


Any 3' or 5' modifications can be specified in the side panel, below the sequence input.


### Deconvoluted Mass Lists

The application expects a .xlsx spreadsheet containing deconvoluted mass lists associated with the input sequence.  Each sheet in the .xlsx file can be named with the associated experimental conditions, and should contain one deconvoluted mass list from a given spectrum.  Within the deconvoluted mass list, the application will ignore any metadata present above the mass and intensity lists, but will look for a "Mass" and "Intensity" header in the left and right columns, respectively, to identify and collate the experimental data.

For each .xlsx sheet input, the theoretical fragment masses generated from the input sequence will be searched with the specified error tolerance provided in the side panel.  Based on the selected sheet name in the side panel, fragment identification maps, annotated spectra, and ppm errror distributions can be generated in their associated tabs.

---
## Outputs

### Parsed Sequence

The initial tab will generate a table of the parsed input sequence, showing the interpretation of any modifications specified in the sequence.  This tab is intended for proofreading prior to matching of theoretical fragment masses to experimenal data.

### Theoretical Fragments

Based on the parsed sequence, theoretical fragment masses will be generated.  This application considers all possible DNA terminal ion types (5' containing ions: a/b/c/d, 3' containing ions: w/x/y/z), alongside neutral loss of the base for a- and x-type ions.  All positions are relative to the terminus contained within the fragment ion.

### Identified Fragments

Once a sequence and experimental mass list are both input, the application will automatically identify theoretical fragment masses within the provided ppm tolerance with a default of 10 ppm.  In the case where multiple masses are within 10 ppm of a given theortical mass, the peak with a lower ppm error will be selected as a match.  All matched theoretical fragments are populated into a table containing data on the ion type, position, theoretical and experimental mass, and associated ppm error, which can be downloaded using the button in the side panel.  The output format of identified fragments follows the same convention as the input .xlsx spreadsheet, where the identified fragments from a given experimental mass list are output on a single sheet with the same name as the input mass lists.  The spreadsheet of identified fragments is also in .xlsx format.

### Plots

After processing fragment identifications for each input experimental mass list, the results can be visualized within the app in three ways: a fragment ion map over the sequence, an annotated spectrum, and a ppm error distribution.  Each of these plots is interactive using the **plotly** package.  Hovering over a point in each plot will display ion type, position, mass, and ppm error information for a given fragment, and each ion series can be toggled using the legend on the right side of each plot.  Using the drop-down menu in the side panel, fragment maps, annotated spectra, and ppm error distributions can be visualized for each experimental mass list provided, and their associated identified fragments.  Plots can be copied by selecting the camera icon when interacting with a plot.

---
## Advanced Features

### Cross-Ring Cleavages

Cross-ring cleavage refers to multiple cleavages of carbon-carbon bonds present in the sugar ring of each nucleoside.  This phenomenon was observed primarily in UVPD of glycan molecules, however, these cleavages can be observed in DNA spectra when  utilizing a high-energy fragmentation method.  To characterize these fragments, masses are calculated for each carbon in a standard 5-member ring present in ribose and deoxyribose, and sequentially subtracted from a- and z- type ions where the sugar ring is at the cleavage terminus for a fragment ion.  These types of fragments can provide in-depth localization of sugar modifications due to the propensity of UVPD and other high-energy fragmentation methods to preserve labile modifications.

### Internal Fragments

Despite the potential to yield in-depth characterization of a particular sequence, the vast complexity presented by searching for internal fragments combined with the extreme degree of overlap in the mass and m/z domain of these fragments precludes confident identification of any internal fragment without prior knowledge of the spectrum.  To try and mitigate the effect of the combinatorial nightmare that is internal fragments and corresponding frame-shifts, this application uses an "anchor" approach to constrain the internal fragment search space to only fragments that *might* be present in the spectrum.  All internal fragment masses are generated, however only those where a corresponding terminal fragment *or* its complement were identified in the spectrum prior to searching for internal fragments.  By searching for internal fragments in this way, only internal fragments where a terminal fragment from which it formed was identified are searched for, possibly lending more confidence to any identifications. Use this feature with an abundance of caution (and time, it takes a while).


