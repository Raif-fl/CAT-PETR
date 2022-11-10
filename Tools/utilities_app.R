#######################################################################
## compare_CyberT
##    dataframes = A list of dataframes that you wish to compare to each other 
##                  (Calculate fold change and log2 fold change and get P-value)
##    control = a vector of strings that correspond to the named elements of 
##              dataframes. The strings in this list will be compared to the 
##              strings in "treatment" on a one-to-one basis and will be treated 
##              as the control in all calculations (e.g. fold change is control/treatment)
##    treatment = a vector of strings that correspond to the named elements of 
##                dataframes. The strings in this list will be compared to the 
##                strings in "control" on a one-to-one basis and will be treated 
##                as the treatment in all calculations (e.g. fold change is control/treatment)
#######################################################################

compare_CyberT = function(dataframes, controls, treatments, analysis = NULL,
                          apo_name = NULL, no_conf = FALSE, var_norm = "vsn") {
  
  # For each dataframe, convert all column names to sentence case
  for (i in 1:length(dataframes)) {
    colnames(dataframes[[i]]) = stringr::str_to_sentence(colnames(dataframes[[i]]))
  }
  
  # Create empty lists which will contain the data from each sample comparison and the name of the comparison
  comparisons = list()
  comp_names = c()
  
  # Use the compare_CyberT function in combination with a loop to get P-values for all comparisons.
  for (i in 1:length(controls)) {
    
    print(paste("Processing Comparison ", toString(i), "/", toString(length(controls)), sep = ""))
   
    # define the control and treatment groups
    cont = dataframes[[controls[i]]]
    treat = dataframes[[treatments[i]]]
    # Remove "-" symbols as they cause trouble downstream. 
    cont$Name = stringr::str_replace_all(cont$Name, "-", "_")
    treat$Name = stringr::str_replace_all(treat$Name, "-", "_")
    # Add a fullname column based on the identifier columns. 
    cont$full_name = suppressWarnings(paste(cont$Name, cont$P_site, cont$Identifier))
    treat$full_name = suppressWarnings(paste(treat$Name, treat$P_site, treat$Identifier))
    
    # Detect if there are duplicates in the full name identifier and throw a warning. 
    if (sum(duplicated(cont$full_name)) > 0 | sum(duplicated(cont$full_name)) > 0) {
      warning(paste("Proteins/genes with identical names and identifiers detected. Suffixes will be
                    added to the name of duplicated Proteins/genes"))
      cont$full_name = make.unique(cont$full_name, sep = "_")
      treat$full_name = make.unique(treat$full_name, sep = "_")
    }
    
    # Merge the control and treatment groups by identifiers.
    suppressWarnings(if ("Identifier" %in% colnames(cont) & "P_site" %in% colnames(cont)) {
      df = dplyr::inner_join(cont, treat, by = c("full_name", "Name", "P_site", "Identifier"))
      check = T
    } else if ("Identifier" %in% colnames(cont)) {
      df = dplyr::inner_join(cont, treat, by = c("full_name", "Name", "Identifier"))
      check = F
    } else if ("P_site" %in% colnames(cont)) {
      df = dplyr::inner_join(cont, treat, by = c("full_name", "Name", "P_site"))
      check = T
    } else {
      df = dplyr::inner_join(cont, treat, by = c("full_name", "Name"))
      check = F
    })
    
    # Find the locations of everything that is a value. 
    cont_locs = which(stringr::str_detect(colnames(df), "\\.x"))
    treat_locs = which(stringr::str_detect(colnames(df), "\\.y"))
    # Change values column names into control and treatment column names. 
    for (j in cont_locs) {
      c_name = stringr::str_replace(colnames(df)[j], "Rep", "Control_Rep")
      c_name = stringr::str_remove(c_name, ".x")
      df = dplyr::rename(df, !!c_name := colnames(df)[j])
    }
    for (j in treat_locs) {
      c_name = stringr::str_replace(colnames(df)[j], "Rep", "Treatment_Rep")
      c_name = stringr::str_remove(c_name, ".y")
      df = dplyr::rename(df, !!c_name := colnames(df)[j])
    } 
    
    # If P_site is specified, then choose to look at either just apo or just phospho data. 
    suppressWarnings(if (!is.null(analysis) == T & check == T) {
      if (analysis == "Phospho-specific") {
        df = df[!grepl(apo_name, df$P_site),]
        rownames(df) = NULL  
      } else if (analysis == "Pan-specific") {
        df = df[grepl(apo_name, df$P_site),]
        rownames(df) = NULL  
      }
    })
    
    # Determine optimal CyberT window size 
    if (length(df$Name) < 50) {
      winSize = 3
    } else if (length(df$Name) < 1000) {
      winSize = 31
    } else if (length(df$Name) > 1000 & length(df$Name) < 2000) {
      winSize = 51
    } else if (length(df$Name) > 2000) {
      winSize = 101
    }
    
    # Determine optimal CyberT confidence
    if (ceiling(mean(length(cont_locs),length(treat_locs))) < 8) {
      conf = 9 - ceiling(mean(length(cont_locs),length(treat_locs)))
    } else {conf = 1}
    if (no_conf == TRUE) {
      conf = 0
    }
    
    # Alter the dataframe into an appropriate input for CyberT. 
    input = subset(df, select = c("full_name", colnames(df)[stringr::str_detect(colnames(df), "Rep")]))
    input = tibble::column_to_rownames(input, var = "full_name")
    
    # Apply a normalization
    if (var_norm == "vsn") {
      input = runVsn(input)
    } else if (var_norm == "logT") {
      input = log2(input)
    }
  
    # Run CyberT.
    results = bayesT(input, length(cont_locs), length(treat_locs), bayes = TRUE, winSize = winSize,
                                    doMulttest=TRUE, conf = conf)
    
    
    # Bind the adjusted P-values from CyberT to the dataframe.
    df$P_value_adjust = results$BH
    
    # Calculate fold changes and the -log10 P-values.
    df$Treated_Average = rowMeans(df[stringr::str_detect(colnames(df), "Treatment")])
    df$Control_Average = rowMeans(df[stringr::str_detect(colnames(df), "Control")])
    df$Fold_change = df$Treated_Average/df$Control_Average
    df$log2_FC = log2(df$Fold_change)
    df$log_p = -log10(df$P_value_adjust)
  
    # append to lists
    comparisons = append(comparisons, list(df))
    comp_names = append(comp_names, paste(controls[i], "vs", treatments[i], sep = "_"))

  }
  
  # Use the comparison names to name the elements of the list of comparisons. 
  names(comparisons) = comp_names
  
  return(comparisons)
  
}

#######################################################################
## volcano_plot_app
##    data = a dataframe. 
##    to_label = a list of names that specify which proteins/genes to label on the 
##               volcano plot in addition to those chosen via top. 
##    top = An integer which specifies how many of the top proteins to label. Every protein 
##          will be sorted by a combination of its log fold change and P-value. The proteins 
##          with the highest value for this combination (top proteins) will be selected from this. 
##    FC_range = Names with log2 fold changes outside of this range may be considered 
##               significant. 
##    P_cutoff = Names with BH adjusted P-values below this cutoff as determined by 
##               CyberT may be considered significant. 
##    mycolors = a list of 4 colors which will be used for the significantly downregulated,
##               significantly upregulated, non-significantly affected, and labeled
##               proteins respectively.
##    text_size = size of text used for labels.
##    axes_text_size = size of text used for the tick marks on the axes.
##    axes_label_size = size of text used for the axes labels.
##    point_size = A vector which contains the sizes of the non-significant and significant
##                 data points respectively. 
##    range = the range shown on the X-axis. 
##    box_pad = amount of padding around the label boxes. 
##    point_pad = amount of padding around the labelled data points. 
##    label_options = The parts of the label that will be included. 
#######################################################################

volcano_plot_app = function (data, to_label = c(), top = 0, FC_range = c(-1,1), P_cutoff = 3,
                             mycolors = c("blue", "red", "black", "green"), text_size = 3,
                             axes_text_size = 2, axes_label_size = 3,  point_sizes = c(1, 2),
                             range = NULL, box_pad = 0.3, point_pad = 0.2, 
                             label_options = c("Name", "P_site", "Identifier"),
                             sig_label = FALSE) {
  
  # Remove infinite and NA values. 
  data = data[is.finite(rowSums(data[colnames(data) == "log2_FC" | colnames(data) == "log_p"])),]
  
  # Define a new column based on what form the labels will take. 
  data$label_form = do.call(paste, c(data[label_options], sep=" "))
  
  # prevent the function from breaking in the presence of null values for top
  if (is.na(top)) {top = 0}
  
  # Create a column which specifies if something is above our cutoffs or of special interest.  
  data$reg = "Non-significant"
  data$reg[data$log2_FC > FC_range[2] & data$log_p > P_cutoff] = "Up-regulated"
  data$reg[data$log2_FC < FC_range[1] & data$log_p > P_cutoff] = "Down-regulated"
  data$reg[data$full_name %in% to_label] = "Additional"
  
  # If sig labels is true, then make top equal to the number of significant genes/proteins
  if (sig_label) {top = length(data$reg[data$reg == "Up-regulated"]) + length(data$reg[data$reg == "Down-regulated"])}
  
  # Specify the colors the correspond to up, down, not significant, and additional to be labeled. 
  names(mycolors) = c("Down-regulated", "Up-regulated", "Non-significant", "Additional")
  
  # Add labels for the proteins of interest.  
  data$delabel = NA
  data$delabel[data$full_name %in% to_label] = data$label_form[data$full_name %in% to_label]
  
  # Create a column with a combined metric for FC and P then choose a set number to label.
  if (top > 0) {
    data$dist = sqrt(data$log_p^2 + data$log2_FC^2)
    data$dist[data$reg == "Non-significant"] = NA
    data = data[order(data$dist, decreasing = T),]
    top_reg = data[is.na(data$dist) == F,][1:top,]$full_name
    data$delabel[data$full_name %in% top_reg] = data$label_form[data$full_name %in% top_reg] 
  } else {top_reg = c()}
  
  # Check to see if the maximum labels have been overcome. 
  if ((length(unique(c(to_label, top_reg))) > 150)) {stop("Exceeded maximum number of labels (150)")}
  
  # See if there are any dupliate labels and add identifier to fix them. 
  dupl_loc = duplicated(data$delabel, incomparables=NA) | duplicated(data$delabel, fromLast = T, incomparables=NA)
  data$delabel[dupl_loc] = data$full_name[dupl_loc]
  
  # Create different point sizes for labelled or colored points. 
  data$size_p = F
  data$size_p[is.na(data$delabel) == F] = T
  data$size_p[data$reg != "Non-significant"] = T
  
  # Set the default range as the symmetrical maximum log2 FC. 
  if (is.null(range)) {
    up_lim = max(data$log2_FC, na.rm = T)
    down_lim = min(data$log2_FC, na.rm = T)
    max_lim = max(c(abs(up_lim), abs(down_lim)))
    range = c(-max_lim-1, max_lim+1)
  }

  # Create the base scatter plot.
  base_plot = ggplot(data=map_df(arrange(data, delabel), rev),
                     aes(x=log2_FC, y=log_p, col=reg, size = size_p)) +
    geom_point() + 
    theme_bw() +
    geom_vline(xintercept=FC_range, col="red") +
    geom_hline(yintercept=P_cutoff, col="red") +
    scale_color_manual(name = "Legend", values = mycolors) + 
    ylab("-Log 10 P value") + ggtitle("") +
    xlab("Log2 Fold Change") + ggtitle("") +
    xlim(range[1], range[2]) + 
    scale_size_manual(values=point_sizes, guide = "none") +
    theme(legend.position="none", axis.text = element_text(size = axes_text_size),
          axis.title = element_text(size = axes_label_size))
  
  # Add the labels.
  plot = base_plot + geom_label_repel(aes(label = delabel),
                       size = text_size,
                       label.size = NA,
                       min.segment.length = 0.000001,
                       max.overlaps = 10000,
                       box.padding   = box_pad, 
                       point.padding = point_pad,
                       na.rm = T,
                       force = 2,
                       segment.size = 0.3,
                       segment.color = 'grey50',
                       arrow = arrow(length = unit(0.07, "inches")),
                       show.legend = FALSE) +
                       scale_y_continuous(expand = expansion(mult = 0.1))
  
  return(list(plot, unique(top_reg[!is.na(top_reg)])))
}

#######################################################################
## one_to_one_app
##    data: A dataframe that contains the log fold changes from two comparisons. 
##    comp1 + comp2: strings that specify the column names in data that contain 
##                   the log fold changes.
##    to_label: A list of gene/protein names that will be labelled on the scatter plot in addition 
##              to those specified by top. 
##    top = An integer which specifies how many of the top proteins to label. Every protein 
##          will be sorted by the difference between its log fold change to one. The proteins 
##          with the highest difference (top proteins) will be selected from this. 
##    quant_int: The interval that will be used for calculating the quantiles from the data. 
##    line_color: The color of the one to one ratio line that will be used. 
##    int_color: The color that will be used for the quantile lines. 
##    axes_text_size = size of text used for the tick marks on the axes.
##    axes_label_size = size of text used for the axes labels.
##    text_size: The size of the text labels. 
##    point_size: The size to be used for the data points. 
##    box_pad = amount of padding around the label boxes. 
##    point_pad = amount of padding around the labelled datapoint. 
##    mycolors = The colors that will be used for points within, above, and below quants along 
##               with the color used for additional labels specified in to_label. 
##    label_options = The parts of the label that will be included. 
#######################################################################

one_to_one_app = function (data, comp1, comp2, to_label = c(), top = 0, quant_int = c(0.01, 0.99),
                           line_color = "blue", int_color = "red", text_size = 2.3,  point_size = 1,
                           axes_text_size = 2, axes_label_size = 3, box_pad = 0.3, point_pad = 0.2,
                           mycolors = c("black", "black", "black", "green"), sig_label = FALSE,
                           label_options = c("Name", "P_site", "Identifier")) {
  
  # Ensure that only valid names are used. 
  comp1 = make.names(comp1)
  comp2 = make.names(comp2)
  
  # Remove infinite and NA values. 
  data = drop_na(data)
  data = data[is.finite(rowSums(data[colnames(data) == comp1 | colnames(data) == comp2])),]
  
  # Define a new column based on what form the labels will take. 
  data$label_form = do.call(paste, c(data[label_options], sep=" "))
  
  # prevent the function from breaking in the presence of null values for top
  if (is.null(top) | is.na(top) | is.nan(top)) {top = 0}
  
  # Find the quantiles from the data. 
  quants = quantile(data$diff, probs = quant_int, na.rm = TRUE)
  
  # Create a row for the labels.
  data$delabel = NA
  
  # Create a column which specifies if something is outside of the quantiles or of special interest.  
  data$reg = "within-quants"
  data$reg[data$diff > quants[2]] = "above-quants"
  data$reg[data$diff < quants[1]] = "below-quants"
  data$reg[data$full_name %in% to_label] = "Additional"
  
  # If sig labels is true, then make top equal to the number of significant genes/proteins
  if (sig_label) {top = length(data$reg[data$reg == "above-quants"]) + length(data$reg[data$reg == "below-quants"])}
  
  
  # Specify the colors the correspond to within quants, outside quants, and of special interest.
  names(mycolors) = c("within-quants", "above-quants", "below-quants", "Additional")
  
  # Label the top proteins with the highest difference. 
  if (top > 0) {
    data$diff2 = data$diff
    data$diff2[data$reg == "within-quants"] = NA
    data = data[order(abs(data$diff2), decreasing = T),]
    top_reg = data[is.na(data$diff2) == F,][1:top,]$full_name
    data$delabel[data$full_name %in% top_reg] = data$label_form[data$full_name %in% top_reg]
  } else {top_reg = c()}
  
  # Check to see if the maximum labels have been overcome. 
  if ((length(unique(c(to_label, top_reg))) > 150)) {stop("Error: Exceeded maximum number of labels (150)")}
  
  # Add the labels from to_label.  
  data$delabel[data$full_name %in% to_label] = data$label_form[data$full_name %in% to_label]
  
  # See if there are any dupliate labels and add identifier names to fix them. 
  dupl_loc = duplicated(data$delabel, incomparables=NA) | duplicated(data$delabel, fromLast = T, incomparables=NA)
  data$delabel[dupl_loc] = data$full_name[dupl_loc]
  
  # Create a plot and then plot a one-to-one ratio line. 
  base_plot = ggplot(data, aes_string(x = comp1, y = comp2, col = "reg")) +
    geom_point(size = point_size) +
    geom_abline(slope = 1, color = line_color, size = 1.5) +
    geom_abline(slope = 1, intercept = (-1*quants[2]), color = int_color, linetype = "dashed") +
    geom_abline(slope = 1, intercept = (-1*quants[1]), color = int_color, linetype = "dashed") +
    scale_color_manual(name = "Legend", values = mycolors) + 
    theme_bw() +
    ylab(comp2) +
    xlab(comp1) + 
    theme(legend.position="none", axis.text = element_text(size = axes_text_size),
          axis.title = element_text(size = axes_label_size))
  
  # find the maximum values for the X and Y axes. 
  y_max = abs_max(data[[comp2]])
  x_max = abs_max(data[[comp1]])
  
  # Add the labels. 
  label_plot = base_plot + geom_label_repel(aes(label = delabel),
                                              size = text_size,
                                              min.segment.length = 0.000001,
                                              label.size = NA,
                                              max.overlaps = 10000,
                                              box.padding   = box_pad, 
                                              point.padding = point_pad,
                                              na.rm = T,
                                              segment.color = 'grey50',
                                              segment.size = 0.3,
                                              arrow = arrow(length = unit(0.07, "inches")),
                                              show.legend = FALSE) + 
      coord_cartesian(xlim = c(-x_max-1, x_max+1),
                      ylim = c(-y_max-1, y_max+1))
  
  return(list(label_plot, unique(top_reg[!is.na(top_reg)])))
}

#######################################################################
## abs_max. A helper function which finds the maximum value in a vector regardless of sign. 
#######################################################################
abs_max = function(vector) {
  # Define the upper and lower limits based on the maximum and minimum. 
  up_lim = max(vector, na.rm = T)
  down_lim = min(vector, na.rm = T)
  max_lim = max(c(abs(up_lim), abs(down_lim)))
  return(max_lim)
}

#######################################################################
##  hmap_prep
##    dataframes = A list of dataframes outputted by the compare_CyberT function
##    name_searchs = A vector containing the names of all names to plot as strings. 
##    title = A string that will be used as the title of the plot. 
##    values = A string which specifies which column you want to use for the 
##             heatmap values. 
##    order = A list of strings that correspond to the named elements in the 
##            dataframes list. The order of this list will become the order of 
##            the columns in the heatmap. 
#######################################################################

hmap_prep = function (dataframes, title = "", order = c(),
                      label_options = c("Name", "P_site", "Identifier")) {
  
  # add the name of the comparison to a column on the respective dataframe
  for (i in 1:length(dataframes)) {
    dataframes[[i]]$comp = names(dataframes[i])
  }
  
  # Take all of the dataframes and bind them together.
  bound = bind_rows(dataframes)
  
  # Define a new column which will be the names of our matrix. 
  bound$Name = do.call(paste, c(bound[label_options], sep=" "))
  
  # subset the bound data in order to focus on the meaningful columns. 
  bound = subset(bound, select = c("Name", "comp", "log2_FC", "full_name"))
  
  # rearrange the order of the columns
  if (length(order) > 0) {
    bound$comp = factor(x = bound$comp, levels = order)
  }
  
  return(bound)
}

#######################################################################
## hmap 
##    bound = a dataframe that was created by the hmap_prep() function.
##    name_search = A vector containing the names of all genes/proteins to plot as strings. 
##    sort_by = An integer which specifies which of the columns the heatmap will sort by (left to right).
##    heat_comps = A list of strings which specify which of the comparisons contained in bound will
##                 make it into the final product. 
##    heat_num = A number which specifies the upper and lower bounds used for the colorbar.
##    height_hmap = The height that will be added to the plot for each new gene/protein being visualized. 
##    text_size = the size of the text used for both the X and Y axes. 
##    lg_title_size = the size of the text used for the legend title. 
##    lg_text_size = the size of the text used for the tick marks on the legend. 
#######################################################################

hmap = function(bound, name_search, sort_by, heat_comps, heat_num, height_hmap = 30,
                text_size = 12, lg_title_size = 11, lg_text_size = 10, color_choice = "RdBu",
                reverse_scale = TRUE){
  
  # Convert to a matrix. 
  bound_mat = data.matrix(bound)
  
  # Create a duplicate dataframe which will be used for the hover information. 
  bound_mat2 = bound_mat
  bound_mat2 = round(bound_mat2, digits = 4)
  
  # Create a matrix where everything is set to right below the maximum values defined by the
  # user. 
  bound_mat[bound_mat > heat_num[2]] = heat_num[2] - 0.001
  bound_mat[bound_mat < heat_num[1]] = heat_num[1] + 0.001
  
  # Get the length of the matrix along the X-axis
  mat_length = length(bound_mat[,1])
  
  # Create the heatmap. 
  heatmap = plot_ly(colors = color_choice) %>%
    add_heatmap(x = colnames(bound_mat), y = rownames(bound_mat), z = bound_mat, reversescale = reverse_scale,
                text = bound_mat2, zmin = heat_num[1], zmax = heat_num[2],
                hovertemplate = paste('Name: %{y}<extra></extra><br>',
                                      'Comparison: %{x}<br>',
                                      'Log2 FC: %{text}'),
                colorbar = list(limits = c(heat_num[1], heat_num[2]),
                                len = (120 + (mat_length*height_hmap)/2), lenmode = "pixels",
                                title = list(text = "log2 FC", font = list(size = lg_title_size)),
                                tickfont = list(size = lg_text_size), yanchor = "middle")) %>%
    layout(
      xaxis = list(tickfont = list(size = text_size)),
      yaxis = list(tickfont = list(size = text_size))) %>%
    config(modeBarButtons = list(list("toImage")),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE, displayModeBar = TRUE)
  
  return(heatmap)
}

########################################################################################
## clean kinexus
##    path = The path to a raw data file provided by Kinexus.
##    col_names_row = The row which contains a set of names to be converted
##                    Into the column names of the dataframe
##    max_error = The maximum percent error that is tolerated between technical replicates
##                where percent error is calculated as:
##                avg(replicate intensity - mean intensity)/avg(intensity) * 100
#########################################################################################

clean_kinexus = function(path, col_names_row = 53, max_error = 50) {
  
  # Read in the raw kinexus KAM-1325 data file. 
  raw_df = suppressMessages(read_tsv(path))
  
  # Make a certain row the new column name and then subsequently delete that row.  
  colnames(raw_df) = raw_df[col_names_row,]
  raw_df = raw_df[-col_names_row,]
  colnames(raw_df)
  
  # select the columns of interest. 
  raw_df = subset(raw_df, select = c("Target Name with alias", "Human P-Site", "Cat. No.",
                                     "Signal Median", "Background Median", "Signal Area", "Flag"))
  
  # rename the columns
  colnames(raw_df) = c("Target_Name", "P_Site", "Antibody",
                       "Signal_Median", "Background_Median", "Spot_Area", "Flag")
  
  # remove the na containing rows. 
  raw_df = drop_na(raw_df)
  
  # convert to numeric 
  raw_df = mutate(raw_df, Signal_Median = as.numeric(Signal_Median),
                  Background_Median = as.numeric(Background_Median),
                  Spot_Area = as.numeric(Spot_Area),
                  Flag = as.numeric(Flag))
  
  # calulate the raw intensity. 
  raw_df = mutate(raw_df, raw_intensity = (Signal_Median - Background_Median)*(Spot_Area/100))
  
  # Calculate the scalar needed for normalization
  summed_raw_intensity = sum(raw_df$raw_intensity)
  scalar = 20000000/summed_raw_intensity
  
  # calculate the normalized intensity. 
  raw_df = mutate(raw_df, Normalized_intensity = raw_intensity*scalar)
  
  # Remove antibodies with flags equal to 1. 
  flagged = filter(raw_df, Flag == 1 | Flag == 2)$"Antibody"
  raw_df = filter(raw_df, !raw_df$"Antibody" %in% flagged)
  
  # using group by to calculate the average Intensity. 
  df_group = group_by(raw_df, Antibody)
  df = suppressMessages(summarise(df_group, Target_Name = Target_Name, P_Site = P_Site,
                                  Average_intensity = mean(Normalized_intensity),
                                  Normalized_intensity = Normalized_intensity))
  
  # Calculate the error ranges
  df$diff = abs(df$Average_intensity - df$Normalized_intensity)
  df = suppressMessages(summarise(df, Target_Name = Target_Name, P_Site = P_Site,
                                  Average_intensity = Average_intensity,
                                  Error_range = mean(diff)))
  
  # remove duplicate rows 
  df = distinct(df, Antibody, .keep_all = TRUE)
  
  # Calculate the percent error. 
  df$percent_error = (df$Error_range/df$Average_intensity)*100
  
  # Remove antibodies with % error above a certain cutoff. 
  df = filter(df, percent_error < max_error)
  
  # Remove everything in brackets from the protein names.
  df$Target_Name = gsub(r"{\s*\([^\)]+\)}","",as.character(df$Target_Name))
  
  # Sort by target name and rename columns. 
  df = df[order(df$Target_Name),]
  df = df[, c("Target_Name", "P_Site", "Antibody", "Average_intensity", "Error_range", "percent_error")]
  colnames(df) = c("Name", "P_Site", "Identifier", "Rep_1", "Error_range", "percent_error")
  
  # Select only the columns of interest
  df = subset(df, select = c("Name", "P_Site", "Identifier", "Rep_1"))
  
  return(df)
}

#######################################################################
## abs_max. A helper function which allows for significantly more flexible rounding. 
#######################################################################
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

