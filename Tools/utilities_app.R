## compare cybert2 #####

compare_CyberT2 = function(infile, controls, treatments, analysis = NULL,
                          apo_name = NULL, no_conf = FALSE, var_norm = "vsn",
                          kd = FALSE, max_error = 50){
  
  if (kd == TRUE) {
    datatables = clean_kinexus(infile, max_error = max_error)
    names(datatables) = stringr::str_remove(infile$name, ".txt")
  } else {
    datatables = purrr::map(infile$datapath, data.table::fread)
    names(datatables) = stringr::str_remove(infile$name, ".csv")
  }
  
  # Create a copy of all of the column options used. 
  col_options = stats::na.omit(stringr::str_extract(stringr::str_to_sentence(colnames(datatables[[1]])),
                                                    pattern = "Name|P_site|Identifier"))
  
  # Clean up the contents of the dataframe.
  for (i in 1:length(datatables)) {
    # Adjust the format of the column names. 
    colnames(datatables[[i]]) = stringr::str_to_sentence(colnames(datatables[[i]]))
    
    # Replace Rep with a unique identifier for each experiment. 
    colnames(datatables[[i]]) = stringr::str_replace(colnames(datatables[[i]]), "Rep", names(datatables[i])[1])
    
    # Remove - signs as they cause problems later. 
    datatables[[i]][, Name := stringr::str_replace_all(Name, "-", "_")]
    
    # If P_site is specified, then choose to look at either just apo or just phospho data. 
    suppressWarnings(if (!is.null(analysis) == T & "P_site" %in% colnames(datatables[[i]])) {
      if (analysis == "Phospho-specific") {
        datatables[[i]] = datatables[[i]][ datatables[[i]][, .I[!grepl(apo_name, P_site)], by = P_site]$V1 ]
      } else if (analysis == "Pan-specific") {
        datatables[[i]] = datatables[[i]][ datatables[[i]][, .I[grepl(apo_name, P_site)], by = P_site]$V1 ]
      }
    })
    
    # Create the full name column based on the identifiers. 
    suppressWarnings(if ("Identifier" %in% colnames(datatables[[i]]) & "P_site" %in% colnames(datatables[[i]])) {
      datatables[[i]][, full_name:=paste(Name, P_site, Identifier, sep = " _@_ ")]
    } else if ("Identifier" %in% base::colnames(datatables[[i]])) {
      datatables[[i]][, full_name:=paste(Name, Identifier, sep = " _@_ ")]
    } else if ("P_site" %in% base::colnames(datatables[[i]])) {
      datatables[[i]][, full_name:=paste(Name, P_site, sep = " _@_ ")]
    } else {
      datatables[[i]][, full_name:=Name]
    })
    
    # Remove the Name, P_site, and Identifier columns to save on space. 
    suppressWarnings(datatables[[i]][, c("Name", "P_site", "Identifier") := NULL])
    
    # check if duplicates are present and if so throw a warning 
    if (sum(duplicated(datatables[[i]]$full_name)) > 0) {
      warning(paste("Proteins/genes with identical names and identifiers detected. Suffixes will be
                    added to the name of duplicated Proteins/genes"))
      datatables[[i]]$full_name = make.unique(datatables[[i]]$full_name, sep = "_")
    }
  }
  
  comparisons = list()
  comp_names = c()
  
  for (i in 1:length(controls)) {
    print(base::paste("Processing Comparison ", toString(i), "/", toString(length(controls)), sep = ""))
    
    # Find the names to be used by cybert for controls and treats. 
    numC = sum(!colnames(datatables[[controls[i]]]) == "full_name")
    numE = sum(!colnames(datatables[[treatments[i]]]) == "full_name")
    
    print
    
    # Set the keys for the merge. 
    data.table::setkey(datatables[[controls[i]]], full_name)
    data.table::setkey(datatables[[treatments[i]]], full_name)
    
    # Perform the merge. 
    dt = datatables[[controls[i]]][datatables[[treatments[i]]], nomatch = 0]
    
    # Save full names for later then remove it. 
    full = dt$full_name
    dt[, full_name := NULL]
    
    # Determine optimal CyberT window size 
    if (nrow(dt) < 50) {
      winSize = 3
    } else if (nrow(dt) < 1000) {
      winSize = 31
    } else if (nrow(dt) > 1000 & length(dt$full_name) < 2000) {
      winSize = 51
    } else if (nrow(dt) > 2000) {
      winSize = 101
    }
    
    # Determine optimal CyberT confidence
    if (ncol(dt)/2 < 8) {
      conf = 9 - ncol(dt)/2
    } else {conf = 1}
    if (no_conf == TRUE) {
      conf = 0
    }
    
    # Apply a normalization
    if (var_norm == "vsn") {
      input = runVsn(dt)
      input = data.table::as.data.table(dt)
    } else if (var_norm == "logT") {
      input = log(dt)
      input = data.table::as.data.table(dt)
    }
    
    # Run CyberT.
    results = bayesT(dt, numC, numE, bayes = TRUE, winSize = winSize,
                     doMulttest=TRUE, conf = conf)
    
    # Extract just the results and convert to data.table
    results = subset(results, select = ("runMulttest.pVal."))
    
    # Create the finished data table. 
    dt[, full_name := stringr::str_replace_all(full, " _@_ ", " ")]
    dt[, c(col_options) := data.table::tstrsplit(full, " _@_ ", fixed = TRUE)]
    dt[, Control_Average := rowMeans(dt[,1:numC])]
    dt[, Treated_Average := rowMeans(dt[,(numC+1):(numC+numE)])]
    dt[, c(1,2) := NULL]
    dt[, Fold_change := Treated_Average/Control_Average]
    dt[, log2_FC := log2(Fold_change)]
    dt[, P_value_adjust := results$runMulttest.pVal.]
    dt[, log_p := -log10(P_value_adjust)]
    
    # append to lists
    comparisons = base::append(comparisons, list(dt))
    comp_names = base::append(comp_names, base::paste(controls[i], "vs", treatments[i], sep = "_"))
    
    # Determine if old datatables can be deleted. 
    if (!names(datatables[controls[i]]) %in% controls[i+1:length(controls)] &&
        names(datatables[controls[i]]) %in% treatments[i+1:length(treatments)]) {
      datatables[controls[i]] = NULL
    }
    if (!names(datatables[treatments[i]]) %in% controls[i+1:length(controls)] &&
        !names(datatables[treatments[i]]) %in% treatments[i+1:length(treatments)]) {
      datatables[treatments[i]] = NULL
    }
  }
  
  # Use the comparison names to name the elements of the list of comparisons. 
  names(comparisons) = comp_names
  return(comparisons)
}



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
    base::colnames(dataframes[[i]]) = stringr::str_to_sentence(base::colnames(dataframes[[i]]))
  }
  
  # Create empty lists which will contain the data from each sample comparison and the name of the comparison
  comparisons = list()
  comp_names = c()
  
  # Use the compare_CyberT function in combination with a loop to get P-values for all comparisons.
  for (i in 1:length(controls)) {
    
    print(base::paste("Processing Comparison ", toString(i), "/", toString(length(controls)), sep = ""))
   
    # define the control and treatment groups
    cont = dataframes[[controls[i]]]
    treat = dataframes[[treatments[i]]]
    # Remove "-" symbols as they cause trouble downstream. 
    cont$Name = stringr::str_replace_all(cont$Name, "-", "_")
    treat$Name = stringr::str_replace_all(treat$Name, "-", "_")
    # Add a fullname column based on the identifier columns. 
    cont$full_name = suppressWarnings(base::paste(cont$Name, cont$P_site, cont$Identifier))
    treat$full_name = suppressWarnings(base::paste(treat$Name, treat$P_site, treat$Identifier))
    
    # Detect if there are duplicates in the full name identifier and throw a warning. 
    if (sum(base::duplicated(cont$full_name)) > 0 | sum(base::duplicated(cont$full_name)) > 0) {
      warning(base::paste("Proteins/genes with identical names and identifiers detected. Suffixes will be
                    added to the name of duplicated Proteins/genes"))
      cont$full_name = make.unique(cont$full_name, sep = "_")
      treat$full_name = make.unique(treat$full_name, sep = "_")
    }
    
    # Merge the control and treatment groups by identifiers.
    suppressWarnings(if ("Identifier" %in% base::colnames(cont) & "P_site" %in% base::colnames(cont)) {
      df1 = dplyr::inner_join(cont, treat, by = c("full_name", "Name", "P_site", "Identifier"))
      check = T
    } else if ("Identifier" %in% base::colnames(cont)) {
      df1 = dplyr::inner_join(cont, treat, by = c("full_name", "Name", "Identifier"))
      check = F
    } else if ("P_site" %in% base::colnames(cont)) {
      df1 = dplyr::inner_join(cont, treat, by = c("full_name", "Name", "P_site"))
      check = T
    } else {
      df1 = dplyr::inner_join(cont, treat, by = c("full_name", "Name"))
      check = F
    })
    
    # Find the locations of everything that is a value. 
    cont_locs = base::which(stringr::str_detect(base::colnames(df1), "\\.x"))
    treat_locs = base::which(stringr::str_detect(base::colnames(df1), "\\.y"))
    # Change values column names into control and treatment column names. 
    for (j in cont_locs) {
      c_name = stringr::str_replace(base::colnames(df1)[j], "Rep", "Control_Rep")
      c_name = stringr::str_remove(c_name, ".x")
      df1 = dplyr::rename(df1, !!c_name := base::colnames(df1)[j])
    }
    for (j in treat_locs) {
      c_name = stringr::str_replace(base::colnames(df1)[j], "Rep", "Treatment_Rep")
      c_name = stringr::str_remove(c_name, ".y")
      df1 = dplyr::rename(df1, !!c_name := base::colnames(df1)[j])
    } 
    
    # If P_site is specified, then choose to look at either just apo or just phospho data. 
    suppressWarnings(if (!is.null(analysis) == T & check == T) {
      if (analysis == "Phospho-specific") {
        df1 = df1[!base::grepl(apo_name, df1$P_site),]
        base::rownames(df1) = NULL  
      } else if (analysis == "Pan-specific") {
        df1 = df1[base::grepl(apo_name, df1$P_site),]
        base::rownames(df1) = NULL  
      }
    })
    
    # Determine optimal CyberT window size 
    if (length(df1$Name) < 50) {
      winSize = 3
    } else if (length(df1$Name) < 1000) {
      winSize = 31
    } else if (length(df1$Name) > 1000 & length(df1$Name) < 2000) {
      winSize = 51
    } else if (length(df1$Name) > 2000) {
      winSize = 101
    }
    
    # Determine optimal CyberT confidence
    if (ceiling(base::mean(length(cont_locs),length(treat_locs))) < 8) {
      conf = 9 - ceiling(base::mean(length(cont_locs),length(treat_locs)))
    } else {conf = 1}
    if (no_conf == TRUE) {
      conf = 0
    }
    
    # Alter the dataframe into an appropriate input for CyberT. 
    input = base::subset(df1, select = c("full_name", base::colnames(df1)[stringr::str_detect(base::colnames(df1), "Rep")]))
    input = tibble::column_to_rownames(input, var = "full_name")
    
    # Apply a normalization
    if (var_norm == "vsn") {
      input = runVsn(input)
    } else if (var_norm == "logT") {
      input = log(input)
    }
  
    # Run CyberT.
    results = bayesT(input, length(cont_locs), length(treat_locs), bayes = TRUE, winSize = winSize,
                                    doMulttest=TRUE, conf = conf)
    
    
    # Bind the adjusted P-values from CyberT to the dataframe.
    df1$P_value_adjust = results$BH
    
    # Calculate fold changes and the -log10 P-values.
    df1$Treated_Average = base::rowMeans(df1[stringr::str_detect(base::colnames(df1), "Treatment")])
    df1$Control_Average = base::rowMeans(df1[stringr::str_detect(base::colnames(df1), "Control")])
    df1$Fold_change = df1$Treated_Average/df1$Control_Average
    df1$log2_FC = log2(df1$Fold_change)
    df1$log_p = -log10(df1$P_value_adjust)
  
    # append to lists
    comparisons = base::append(comparisons, list(df1))
    comp_names = base::append(comp_names, base::paste(controls[i], "vs", treatments[i], sep = "_"))

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
##    range = The range shown on the X-axis. 
##    box_pad = amount of padding around the label boxes. 
##    point_pad = amount of padding around the labelled data points. 
##    label_options = The parts of the label that will be included. 
##    sig_label = A boolean which specifies whether all genes/proteins outside of cutoff 
##                should be labelled. 
#######################################################################

volcano_plot_app = function (data, to_label = c(), top = 0, FC_range = c(-1,1), P_cutoff = 3,
                             mycolors = c("blue", "red", "black", "green"), text_size = 3,
                             axes_text_size = 2, axes_label_size = 3,  point_sizes = c(1, 2),
                             range = NULL, box_pad = 0.3, point_pad = 0.2, 
                             label_options = c("Name", "P_site", "Identifier"),
                             sig_label = FALSE) {
  
  data = as.data.frame(data)
  
  # Remove infinite and NA values. 
  data = data[is.finite(rowSums(data[base::colnames(data) == "log2_FC" | base::colnames(data) == "log_p"])),]
  
  # Define a new column based on what form the labels will take. 
  data$label_form = do.call(base::paste, c(data[label_options], sep=" "))
  
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
  dupl_loc = base::duplicated(data$delabel, incomparables=NA) | base::duplicated(data$delabel, fromLast = T, incomparables=NA)
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
  base_plot = ggplot2::ggplot(data=purrr::map_df(dplyr::arrange(data, delabel), rev),
                     ggplot2::aes(x=log2_FC, y=log_p, col=reg, size = size_p)) +
    ggplot2::geom_point() + 
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept=FC_range, col="red") +
    ggplot2::geom_hline(yintercept=P_cutoff, col="red") +
    ggplot2::scale_color_manual(name = "Legend", values = mycolors) + 
    ggplot2::ylab("-Log 10 P value") + ggplot2::ggtitle("") +
    ggplot2::xlab("Log2 Fold Change") + ggplot2::ggtitle("") +
    ggplot2::xlim(range[1], range[2]) + 
    ggplot2::scale_size_manual(values=point_sizes, guide = "none") +
    ggplot2::theme(legend.position="none", axis.text = ggplot2::element_text(size = axes_text_size),
          axis.title = ggplot2::element_text(size = axes_label_size))
  
  # Add the labels.
  plot = base_plot + ggrepel::geom_label_repel(ggplot2::aes(label = delabel),
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
                       arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches")),
                       show.legend = FALSE) +
                       ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.1))
  
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
##    sig_label = A boolean which specifies whether all genes/proteins outside of the quantiles 
##                should be labelled. 
##    label_options = The parts of the label that will be included. 
#######################################################################

one_to_one_app = function (data, comp1, comp2, to_label = c(), top = 0, quant_int = c(0.01, 0.99),
                           line_color = "blue", int_color = "red", text_size = 2.3,  point_size = 1,
                           axes_text_size = 2, axes_label_size = 3, box_pad = 0.3, point_pad = 0.2,
                           mycolors = c("black", "black", "black", "green"), sig_label = FALSE,
                           label_options = c("Name", "P_site", "Identifier")) {
  
  data = as.data.frame(data)
  
  # Ensure that only valid names are used. 
  comp1 = make.names(comp1)
  comp2 = make.names(comp2)
  
  # Remove infinite and NA values. 
  data = tidyr::drop_na(data)
  data = data[is.finite(base::rowSums(data[base::colnames(data) == comp1 | base::colnames(data) == comp2])),]
  
  # Define a new column based on what form the labels will take. 
  data$label_form = do.call(base::paste, c(data[label_options], sep=" "))
  
  # prevent the function from breaking in the presence of null values for top
  if (is.null(top) | is.na(top) | is.nan(top)) {top = 0}
  
  # Find the quantiles from the data. 
  quants = stats::quantile(data$diff, probs = quant_int, na.rm = TRUE)
  
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
  dupl_loc = base::duplicated(data$delabel, incomparables=NA) | base::duplicated(data$delabel, fromLast = T, incomparables=NA)
  data$delabel[dupl_loc] = data$full_name[dupl_loc]
  
  # Create a plot and then plot a one-to-one ratio line. 
  base_plot = ggplot2::ggplot(data, ggplot2::aes_string(x = comp1, y = comp2, col = "reg")) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_abline(slope = 1, color = line_color, linewidth = 1.5) +
    ggplot2::geom_abline(slope = 1, intercept = (-1*quants[2]), color = int_color, linetype = "dashed") +
    ggplot2::geom_abline(slope = 1, intercept = (-1*quants[1]), color = int_color, linetype = "dashed") +
    ggplot2::scale_color_manual(name = "Legend", values = mycolors) + 
    ggplot2::theme_bw() +
    ggplot2::ylab(comp2) +
    ggplot2::xlab(comp1) + 
    ggplot2::theme(legend.position="none", axis.text = ggplot2::element_text(size = axes_text_size),
          axis.title = ggplot2::element_text(size = axes_label_size))
  
  # find the maximum values for the X and Y axes. 
  y_max = abs_max(data[[comp2]])
  x_max = abs_max(data[[comp1]])
  
  # Add the labels. 
  label_plot = base_plot + ggrepel::geom_label_repel(ggplot2::aes(label = delabel),
                                              size = text_size,
                                              min.segment.length = 0.000001,
                                              label.size = NA,
                                              max.overlaps = 10000,
                                              box.padding   = box_pad, 
                                              point.padding = point_pad,
                                              na.rm = T,
                                              segment.color = 'grey50',
                                              segment.size = 0.3,
                                              arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches")),
                                              show.legend = FALSE) + 
      ggplot2::coord_cartesian(xlim = c(-x_max-1, x_max+1),
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
##    title = A string that will be used as the title of the plot. 
##    order = A list of strings that correspond to the named elements in the 
##            dataframes list. The order of this list will become the order of 
##            the columns in the heatmap. 
##    label_options = The parts of the label that will be included. 
#######################################################################

hmap_prep = function (dataframes, title = "", order = c(),
                      label_options = c("Name", "P_site", "Identifier")) {
  
  # add the name of the comparison to a column on the respective dataframe
  for (i in 1:length(dataframes)) {
    dataframes[[i]]$comp = names(dataframes[i])
  }
  
  # Take all of the dataframes and bind them together.
  bound = dplyr::bind_rows(dataframes)
  bound = data.frame(bound)
  
  # Define a new column which will be the names of our matrix. 
  bound$Name = do.call(base::paste, c(bound[label_options], sep=" "))
  
  # subset the bound data in order to focus on the meaningful columns. 
  bound = base::subset(bound, select = c("Name", "comp", "log2_FC", "full_name"))
  
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
##    text_size = The size of the text used for both the X and Y axes. 
##    lg_title_size = The size of the text used for the legend title. 
##    lg_text_size = The size of the text used for the tick marks on the legend. 
##    color choice = The RColorBrewer palette to be used for the color bar. 
##    reverse_scale = A boolean which specifies if the color bar should be reversed. 
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
  heatmap = plotly::plot_ly(colors = color_choice) %>%
    plotly::add_heatmap(x = base::colnames(bound_mat), y = base::rownames(bound_mat), z = bound_mat, reversescale = reverse_scale,
                text = bound_mat2, zmin = heat_num[1], zmax = heat_num[2],
                hovertemplate = base::paste('Name: %{y}<extra></extra><br>',
                                      'Comparison: %{x}<br>',
                                      'Log2 FC: %{text}'),
                colorbar = list(limits = c(heat_num[1], heat_num[2]),
                                len = (120 + (mat_length*height_hmap)/2), lenmode = "pixels",
                                title = list(text = "log2 FC", font = list(size = lg_title_size)),
                                tickfont = list(size = lg_text_size), yanchor = "middle")) %>%
    plotly::layout(
      xaxis = list(tickfont = list(size = text_size)),
      yaxis = list(tickfont = list(size = text_size))) %>%
    plotly::config(modeBarButtons = list(list("toImage")),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE, displayModeBar = TRUE)
  
  return(heatmap)
}

########################################################################################
## clean kinexus
##    path = The path to a raw data file provided by Kinexus.
##    max_error = The maximum percent error that is tolerated between technical replicates
##                where percent error is calculated as:
##                avg(replicate intensity - mean intensity)/avg(intensity) * 100
#########################################################################################

clean_kinexus = function(infile, max_error = 50) {
  
  raw_datatables = suppressMessages(purrr::map(infile$datapath, readr::read_tsv))
  
  datatables = list()
  
  for (i in 1:length(raw_datatables)) {
    # Make a certain row the new column name and then subsequently delete that row.  
    colnames(raw_datatables[[i]]) = as.character(raw_datatables[[i]][53,])
    raw_datatables[[i]] = raw_datatables[[i]][-53,]
    colnames(raw_datatables[[i]])
    
    # select the columns of interest. 
    raw_datatables[[i]] = subset(raw_datatables[[i]], select = c("Target Name with alias", "Human P-Site", "Cat. No.",
                                                                 "Signal Median", "Background Median", "Signal Area", "Flag"))
    
    # rename the columns
    colnames(raw_datatables[[i]]) = c("Target_Name", "P_Site", "Antibody",
                                      "Signal_Median", "Background_Median", "Spot_Area", "Flag")
    
    # remove the na containing rows. 
    raw_datatables[[i]] = tidyr::drop_na(raw_datatables[[i]])
    
    # convert to numeric 
    raw_datatables[[i]] = suppressWarnings(dplyr::mutate(raw_datatables[[i]], Signal_Median = as.numeric(Signal_Median),
                                                         Background_Median = as.numeric(Background_Median),
                                                         Spot_Area = as.numeric(Spot_Area),
                                                         Flag = as.numeric(Flag)))
    
    # calulate the raw intensity. 
    raw_datatables[[i]] = dplyr::mutate(raw_datatables[[i]], raw_intensity = (Signal_Median - Background_Median)*(Spot_Area/100))
    
    # Calculate the scalar needed for normalization
    summed_raw_intensity = sum(raw_datatables[[i]]$raw_intensity)
    scalar = 20000000/summed_raw_intensity
    
    # calculate the normalized intensity. 
    raw_datatables[[i]] = dplyr::mutate(raw_datatables[[i]], Normalized_intensity = raw_intensity*scalar)
    
    # Remove antibodies with flags equal to 1. 
    flagged = dplyr::filter(raw_datatables[[i]], Flag == 1 | Flag == 2)$"Antibody"
    raw_datatables[[i]] = dplyr::filter(raw_datatables[[i]], !raw_datatables[[i]]$"Antibody" %in% flagged)
    
    # use group by to calculate the average Intensity. 
    dt_group = dplyr::group_by(raw_datatables[[i]], Antibody)
    dt = suppressMessages(dplyr::summarise(dt_group, Target_Name = Target_Name, P_Site = P_Site,
                                           Average_intensity = mean(Normalized_intensity),
                                           Normalized_intensity = Normalized_intensity))
    
    # Calculate the error ranges
    dt$diff = abs(dt$Average_intensity - dt$Normalized_intensity)
    dt = suppressMessages(dplyr::summarise(dt, Target_Name = Target_Name, P_Site = P_Site,
                                           Average_intensity = Average_intensity,
                                           Error_range = mean(diff)))
    
    # remove duplicate rows 
    dt = dplyr::distinct(dt, Antibody, .keep_all = TRUE)
    
    # Calculate the percent error. 
    dt$percent_error = (dt$Error_range/dt$Average_intensity)*100
    
    # Remove antibodies with % error above a certain cutoff. 
    dt = dplyr::filter(dt, percent_error < max_error)
    
    # Remove everything in brackets from the protein names.
    dt$Target_Name = gsub(r"{\s*\([^\)]+\)}","",as.character(dt$Target_Name))
    
    # Sort by target name and rename columns. 
    dt = dt[order(dt$Target_Name),]
    dt = dt[, c("Target_Name", "P_Site", "Antibody", "Average_intensity", "Error_range", "percent_error")]
    colnames(dt) = c("Name", "P_Site", "Identifier", "Rep_1", "Error_range", "percent_error")
    
    # Select only the columns of interest
    dt = subset(dt, select = c("Name", "P_Site", "Identifier", "Rep_1"))
    datatables = append(datatables, list(data.table::as.data.table(dt)))
  }
  
  return(datatables)
}

#######################################################################
## abs_max. A helper function which allows for significantly more flexible rounding. 
#######################################################################
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

