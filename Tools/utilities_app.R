##########################compare_CyberT ################################################
##    infile = a shiny fileInput dataframe. 
##    controls = a vector of strings that correspond to the named elements of 
##              dataframes. The strings in this list will be compared to the 
##              strings in "treatment" on a one-to-one basis and will be treated 
##              as the control in all calculations (e.g. fold change is control/treatment)
##    treatments = a vector of strings that correspond to the named elements of 
##                dataframes. The strings in this list will be compared to the 
##                strings in "control" on a one-to-one basis and will be treated 
##                as the treatment in all calculations (e.g. fold change is control/treatment)
##    analysis = a string which specifies if the analysis should look at only phospho site specific
##               data or only Pan specific data. Options are "Phospho-specific" and "Pan-specific"
##    apo_name = The name used in the P_site column to specify that a datapoint is pan-specific. 
##    no_conf = A boolean which specifies whether or not to set the confidence in CyberT to
##              zero. If made True, then a normal t-test is used instead of CyberT t-test.
##    var_norm = The normalization method to be used. options are "vsn", "logT", and "none"
##    kd = A boolean which specifies whether the data uploaded are kinexus data
##    fmd = A boolean which specifies whether the data uploaded are full moon biosystems data
##    max_error = The maximum acceptable error when filtering kinexus data. 
#########################################################################################

compare_CyberT = function(infile, controls, treatments, analysis = NULL,
                          apo_name = NULL, no_conf = FALSE, var_norm = "vsn",
                          kd = FALSE, fmd = FALSE, max_error = 50, phospho = F,
                          pho_spec = F){

  # Look at the datatype being used and decide how to upload the data.
  if (kd) {
    datatables = clean_kinexus(infile, max_error = max_error)
    names(datatables) = stringr::str_remove(infile$name, ".txt")
  } else if (fmd) {
    datatables = clean_full_moon(infile, phospho = phospho, pho_spec = pho_spec)
  } else {
    datatables = purrr::map(infile$datapath, data.table::fread)
    names(datatables) = stringr::str_remove(infile$name, ".csv")
  }

  # Create a copy of all of the column options used.
  col_options = na.omit(stringr::str_extract(stringr::str_to_sentence(colnames(datatables[[1]])),
                                             pattern = "Name|P_site|Identifier"))
  col_options = factor(col_options, levels = c("Name", "P_site", "Identifier"))
  col_options = col_options[order(col_options)]
  col_options = as.character(col_options)

  print(class(datatables[[1]]))
  #Sys.sleep(30)

  # Clean up the contents of the datatables.
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

  # Create empty lists to store the data comparisons.
  comparisons = list()
  comp_names = c()

  for (i in 1:length(controls)) {
    print(base::paste("Processing Comparison ", toString(i), "/", toString(length(controls)), sep = ""))

    # Find the names to be used by cybert for controls and treats.
    numC = sum(!colnames(datatables[[controls[i]]]) == "full_name")
    numE = sum(!colnames(datatables[[treatments[i]]]) == "full_name")

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
    dt[, c(1:(numC+numE)) := NULL]
    dt[, Fold_change := Treated_Average/Control_Average]
    dt[, log2_FC := log2(Fold_change)]
    dt[, P_value_adjust := results$runMulttest.pVal.]
    dt[, log_p := -log10(P_value_adjust)]

    # append to lists
    comparisons = base::append(comparisons, list(dt))
    comp_names = base::append(comp_names, base::paste(controls[i], "vs", treatments[i], sep = "_"))

    # Determine if old datatables can be deleted.
    if (!names(datatables[controls[i]]) %in% controls[i+1:length(controls)] &&
        !names(datatables[controls[i]]) %in% treatments[i+1:length(treatments)]) {
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

########################## volcano_plot_app #############################################
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
########################################################################################

volcano_plot_app = function (data, to_label = c(), top = 0, FC_range = c(-1,1), P_cutoff = 3,
                             mycolors = c("blue", "red", "black", "green"), text_size = 3,
                             axes_text_size = 2, axes_label_size = 3,  point_sizes = c(1, 2),
                             range = NULL, box_pad = 0.3, point_pad = 0.2, 
                             label_options = c("Name", "P_site", "Identifier"),
                             sig_label = FALSE) {
  
  # Make sure the data uploaded is a dataframe. 
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
  
  # Create a column with the euclidean distances to the originn and choose a set number to label.
  if (top > 0) {
    data$dist = sqrt(data$log_p^2 + data$log2_FC^2)
    data$dist[data$reg == "Non-significant"] = NA
    data = data[order(data$dist, decreasing = T),]
    top_reg = data[is.na(data$dist) == F,][1:top,]$full_name
    data$delabel[data$full_name %in% top_reg] = data$label_form[data$full_name %in% top_reg] 
  } else {top_reg = c()}
  
  # Check to see if the maximum labels have been overcome. 
  if ((length(unique(c(to_label, top_reg))) > 150)) {stop(shiny::safeError("Exceeded maximum number of labels (150)"))}
  
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

######################### scatter_plot_app #############################################
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
######################################################################################

scatter_plot_app = function (data, comp1, comp2, to_label = c(), top = 0, quant_int = c(0.01, 0.99),
                           line_color = "blue", int_color = "red", text_size = 2.3,  point_size = 1,
                           axes_text_size = 2, axes_label_size = 3, box_pad = 0.3, point_pad = 0.2,
                           mycolors = c("black", "black", "black", "green"), sig_label = FALSE,
                           label_options = c("Name", "P_site", "Identifier")) {
  
  # Make sure the data uploaded is a dataframe.
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
  if ((length(unique(c(to_label, top_reg))) > 150)) {stop(shiny::safeError("Exceeded maximum number of labels (150)"))}
  
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

######################### abs_max #######################################################
## A helper function which finds the maximum value in a vector regardless of sign. 
#########################################################################################

abs_max = function(vector) {
  # Define the upper and lower limits based on the maximum and minimum. 
  up_lim = max(vector, na.rm = T)
  down_lim = min(vector, na.rm = T)
  max_lim = max(c(abs(up_lim), abs(down_lim)))
  return(max_lim)
}

######################### hmap_prep #####################################################
##  hmap_prep
##    dataframes = A list of dataframes outputted by the compare_CyberT function
##    title = A string that will be used as the title of the plot. 
##    order = A list of strings that correspond to the named elements in the 
##            dataframes list. The order of this list will become the order of 
##            the columns in the heatmap. 
##    label_options = The parts of the label that will be included. 
#########################################################################################

  hmap_prep = function (dataframes, title = "", order = c(),
                        label_options = c("Name", "P_site", "Identifier")) {
    
    # add the name of the comparison to a column on the respective dataframe
    for (i in 1:length(dataframes)) {
      dataframes[[i]]$comp = names(dataframes[i])
    }
    
    # Take all of the dataframes and bind them together.
    bound = data.table::rbindlist(dataframes)
    bound = as.data.frame(bound)
    
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

######################### hmap ############################################################
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
#########################################################################################

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

######################### clean_kinexus #################################################
##    path = The path to the raw data files provided by Kinexus.
##    max_error = The maximum percent error that is tolerated between technical replicates
##                where percent error is calculated as:
##                avg(replicate intensity - mean intensity)/avg(intensity) * 100
#########################################################################################

clean_kinexus = function(infile, max_error = 50) {
  
  # Load up the raw Kinexus files. 
  raw_datatables = purrr::map(infile$datapath, data.table::fread)
  
  # Create a list of for_loops to be populated during the analysis. 
  datatables = list()
  for (i in 1:length(raw_datatables)) {
    print(base::paste("Cleaning Kinexus File ", toString(i), "/", toString(length(raw_datatables)), sep = ""))
    
    # Make a certain row the new column name and then subsequently delete that row.  
    colnames(raw_datatables[[i]]) = as.character(raw_datatables[[i]][57,])
    raw_datatables[[i]] = raw_datatables[[i]][-57,]
    
    # select the columns of interest. 
    raw_datatables[[i]] = raw_datatables[[i]][, c("Target Name with alias", "Human P-Site", "Cat. No.",
                                                  "Signal Median", "Background Median", "Signal Area", "Flag")]
    
    # rename the columns
    colnames(raw_datatables[[i]]) = c("Target_Name", "P_Site", "Antibody",
                                      "Signal_Median", "Background_Median", "Spot_Area", "Flag")
    
    # remove the na containing rows. 
    raw_datatables[[i]] = na.omit(raw_datatables[[i]])
    
    # convert to numeric 
    raw_datatables[[i]][, Signal_Median := as.numeric(Signal_Median)]
    raw_datatables[[i]][, Background_Median := as.numeric(Background_Median)]
    raw_datatables[[i]][, Spot_Area := as.numeric(Spot_Area)]
    raw_datatables[[i]][, Flag := as.numeric(Flag)]
    
    # calulate the raw intensity. 
    raw_datatables[[i]][, raw_intensity := (Signal_Median - Background_Median)*(Spot_Area/100)]
    
    # Calculate the scalar needed for normalization
    summed_raw_intensity = sum(raw_datatables[[i]]$raw_intensity, na.rm = T)
    scalar = 20000000/summed_raw_intensity
    
    # calculate the normalized intensity. 
    raw_datatables[[i]][, Normalized_intensity := raw_intensity*scalar]
    
    # Remove antibodies with flags equal to 1. 
    flagged = raw_datatables[[i]][Flag == 1 | Flag == 2,]$"Antibody"
    raw_datatables[[i]] = raw_datatables[[i]][!raw_datatables[[i]]$"Antibody" %in% flagged,]
    
    # Calculate the average Intensity. 
    raw_datatables[[i]][, Average_intensity:=mean(Normalized_intensity), by=list(Antibody)]
    
    # Calculate the error ranges
    raw_datatables[[i]][, Error_range:=mean(abs(Average_intensity - Normalized_intensity)), by = list(Antibody)]
    
    # remove duplicate rows 
    raw_datatables[[i]] = unique(raw_datatables[[i]], by = "Antibody")
    
    # Calculate the percent error. 
    raw_datatables[[i]][, percent_error := (Error_range/Average_intensity)*100]
    
    # Remove antibodies with % error above a certain cutoff. 
    raw_datatables[[i]] = raw_datatables[[i]][percent_error < max_error,]
    
    # Remove everything in brackets from the protein names.
    raw_datatables[[i]][, Target_Name := gsub(r"{\s*\([^\)]+\)}","",as.character(Target_Name))]
    
    # Sort by target name and rename columns. 
    raw_datatables[[i]] = raw_datatables[[i]][order(Target_Name),]
    raw_datatables[[i]] = raw_datatables[[i]][, c("Target_Name", "P_Site", "Antibody", "Average_intensity", "Error_range", "percent_error")]
    colnames(raw_datatables[[i]]) = c("Name", "P_Site", "Identifier", "Rep_1", "Error_range", "percent_error")
    
    # Select only the columns of interest
    raw_datatables[[i]] = raw_datatables[[i]][,c("Name", "P_Site", "Identifier", "Rep_1")]
    datatables = append(datatables, list(raw_datatables[[i]]))
    raw_datatables[[i]] = data.frame(c())
  }
  
  return(datatables)
}

######################### clean_full_moon ###############################################
##    path = The path to a raw data file provided by Kinexus.
##    phospho = A boolean which specifies if working with phosphorylation data or
##              expression data. 
##    pho_spec = A boolean which specifies if looking at site specific Ab data or 
##               phosphorylation specific Phospho- data. 
#########################################################################################

clean_full_moon = function(infile, phospho = FALSE, pho_spec = FALSE) {
  # Load in the excel file. 
  fm_data = readxl::read_excel(infile$datapath, sheet = "Assay Data")
  
  # Extract just the part around the Average Signal Data.  
  as_col = which(stringr::str_detect(fm_data, "Average Signal"))
  fm_data = fm_data[(as_col - 1):length(fm_data)]
  
  # Remove the empty rows.
  as_i = which(stringr::str_detect(c(fm_data[[2]]), "Average Signal"))
  fm_data = fm_data[-(1:(as_i-1)),]
  
  # Remove the wealth of unnecessary in-between data. 
  as_col = which(stringr::str_detect(fm_data, "Average Signal"))
  dn_col = which(stringr::str_detect(fm_data[1,], "Data Normalized"))
  fm_data = fm_data[-(as_col:(dn_col-1))]
  
  # Remove the needless bottom 4 rows. 
  fm_data = fm_data[-((nrow(fm_data)-3):nrow(fm_data)),]
  
  # Extract the proper names from the data. 
  if (phospho) {
    # Save the antibody names
    Antibodies = fm_data[1][-(1:2),]
    names(Antibodies) = "Name"
    
    # Extract everything within and outside of parenthesis. 
    phospho = apply(Antibodies, 1, extract_parenthesis)
    Name = apply(Antibodies, 1, exclude_parenthesis)
    
    # Extract the site itself. 
    if (pho_spec) {choice = "Phospho-"} else {choice = "Ab-"}
    print(choice)
    site = c()
    for (i in 1:length(phospho)) {
      if (stringr::str_detect(phospho[i], choice)) {
        site[i] = stringr::str_split(phospho[i], choice, simplify = TRUE)[,2]
      } else {site[i] = NA}
    }
    
    # Create the new antibodies dataframe. 
    Antibodies = data.frame(Name = Name, P_site = site)
    fm_data = fm_data[-1]
  } else {
    Antibodies = fm_data[1][-(1:2),]
    names(Antibodies) = "Name"
    fm_data = fm_data[-1]
  }
  
  # Convert dataframe into list of dataframes. 
  if ("Group Mean" %in% fm_data[1,]) {
    
    # Recalculate the indexes of important sections. 
    dn_col = which(stringr::str_detect(fm_data[1,], "(?i)Data Normalized"))
    gm_col = which(fm_data[1,] == "Group Mean")
    gc_col = which(fm_data[1,] == "Group CV")
    
    # Define the filenames based on the group names. 
    file_names = as.character(fm_data[gm_col:(gc_col - 2)][2,], na.rm = TRUE)
    
    # Remove the group data. 
    fm_data = fm_data[-((gm_col-1): length(fm_data))]
    
    # Find the number of replicates per group. 
    N_reps = length(fm_data[dn_col:(gm_col - 2)])/length(file_names)
    
    # Remove the unneeded rows. 
    fm_data = fm_data[-(1:2),]
    
    files = list()
    for (i in 1:length(file_names)) {
      temp = fm_data[i:(i + N_reps - 1)]
      rep_names = sprintf("Rep_%s",seq(1:N_reps))
      colnames(temp) = c(rep_names)
      temp = sapply(temp,as.numeric)
      temp = cbind(Antibodies, temp)
      temp = data.table::as.data.table(temp)
      files = append(files, list(temp))
      temp = NULL
    }
    names(files) = file_names
    Antibodies = NULL
    fm_data = NULL
    
  } else {
    # Define the filenames based on the sample names. 
    file_names = as.character(fm_data[2,])
    
    # Remove extra columns if necessary. 
    if (any(file_names == "NA")) {
      x = which(file_names == "NA")
      fm_data = fm_data[-(x:length(fm_data))]
      file_names = file_names[-(x:length(file_names))]
    }
    
    # Find the number of samples.  
    N_samps = length(file_names)
    
    # Remove the unneeded rows. 
    fm_data = fm_data[-(1:2),]
    
    # Save all samples
    files = list()
    for (i in 1:length(file_names)) {
      temp = fm_data[i]
      colnames(temp) = c("Rep_1")
      temp = sapply(temp,as.numeric)
      temp = cbind(Antibodies, temp)
      temp = tidyr::drop_na(temp)
      temp = data.table::as.data.table(temp)
      files = append(files, list(temp))
      temp = NULL
    }
    names(files) = file_names
    Antibodies = NULL
    fm_data = NULL
  }
  return(files)
}

get_full_moon_names = function(infile) {
  # Load in the excel file. 
  fm_data = readxl::read_excel(infile$datapath, sheet = "Assay Data")
  
  # Extract just the part around the Average Signal Data.  
  as_col = which(stringr::str_detect(fm_data, "Average Signal"))
  fm_data = fm_data[(as_col - 1):length(fm_data)]
  
  # Remove the empty rows.
  as_i = which(stringr::str_detect(c(fm_data[[2]]), "Average Signal"))
  fm_data = fm_data[-(1:(as_i-1)),]
  
  # Remove the wealth of unnecessary in-between data. 
  as_col = which(stringr::str_detect(fm_data, "Average Signal"))
  dn_col = which(stringr::str_detect(fm_data[1,], "Data Normalized"))
  fm_data = fm_data[-(as_col:(dn_col-1))]
  fm_data = fm_data[-1]
  
  # Convert dataframe into list of dataframes. 
  if ("Group Mean" %in% fm_data[1,]) {
    
    # Recalculate the indexes of important sections. 
    gm_col = which(fm_data[1,] == "Group Mean")
    gc_col = which(fm_data[1,] == "Group CV")
    
    # Define the filenames based on the group names. 
    file_names = as.character(fm_data[gm_col:(gc_col - 2)][2,], na.rm = TRUE)
    
  } else {
    # Define the filenames based on the sample names. 
    file_names = as.character(fm_data[2,])
    
    # Remove extra columns if necessary. 
    if (any(file_names == "NA")) {
      x = which(file_names == "NA")
      fm_data = fm_data[-(x:length(fm_data))]
      file_names = file_names[-(x:length(file_names))]
    }
  }
  return(file_names)
}

######################### round_any #####################################################
## A helper function that lets one round to any decimal position.  
#########################################################################################
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

######################### extract/exclude parentheses ###################################
## Helper functions that allow for part of a character string within/outside of parentheses
## to be isolated. 
#########################################################################################

extract_parenthesis = function(x) {
  y = gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", as.character(x), perl=T)
  return(y)
}

exclude_parenthesis = function(x) {
  y = gsub(r"{\s*\([^\)]+\)}","",as.character(x), perl=T)
  return(y)
}

