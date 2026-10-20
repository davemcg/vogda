mega_maker <- function(jid_of_interest,
                       up_delta = 1000,
                       down_delta = 1000,
                       logFC_gt = -10){

  # Gene Expression Plots

  ## Across the body
  g <- diff$sj_topsplice_exon |> filter(jid == jid_of_interest) |> pull(SYMBOL)

  res <- geyser::geyser_boxplot(
    rse = cts,
    features = g,
    groupings = c("Tissue"),
    slot = "cpm",
    feature_col = "gene_name",
    color_by = "Tissue",
    color_palette = "okabe",
    expression_scale = TRUE,
    show_points = TRUE,   point_size = 0.5,

    data_source_name = "EiaD"
  )

  body_plot <- res$plot +coord_flip() +
    theme(plot.margin = unit(c(1,1,1,20), "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(caption = '')



  ## Across adult and fetal retina

  retina_cts <- cts[,(cts$Tissue == 'Retina' & !grepl("AMD", cts$Perturbation))]
  retina_cts$Age[is.na(retina_cts$Age)] <- "Fetal"


  res <- geyser::geyser_boxplot(
    rse = retina_cts,
    features = g,
    groupings = c("Age"),
    slot = "cpm",
    feature_col = "gene_name",
    color_by = "Source",
    color_palette = "okabe",
    expression_scale = TRUE,
    show_points = TRUE,   point_size = 0.5,
    data_source_name = "EiaD"
  )

  adult_fetal_retina_plot <- res$plot +coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.1) # Coordinates from 0 to 1
    ) +
    labs(caption = '')

  #
  #
  # <!-- ## Exon Counts (Retina only) -->
  #   <!--
  #   <!-- retina_exon <- exon[,(exon$Tissue == 'Retina' & !grepl("AMD", exon$Perturbation))] -->
  #   <!-- res <- geyser::geyser_boxplot( -->
  #                                         <!--   rse = retina_exon, -->
  #                                         <!--   features = gene, -->
  #                                         <!--   groupings = c("Tissue", "Age", "Source"), -->
  #                                         <!--   slot = "cpm", -->
  #                                         <!--   feature_col = "gene_name", -->
  #                                         <!--   color_by = "Tissue", -->
  #                                         <!--   color_palette = "okabe", -->
  #                                         <!--   expression_scale = TRUE, -->
  #                                         <!--   show_points = FALSE, -->
  #                                         <!--   data_source_name = "EiaD", -->
  #                                         <!--   facet_ncol = 5 -->
  #                                         <!-- ) -->
  #
  #   <!-- res$plot + coord_flip() +  -->
  #   <!--   theme(axis.text.x = element_text(angle = 90, hjust = 1), -->
  #                  <!--         legend.position = "none") -->
  #   <!--
  #
  #   <!-- ## Junction Counts (Retina only) -->
  #   <!--
  #   <!-- retina_sj <- sj[,(sj$Tissue == 'Retina' & !grepl("AMD", sj$Perturbation))] -->
  #   <!-- res <- geyser::geyser_boxplot( -->
  #                                         <!--   rse = retina_sj, -->
  #                                         <!--   features =gene, -->
  #                                         <!--   groupings = c("Tissue", "Age", "Source"), -->
  #                                         <!--   slot = "cpm", -->
  #                                         <!--   feature_col = "SYMBOL", -->
  #                                         <!--   color_by = "Tissue", -->
  #                                         <!--   color_palette = "okabe", -->
  #                                         <!--   expression_scale = TRUE, -->
  #                                         <!--   show_points = FALSE, -->
  #                                         <!--   data_source_name = "EiaD", -->
  #                                         <!--   facet_ncol = 5 -->
  #                                         <!-- ) -->
  #
  #   <!-- res$plot +coord_flip() +  -->
  #   <!--   theme(axis.text.x = element_text(angle = 90, hjust = 1), -->
  #                  <!--         legend.position = "none") -->
  #   <!--
  #

  # Plot Upstream Exon, Junction, and Downstream Exons Counts


  retina_exon <- exon[,(exon$Tissue == 'Retina' & !grepl("AMD", exon$Perturbation))]
  retina_exon$Age[is.na(retina_exon$Age)] <- "Fetal"
  retina_sj <- sj[,(sj$Tissue == 'Retina' & !grepl("AMD", sj$Perturbation))]
  retina_sj$Age[is.na(retina_sj$Age)] <- "Fetal"
  exon_rd <- as.data.table(as.data.frame(rowData(exon)))
  exon_rd[,'id'] <- row.names(rowData(exon))

  sj_rd <- as.data.table(as.data.frame(rowData(sj)))
  junc_gr <- GRanges(
    seqnames = sj_rd$chr,
    ranges   = IRanges(start = sj_rd$start, end = sj_rd$end),
    strand   = sj_rd$strand,
    jid      = sj_rd$jid,
    SYMBOL   = sj_rd$SYMBOL
  )

  exon_gr <- GRanges(
    seqnames = exon_rd$Chr,
    ranges   = IRanges(start = exon_rd$start, end = exon_rd$end),
    strand   = exon_rd$strand,
    exon_id  = exon_rd$id,
    SYMBOL   = exon_rd$gene_name
  )

  upstream_hits   <- follow(junc_gr, exon_gr)
  downstream_hits <- precede(junc_gr, exon_gr)


  j2e_map <- data.table::data.table(
    jid             = mcols(junc_gr)$jid,
    SYMBOL          = mcols(junc_gr)$SYMBOL,
    upstream_exon   = mcols(exon_gr)$exon_id[upstream_hits],
    downstream_exon = mcols(exon_gr)$exon_id[downstream_hits]
  )


  j2e_map |>
    filter(jid == diff$sj_topsplice_exon |> filter(SYMBOL == g) |>
             row.names())








  upstream <- j2e_map |> filter(jid == jid_of_interest) |> pull(upstream_exon) |> gsub("_",":", x= _)
  downstream <- j2e_map |> filter(jid ==jid_of_interest) |> pull(downstream_exon) |> gsub("_",":", x= _)
  res <- geyser::geyser_boxplot(
    rse = retina_exon[upstream,],
    features = upstream,
    groupings = c("Age"),
    slot = "cpm",
    feature_col = "row names",
    color_by = "Source",
    color_palette = "okabe",
    expression_scale = TRUE,
    show_points = TRUE,   point_size = 0.5,
    data_source_name = "EiaD"
  )

  up_plot <- res$plot +coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") + ggtitle("Upstream Exon") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1, r = 1, b = -35, l = 1, unit = "pt")) +
    labs(y = '', caption = '')

  res <- geyser::geyser_boxplot(
    rse = retina_sj[jid_of_interest,],
    features = jid_of_interest,
    groupings = c("Age"),
    slot = "cpm",
    feature_col = "row names",
    color_by = "Source",
    color_palette = "okabe",
    expression_scale = TRUE,
    show_points = TRUE,   point_size = 0.5,
    data_source_name = "EiaD"
  )

  sj_plot <- res$plot +coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") + ggtitle("Junction") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1, r = 1, b = -35, l = 1, unit = "pt")) +
    labs(y = '', caption = '')



  res <- geyser::geyser_boxplot(
    rse = retina_exon[downstream,],
    features = downstream,
    groupings = c("Age"),
    slot = "cpm",
    feature_col = "row names",
    color_by = "Source",
    color_palette = "okabe",
    expression_scale = TRUE,
    show_points = TRUE,   point_size = 0.5,
    data_source_name = "EiaD"
  )

  down_plot <- res$plot +coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") + ggtitle("Downstream Exon") +
    labs(caption = '')


  adult_fetal_junc <- up_plot / plot_spacer() / sj_plot / plot_spacer() / down_plot + plot_layout(heights = c(1, -0.3, 1, -0.4, 1))





  # supeRwig Coverage


  ctx <- supeRwig:::build_context(
    se_dir = '~/data/vogda/eiad/EiaD_2026_hdf5/',
    sj_se_dir  = '~/data/vogda/eiad/se_h5/se_sj_human.20261007.hdf5/',
    bigwig_base = '~/data/vogda/eiad/bigwig/',
    anno_fst = '~/git/EiaD_build/data/gencode.v48.primary_assembly.annotation.table.fst',
    bigwig_ext = '.bw',
    n_workers = 10
  )
  # only keep retina
  ctx$eiad_meta <- ctx$eiad_meta |> dplyr::filter(Tissue == 'Retina',
                                                  Sub_Tissue != 'Endothelial' | is.na(Sub_Tissue),
                                                  study_accession != 'SRP287152')
  # SRP287152 has some poor quality samples
  # only keep ref junctions
  ctx$sj_row_meta <- ctx$sj_row_meta[SummarizedExperiment::rowData(ctx$sj_se)$annot == 1,]
  # rename fetal data
  ctx$eiad_meta$Age[is.na(ctx$eiad_meta$Age)] <- ctx$eiad_meta$Source[is.na(ctx$eiad_meta$Age)]



  chr <- j2e_map |> filter(jid == diff$sj_topsplice_gene |> filter(SYMBOL == g) |> row.names()) |> pull(upstream_exon) |> str_split(":|-|_") |> unlist() |> (\(x) x[1])()
  up_coord <- upstream |> str_split(":|-|_") |> unlist() |> (\(x) x[2])() |> as.integer()
  down_coord <- downstream |> str_split(":|-|_") |> unlist() |> (\(x) x[3])() |> as.integer()

  up_coord2 <- min(up_coord, down_coord)
  down_coord2 <- max(up_coord, down_coord)
  #2. Script plot generation directly to self-contained HTML files
  # a <- supeRwig::plot_wiggle(
  #   ctx,
  #   chr           = chr,
  #   start         = up_coord2-1000,
  #   end           = down_coord2+1000,
  #   facet_cols      = c("Tissue", "Age"),
  #   color_by        = "Tissue",
  #   show_junctions  = TRUE,
  #   #target_gene   = "ANK2",
  #   min_psi5        = 10,
  #   min_psi3        = 10,
  #   overlap_factor = 0.8,
  #   mean_average = FALSE,
  #   interactive = TRUE
  # )

  b <- supeRwig:::plot_region_heatmap(
    ctx           = ctx,
    max_samples = 5,
    #   min_expr = 100,
    chr           = chr,

    start         = up_coord2 - up_delta,
    end           = down_coord2+ down_delta,
    facet_cols    = c("Tissue", "Age"),
    #target_gene   = "DPY19L2P2",
    palette       = "magma",
    show_transcripts = TRUE,        # show tx models over heatmap
    show_column_axis = FALSE
  )
  #b
  #a$plot



  # Super Transcript



  # 1. Extract junction metadata and target exon IDs
  junc_row <- sj_rd |>
    filter(jid == jid_of_interest) |>
    head(1)

  junc_diff_info <- diff$sj_topsplice_exon |>
    filter(jid == jid_of_interest) |>
    head(1)

  fc_label <- if (nrow(junc_diff_info) > 0 && "logFC" %in% colnames(junc_diff_info)) {
    paste0(" (log2FC: ", round(junc_diff_info$logFC, 2), ")")
  } else {
    ""
  }

  upstream_raw <- j2e_map |> filter(jid == jid_of_interest) |> pull(upstream_exon) |> head(1)
  downstream_raw <- j2e_map |> filter(jid == jid_of_interest) |> pull(downstream_exon) |> head(1)

  # 2. Build exon table
  tx_exons <- exon_rd |>
    filter(gene_name == g) |>
    as.data.table()

  tx_exons[, norm_id := gsub("_", ":", id)]
  tx_exons[, exon_status := case_when(
    norm_id == gsub("_", ":", upstream) | id %in% c(upstream, upstream_raw) ~ "Upstream Exon",
    norm_id == gsub("_", ":", downstream) | id %in% c(downstream, downstream_raw) ~ "Downstream Exon",
    TRUE ~ "Other Exons"
  )]
  tx_exons[, tx_idx := 1]

  # De-duplicate overlapping exon definitions
  tx_exons <- unique(tx_exons, by = c("start", "end"))
  setorder(tx_exons, start)

  # 3. Compute midpoints for target upstream/downstream exons to place dots
  marked_exons <- tx_exons[exon_status %in% c("Upstream Exon", "Downstream Exon")]
  marked_exons[, mid_x := (start + end) / 2]

  # 4. Build transcript spine base (`tx_base`)
  w_start <- min(tx_exons$start) - 500
  w_end <- max(tx_exons$end) + 500
  strand_val <- tx_exons$strand[1]

  tx_base <- data.table::data.table(
    tx_idx     = 1,
    start      = min(tx_exons$start),
    end        = max(tx_exons$end),
    strand     = strand_val,
    tx_label   = '',
    draw_start = if (strand_val == "-") max(tx_exons$end) + 0.5 else min(tx_exons$start) - 0.5,
    draw_end   = if (strand_val == "-") min(tx_exons$start) - 0.5 else max(tx_exons$end) + 0.5,
    label_x    = min(tx_exons$start)
  )

  # 5. Junction layer data
  junction_plot_df <- data.table::data.table(
    start           = as.numeric(junc_row$start),
    end             = as.numeric(junc_row$end),
    strand          = as.character(junc_row$strand),
    transcript_name = g,
    tx_idx          = 1,
    label           = jid_of_interest
  )

  # 6. Build plot
  super_tx_plot <- ggplot2::ggplot() +
    # Intron / Transcript Spine Segment
    ggplot2::geom_segment(
      data = tx_base,
      ggplot2::aes(x = draw_start, xend = draw_end, y = tx_idx, yend = tx_idx),
      color = "grey40",
      linewidth = 0.6,
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.08, "inches"),
        type = "closed",
        ends = if (strand_val == "-") "first" else "last"
      )
    ) +
    # Uniform Exon Rectangles (all grey75)
    ggplot2::geom_rect(
      data = tx_exons,
      ggplot2::aes(
        xmin = start - 0.5, xmax = end + 0.5,
        ymin = tx_idx - 0.2, ymax = tx_idx + 0.2
      ),
      fill = "grey75",
      color = "black",
      linewidth = 0.3
    ) +
    # Marker Dots below Upstream and Downstream Exons
    ggplot2::geom_point(
      data = marked_exons,
      ggplot2::aes(x = mid_x, y = tx_idx - 0.35, color = exon_status),
      size = 3
    ) +
    # Splice Junction Arc
    ggtranscript::geom_junction(
      data = junction_plot_df,
      ggplot2::aes(xstart = start, xend = end, y = tx_idx),
      color = "#D55E00",
      linewidth = 0.8,
      junction.y.max = 1
    ) +
    # Junction Label Annotation
    ggplot2::geom_label(
      data = junction_plot_df,
      ggplot2::aes(x = (start + end) / 2, y = 1.35, label = label),
      color = "#D55E00",
      fill = "white",
      size = 3,
      label.padding = ggplot2::unit(0.15, "lines"),
      label.size = 0.2
    ) +
    # Transcript Name Label
    ggplot2::geom_text(
      data = tx_base,
      ggplot2::aes(x = label_x, y = tx_idx - 0.35, label = tx_label),
      hjust = 0, vjust = 1, size = 3.2, fontface = "bold", color = "grey15"
    ) +
    # Color Scale for Marker Dots
    ggplot2::scale_color_manual(
      values = c(
        "Upstream Exon"   = "#0072B2",  # Okabe-Ito Blue
        "Downstream Exon" = "#009E73"   # Okabe-Ito Bluish Green
      ),
      name = "Marked Exons"
    ) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(add = c(0.5, 0.6))) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::coord_cartesian(xlim = c(w_start, w_end), clip = "off") +
    ggplot2::labs(
      x = paste0("Genomic Coordinates (", tx_exons$Chr[1], ")"),
      y = NULL,
      title = paste0(g, " Super Transcript")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = 10)
    )

  layout <- "
EEEEEEEEEE
AAAAACCDDD
AAAAACCDDD
BBBBBCCDDD
BBBBBCCDDD
"

  plot <- suppressWarnings({
    body_plot + adult_fetal_retina_plot + adult_fetal_junc + grid::grid.grabExpr(ComplexHeatmap::draw(b)) + super_tx_plot +
      plot_layout(design = layout)
  })
  return(plot)

}
