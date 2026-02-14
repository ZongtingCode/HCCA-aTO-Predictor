# ---------- Libraries ----------
library(shiny)
library(h2o)
library(ggplot2)

# ---------- Settings ----------
THRESHOLD_F1 <- 0.175312
N_EFF <- 6041
MAX_ROWS <- 11

# ---------- Robust app directory (local + shinyapps.io) ----------
get_app_dir <- function() {
  # RStudio "Run App" sometimes provides ofile
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of) && !is.na(of) && nzchar(of)) {
    return(dirname(normalizePath(of, winslash = "/", mustWork = FALSE)))
  }
  # fallback
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

APP_DIR <- get_app_dir()

MODEL_PATH  <- normalizePath(file.path(APP_DIR, "DRF_FINAL_fullData"), winslash = "/", mustWork = FALSE)
FEAT_PATH   <- normalizePath(file.path(APP_DIR, "final_feature_names.rds"), winslash = "/", mustWork = FALSE)
LEVELS_PATH <- normalizePath(file.path(APP_DIR, "final_factor_levels.rds"), winslash = "/", mustWork = FALSE)

if (!file.exists(MODEL_PATH))  stop("Model file not found at: ", MODEL_PATH)
if (!file.exists(FEAT_PATH))   stop("Feature file not found at: ", FEAT_PATH)
if (!file.exists(LEVELS_PATH)) stop("Factor level file not found at: ", LEVELS_PATH)

# ---------- H2O init (safe for shinyapps.io) ----------
h2o.init(nthreads = 1, max_mem_size = "1G")

# ---------- Load artifacts ONCE ----------
model <- h2o.loadModel(MODEL_PATH)
feature_names <- readRDS(FEAT_PATH)
levels_map <- readRDS(LEVELS_PATH)

# ---------- Helpers ----------
`%||%` <- function(a, b) if (!is.null(a)) a else b

mk_codebook <- function(labels, values) {
  setNames(as.character(values), labels)  # shown label -> submitted code
}

codebook <- list(
  sex = mk_codebook(c("Male","Female"), c(0,1)),
  ajcc8_M = mk_codebook(c("M0","M1"), c(0,1)),
  ajcc8_N = mk_codebook(c("N0","N1","N2"), c(0,1,2)),
  ajcc8_T = mk_codebook(c("T1","T2","T3","T4"), c(1,2,3,4)),
  chemotherapy = mk_codebook(c("No","Yes"), c(0,1)),
  summary_stage3 = mk_codebook(c("Localized","Regional","Distant"), c(1,2,3)),
  rural_urban_code = mk_codebook(c("Metro","Urban","Suburban","Town","Rural","Remote"), c(1,2,3,4,5,6)),
  marital_status = mk_codebook(c("Single","Married","Separated","Divorced","Widowed","Other","Unknown"), c(1,2,3,4,5,6,7)),
  race_origin = mk_codebook(c("White","Black","Asian","AI/AN","Pacific","Other"), c(1,2,3,4,5,6)),
  mets_liver = mk_codebook(c("No","Yes"), c(0,1)),
  mets_lung  = mk_codebook(c("No","Yes"), c(0,1)),
  mets_bone  = mk_codebook(c("No","Yes"), c(0,1)),
  mets_brain = mk_codebook(c("No","Yes"), c(0,1))
)

label_map <- list(
  age = "Age (years)",
  sex = "Sex",
  race_origin = "Race / Origin",
  marital_status = "Marital status",
  summary_stage3 = "Summary stage",
  ajcc8_T = "AJCC 8th T",
  ajcc8_N = "AJCC 8th N",
  ajcc8_M = "AJCC 8th M",
  chemotherapy = "Chemotherapy",
  income_10k = "Median household income (per $10k)",
  rural_urban_code = "Rural-Urban continuum",
  mets_liver = "Liver metastasis",
  mets_lung = "Lung metastasis",
  mets_bone = "Bone metastasis",
  mets_brain = "Brain metastasis"
)
pretty_label <- function(var) label_map[[var]] %||% var

align_levels <- function(df) {
  for (nm in names(levels_map)) {
    if (nm %in% names(df) && !is.null(levels_map[[nm]])) {
      df[[nm]] <- factor(df[[nm]], levels = levels_map[[nm]])
    }
  }
  df
}

wilson_ci <- function(p, n, z = 1.96) {
  p <- max(min(p, 1), 0)
  n <- max(as.numeric(n), 1)
  denom <- 1 + (z^2)/n
  center <- (p + (z^2)/(2*n)) / denom
  half <- (z * sqrt((p*(1-p) + (z^2)/(4*n))/n)) / denom
  c(lower = max(0, center - half), upper = min(1, center + half))
}

# a nicer discrete palette (no need to set ggplot colors elsewhere)
palette_vec <- c("#0066CC","#E41A1C","#54A552","#FF8000","#BA55D3",
                 "#00BFFF","#006400","#994C00","#F781BF","#A9A9A9","#0E0000")

# build a short readable “covariates string”
covariate_string <- function(input_row_named, show_labels = TRUE) {
  # input_row_named: named list of "pretty labels" already
  paste(paste0(names(input_row_named), ": ", unname(input_row_named)), collapse = " | ")
}

# ---------- UI ----------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .main-panel {padding-top: 8px;}
    .plot-wrap {margin-bottom: 14px;}
    .tbl-wrap  {margin-top: 10px;}
    .btn-row {display:flex; gap:10px; align-items:center;}
  "))),
  titlePanel("Perihilar Cholangiocarcinoma a-TO Predictor (H2O-DRF)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("age", pretty_label("age"), value = 60, min = 0, max = 120),
      
      selectInput("sex", pretty_label("sex"), choices = codebook$sex),
      selectInput("race_origin", pretty_label("race_origin"), choices = codebook$race_origin),
      selectInput("marital_status", pretty_label("marital_status"), choices = codebook$marital_status),
      selectInput("summary_stage3", pretty_label("summary_stage3"), choices = codebook$summary_stage3),
      selectInput("ajcc8_T", pretty_label("ajcc8_T"), choices = codebook$ajcc8_T),
      selectInput("ajcc8_N", pretty_label("ajcc8_N"), choices = codebook$ajcc8_N),
      selectInput("ajcc8_M", pretty_label("ajcc8_M"), choices = codebook$ajcc8_M),
      selectInput("chemotherapy", pretty_label("chemotherapy"), choices = codebook$chemotherapy),
      
      numericInput("income_10k", pretty_label("income_10k"), value = 5, min = 0),
      
      selectInput("rural_urban_code", pretty_label("rural_urban_code"), choices = codebook$rural_urban_code),
      selectInput("mets_liver", pretty_label("mets_liver"), choices = codebook$mets_liver),
      selectInput("mets_lung", pretty_label("mets_lung"), choices = codebook$mets_lung),
      selectInput("mets_bone", pretty_label("mets_bone"), choices = codebook$mets_bone),
      selectInput("mets_brain", pretty_label("mets_brain"), choices = codebook$mets_brain),
      
      div(class="btn-row",
          actionButton("go", "Predict"),
          actionButton("clear", "Clear history")
      ),
      helpText(sprintf("Note: Plot shows an approximate 95%% CI (Wilson). History keeps last %d predictions.", MAX_ROWS))
    ),
    
    mainPanel(
      class = "main-panel",
      div(class="plot-wrap",
          plotOutput("probPlot", height = "380px")
      ),
      div(class="tbl-wrap",
          h4("Prediction history (input + output)"),
          tableOutput("historyTable")
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  # store history: each row is one prediction
  hist <- reactiveVal(data.frame())
  
  observeEvent(input$clear, {
    hist(data.frame())
  })
  
  observeEvent(input$go, {
    
    # (A) Collect raw numeric codes for model
    nd <- data.frame(
      age = as.numeric(input$age),
      sex = as.numeric(input$sex),
      race_origin = as.numeric(input$race_origin),
      marital_status = as.numeric(input$marital_status),
      summary_stage3 = as.numeric(input$summary_stage3),
      ajcc8_T = as.numeric(input$ajcc8_T),
      ajcc8_N = as.numeric(input$ajcc8_N),
      ajcc8_M = as.numeric(input$ajcc8_M),
      chemotherapy = as.numeric(input$chemotherapy),
      income_10k = as.numeric(input$income_10k),
      rural_urban_code = as.numeric(input$rural_urban_code),
      mets_liver = as.numeric(input$mets_liver),
      mets_lung = as.numeric(input$mets_lung),
      mets_bone = as.numeric(input$mets_bone),
      mets_brain = as.numeric(input$mets_brain),
      check.names = FALSE
    )
    
    # ensure exact model feature order
    nd <- nd[, feature_names, drop = FALSE]
    nd <- align_levels(nd)
    
    # (B) Predict
    pred <- as.data.frame(h2o.predict(model, as.h2o(nd)))
    p1 <- as.numeric(pred$p1[1])
    ci <- wilson_ci(p1, N_EFF)
    
    # (C) Build a human-readable covariate string (using labels, not codes)
    # Map submitted codes back to labels (reverse lookup)
    rev_label <- function(var, code_chr) {
      cb <- codebook[[var]]
      if (is.null(cb)) return(code_chr)
      nm <- names(cb)[match(code_chr, cb)]
      ifelse(is.na(nm), code_chr, nm)
    }
    
    pretty_inputs <- list(
      "Age" = input$age,
      "Sex" = rev_label("sex", input$sex),
      "Race/Origin" = rev_label("race_origin", input$race_origin),
      "Marital" = rev_label("marital_status", input$marital_status),
      "Stage" = rev_label("summary_stage3", input$summary_stage3),
      "T" = rev_label("ajcc8_T", input$ajcc8_T),
      "N" = rev_label("ajcc8_N", input$ajcc8_N),
      "M" = rev_label("ajcc8_M", input$ajcc8_M),
      "Chemo" = rev_label("chemotherapy", input$chemotherapy),
      "Income($10k)" = input$income_10k,
      "RU" = rev_label("rural_urban_code", input$rural_urban_code),
      "Liver mets" = rev_label("mets_liver", input$mets_liver),
      "Lung mets"  = rev_label("mets_lung", input$mets_lung),
      "Bone mets"  = rev_label("mets_bone", input$mets_bone),
      "Brain mets" = rev_label("mets_brain", input$mets_brain)
    )
    
    cov_text <- covariate_string(pretty_inputs)
    
    # (D) Append to history
    h <- hist()
    new_row <- data.frame(
      id = if (nrow(h) == 0) 1 else max(h$id) + 1,
      p1 = p1,
      lower = ci["lower"],
      upper = ci["upper"],
      class_0.50 = ifelse(p1 >= 0.5, 1, 0),
      class_F1 = ifelse(p1 >= THRESHOLD_F1, 1, 0),
      covariates = cov_text,
      stringsAsFactors = FALSE
    )
    
    # also store the coded inputs (optional columns in table)
    for (v in feature_names) new_row[[v]] <- nd[[v]][1]
    
    h2 <- rbind(h, new_row)
    
    # keep last MAX_ROWS
    if (nrow(h2) > MAX_ROWS) h2 <- h2[(nrow(h2) - MAX_ROWS + 1):nrow(h2), ]
    
    hist(h2)
  })
  
  output$probPlot <- renderPlot({
    h <- hist()
    validate(need(nrow(h) > 0, "Click Predict to add results."))
    
    # y positions (top = newest)
    h$y <- seq_len(nrow(h))
    # color id cycles through palette
    h$col <- factor((h$id - 1) %% length(palette_vec) + 1)
    
    ggplot(h, aes(y = y)) +
      geom_segment(aes(x = lower, xend = upper, yend = y, color = col), linewidth = 1.3) +
      geom_point(aes(x = p1, color = col), size = 3) +
      geom_vline(xintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = THRESHOLD_F1, linetype = "dotted") +
      scale_color_manual(values = palette_vec, guide = "none") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      scale_y_continuous(NULL, breaks = h$y, labels = paste0("Case ", h$id)) +
      labs(
        title = "95% Confidence Interval for Response (a-TO = 1)",
        subtitle = sprintf("Dashed: 0.50 | Dotted: F1-opt = %.3f | Showing last %d predictions", THRESHOLD_F1, MAX_ROWS),
        x = "Probability"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  })
  
  output$historyTable <- renderTable({
    h <- hist()
    validate(need(nrow(h) > 0, ""))
    
    # show a compact table: id + p1/CI + classes + readable covariates
    out <- h[, c("id","p1","lower","upper","class_0.50","class_F1","covariates")]
    names(out) <- c("Case","p(a-TO=1)","CI_lower","CI_upper","Class@0.50","Class@F1","Inputs")
    out
  })
}

shinyApp(ui, server)
