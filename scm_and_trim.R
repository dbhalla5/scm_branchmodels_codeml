#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(ape)
  library(phytools)
  library(geiger)
  library(readr)
  library(dplyr)
})

###  How to use the script (2 modes)

### ---------------------------------------  ### 

##   1) Generate full-tree SCMs (once per trait)

# Rscript scm_and_trim.R --mode scm \
#  --trait Primary.Lifestyle \
#  --baseline Terrestrial \
#  --nsim 50

##   2) Trim for a specific CYP gene (as needed)

# Rscript scm_and_trim.R --mode trim \
#  --trait Primary.Lifestyle \
#  --gene CYP4A5 \
#  --cyp_tree CYP4A5_subset_tree.nh \
#  --take_n 20 \
#  --random


# ----------------------------- USER PATHS -----------------------------
FULL_TREE <- "/path/to/SCM/363-avian-2020-phast_shortname.nh"
TRAIT_TABLE <- "/path/to/SCM/avian_trait_mapping.csv"

# Output base
OUT_ROOT <- "/path/to/SCM"
# ---------------------------------------------------------------------

# ------------------------------ HELPERS -------------------------------
fit_best_mk <- function(tree, states) {
  models <- c("ER","SYM","ARD")
  fits <- list(); AICs <- setNames(rep(Inf, length(models)), models)

  for (m in models) {
    cat("[INFO] Fitting", m, "model...\n")
    fi <- tryCatch(phytools::fitMk(tree, states, model = m), error = function(e) NULL)
    if (!is.null(fi)) {
      fits[[m]] <- fi
      AICs[m] <- AIC(fi)
    } else {
      cat("[WARN]", m, "model failed; skipping.\n")
    }
  }
  cat("[INFO] AIC values:\n"); print(AICs)
  best <- names(which.min(AICs))
  if (!is.finite(AICs[best])) stop("All Mk fits failed.")
  list(best_model = best, best_fit = fits[[best]], AICs = AICs)
}

construct_Q_if_needed <- function(best_fit, best_model, k) {
  constructed <- FALSE
  if (is.null(best_fit$Q)) {
    constructed <- TRUE
    r <- if (!is.null(best_fit$rates)) as.numeric(best_fit$rates) else 1.0
    Q <- matrix(0, k, k)
    Q[row(Q) != col(Q)] <- r
    diag(Q) <- -rowSums(Q)
  } else {
    Q <- best_fit$Q
  }
  list(Q = Q, constructed = constructed)
}

label_for_paml <- function(tree, edge_states, bg, class_ids) {
  # Append '#k' to non-baseline branches (tips or internal nodes)
  if (is.null(tree$node.label)) tree$node.label <- rep("", tree$Nnode)
  for (i in seq_len(nrow(tree$edge))) {
    child <- tree$edge[i, 2]
    st <- edge_states[i]
    if (!is.na(st) && st != bg) {
      tag <- paste0("#", class_ids[[as.character(st)]])
      if (child <= length(tree$tip.label)) {
        tree$tip.label[child] <- paste0(tree$tip.label[child], tag)
      } else {
        idx <- child - length(tree$tip.label)
        tree$node.label[idx] <- paste0(tree$node.label[idx], tag)
      }
    }
  }
  tree
}

read_states_from_table <- function(tr, trait_col) {
  # Expect a "short_name" column that matches tree tip labels (e.g., APTOWE)
  tab <- read_csv(TRAIT_TABLE, show_col_types = FALSE)

  if (!("short_name" %in% names(tab))) stop("Trait table lacks 'short_name' column.")
  if (!(trait_col %in% names(tab))) stop("Trait table lacks '", trait_col, "' column.")

  df <- tab %>% select(short_name, !!trait_col)
  names(df) <- c("Species","Trait")

  df$Species <- gsub(" ", "_", df$Species)
  # Keep only species present in tree
  keep <- intersect(tr$tip.label, df$Species)
  if (length(keep) == 0) stop("No species names match between tree and table.")
  missing <- setdiff(tr$tip.label, df$Species)
  if (length(missing)) cat("[WARN]", length(missing), "tree tips lack trait; they will be dropped in SCM.\n")

  states <- setNames(df$Trait, df$Species)[tr$tip.label]
  states <- droplevels(factor(states))
  states
}

clean_tree_for_mk <- function(tr) {
  if (!ape::is.binary(tr)) {
    cat("[WARN] Resolving polytomies with multi2di().\n")
    tr <- multi2di(tr, random = TRUE)
  }
  if (!is.null(tr$edge.length)) {
    z <- which(tr$edge.length <= 1e-8)
    if (length(z)) {
      tr$edge.length[z] <- 1e-6
      cat("[WARN] Adjusted", length(z), "zero-length branches.\n")
    }
  }
  if (!is.rooted(tr)) {
    cat("[WARN] Midpoint-rooting tree for SCM.\n")
    tr <- phytools::midpoint.root(tr)
  }
  tr
}
# ---------------------------------------------------------------------

# ------------------------------ ARGPARSE ------------------------------
parser <- ArgumentParser(description = "SCM on full tree and trimming to CYP subsets for codeml.")
parser$add_argument("--mode", required = TRUE, choices = c("scm","trim"),
                    help = "scm = run SCM on full tree; trim = trim existing SCMs to a CYP subset")
parser$add_argument("--trait", required = TRUE, help = "Trait column name (e.g., Primary.Lifestyle)")
parser$add_argument("--baseline", help = "Baseline trait state (for SCM labeling; left unlabeled)")
parser$add_argument("--nsim", type = "integer", default = 50, help = "Number of SCM replicates (default 50)")

# TRIM mode args
parser$add_argument("--gene", help = "CYP gene name for trimming (e.g., CYP4A5)")
parser$add_argument("--cyp_tree", help = "Path to CYP subset tree (.nh) for trimming")
parser$add_argument("--take_n", type = "integer", default = 20, help = "How many SCM trees to trim (default 20)")
parser$add_argument("--random", action = "store_true", help = "Randomly sample SCM trees (default: take first N)")

args <- parser$parse_args()
MODE <- args$mode
TRAIT <- args$trait
BASELINE <- args$baseline
NSIM <- args$nsim

cat("[INFO] Mode:", MODE, "\n")
cat("[INFO] Trait:", TRAIT, "\n")

# ------------------------------- SCM MODE ----------------------------
if (MODE == "scm") {
  if (is.null(BASELINE)) stop("--baseline is required in --mode scm")

  if (!file.exists(FULL_TREE)) stop("Full tree not found: ", FULL_TREE)
  tr0 <- read.tree(FULL_TREE)
  cat("[INFO] Full-tree tips:", length(tr0$tip.label), "\n")

  states0 <- read_states_from_table(tr0, TRAIT)
  cat("[INFO] Trait levels:", paste(levels(states0), collapse = ", "), "\n")
  cat("[INFO] State counts:\n"); print(table(states0))
  if (nlevels(states0) < 2) stop("Only one trait state present; SCM not possible.")

  # Drop tips lacking trait
  keep <- names(states0)[!is.na(states0)]
  tr <- drop.tip(tr0, setdiff(tr0$tip.label, keep))
  states <- droplevels(states0[tr$tip.label])

  # Clean tree for Mk fitting
  cat("[INFO] Cleaning tree...\n")
  tr <- clean_tree_for_mk(tr)

  # Fit Mk models and pick best by AIC
  mk <- fit_best_mk(tr, states)
  best_model <- mk$best_model
  best_fit <- mk$best_fit
  AICs <- mk$AICs

  # Ensure Q
  k <- nlevels(states)
  qc <- construct_Q_if_needed(best_fit, best_model, k)
  Q <- as.matrix(qc$Q)
  constructed_flag <- qc$constructed
  levs <- levels(states)
  rownames(Q) <- colnames(Q) <- levs

  # Output dirs
  OUT_BASE <- file.path(OUT_ROOT, "full_tree", TRAIT)
  TREES_DIR <- file.path(OUT_BASE, "trees_paml")
  dir.create(TREES_DIR, recursive = TRUE, showWarnings = FALSE)

  # Save model info and Q
  model_info_path <- file.path(OUT_BASE, sprintf("%s_best_model.txt", TRAIT))
  writeLines(
    c(
      paste("Trait:", TRAIT),
      paste("Best_Model:", best_model),
      paste("AICs:", paste(names(AICs), round(AICs, 2), collapse = " | ")),
      paste("Constructed_Q:", constructed_flag),
      paste("Date:", Sys.time())
    ),
    con = model_info_path
  )
  cat("[INFO] Saved model info to", model_info_path, "\n")

  Q_suffix <- if (constructed_flag) "constructed" else "estimated"
  Q_csv_path <- file.path(OUT_BASE, sprintf("%s_Qmatrix_%s_%s.csv", TRAIT, best_model, Q_suffix))
  write.csv(round(Q, 6), Q_csv_path, row.names = TRUE)
  cat("[INFO] Saved Q matrix to", Q_csv_path, "\n")

  # Baseline & class labels
  if (!(BASELINE %in% levs)) {
    stop(paste0("Baseline '", BASELINE, "' not among trait levels: ", paste(levs, collapse = ", ")))
  }
  class_ids <- setNames(seq_along(levs), levs)
  writeLines(sprintf("%s->#%d", levs, class_ids), con = file.path(TREES_DIR, "class_labels.txt"))
  cat("[INFO] PAML classes:", paste(sprintf("%s->#%d", levs, class_ids), collapse=", "), "\n")

  # SCM
  set.seed(2000)
  scm_rds <- file.path(OUT_BASE, sprintf("%s_simmap_%03d.rds", TRAIT, NSIM))
  cat("[INFO] Running make.simmap (best model =", best_model, "), nsim =", NSIM, "\n")
  simlist <- phytools::make.simmap(tr, states, model = best_model, nsim = NSIM, Q = Q, pi = "estimated")
  saveRDS(simlist, scm_rds)
  cat("[INFO] Saved SCMs:", scm_rds, "\n")

  # Write PAML-tagged full trees (rooting kept; trimming step will unroot)
  for (i in seq_len(NSIM)) {
    sm <- simlist[[i]]
    edge_states <- tryCatch(phytools::getStates(sm, type = "edges"),
                            error = function(e) {
                              if (!is.null(sm$maps)) sapply(sm$maps, function(x) names(x)[length(x)]) else stop(e)
                            })
    tagged <- label_for_paml(sm, edge_states, bg = BASELINE, class_ids = class_ids)
    fn <- file.path(TREES_DIR, sprintf("%s_SCM_%03d_paml.nh", TRAIT, i))
    write.tree(tagged, file = fn)
    cat("[INFO] Wrote", fn, "\n")
  }
  cat("[INFO] Done. Full-tree SCM Newicks in:", TREES_DIR, "\n")
  quit(save="no")
}

# ------------------------------ TRIM MODE ----------------------------
if (MODE == "trim") {
  if (is.null(args$gene) || is.null(args$cyp_tree)) {
    stop("--gene and --cyp_tree are required in --mode trim")
  }
  GENE <- args$gene
  CypTreePath <- args$cyp_tree
  TAKE_N <- args$take_n
  RANDOM <- isTRUE(args$random)

  # Where to read full-tree SCMs
  FULL_TREES_DIR <- file.path(OUT_ROOT, "full_tree", TRAIT, "trees_paml")
  if (!dir.exists(FULL_TREES_DIR)) stop("Full-tree SCM directory not found: ", FULL_TREES_DIR)

  class_map_file <- file.path(FULL_TREES_DIR, "class_labels.txt")
  if (!file.exists(class_map_file)) {
    cat("[WARN] class_labels.txt not found; proceeding without explicit mapping file.\n")
  }

  # Where to write trimmed trees
  OUT_BASE <- file.path(OUT_ROOT, GENE, TRAIT)
  TRIM_DIR <- file.path(OUT_BASE, "trees_paml_trimmed")
  dir.create(TRIM_DIR, recursive = TRUE, showWarnings = FALSE)

  # Load CYP subset tree (use the same one you feed to codeml site models)
  if (!file.exists(CypTreePath)) stop("CYP tree not found: ", CypTreePath)
  cyp_tr <- read.tree(CypTreePath)
  cat("[INFO] CYP tree tips:", length(cyp_tr$tip.label), "\n")

  # Collect available full SCM trees
  all_scm <- list.files(FULL_TREES_DIR, pattern = paste0("^", TRAIT, "_SCM_\\d+_paml\\.nh$"), full.names = TRUE)
  if (length(all_scm) == 0) stop("No full-tree SCM Newicks found for trait: ", TRAIT)

  if (RANDOM) {
    set.seed(2025)
    sel <- sample(all_scm, min(TAKE_N, length(all_scm)))
  } else {
    sel <- head(all_scm, min(TAKE_N, length(all_scm)))
  }

  writeLines(basename(sel), con = file.path(TRIM_DIR, "selected_full_scm_trees.txt"))
  cat("[INFO] Selected", length(sel), "SCM trees. List saved to selected_full_scm_trees.txt\n")

  # Trim and unroot
  for (f in sel) {
    tr <- read.tree(f)

    # Intersect tips with CYP tree
    keep <- intersect(cyp_tr$tip.label, tr$tip.label)
    if (length(keep) == 0) {
      cat("[WARN] No overlapping tips between", basename(f), "and CYP tree; skipping.\n")
      next
    }
    drop_from_tr <- setdiff(tr$tip.label, keep)
    if (length(drop_from_tr)) {
      tr_trim <- drop.tip(tr, drop_from_tr)
    } else {
      tr_trim <- tr
    }

    # Unroot for codeml
    tr_trim <- unroot(tr_trim)

    out_fn <- file.path(TRIM_DIR,
                        sprintf("%s_%s.nh",
                                GENE,
                                sub("\\.nh$", "", basename(f))))
    write.tree(tr_trim, file = out_fn)
    cat("[INFO] Wrote trimmed/unrooted:", out_fn, "\n")
  }

  cat("[INFO] Done. Trimmed codeml-ready trees in:", TRIM_DIR, "\n")
  quit(save="no")
}

stop("Unknown mode: ", MODE)

