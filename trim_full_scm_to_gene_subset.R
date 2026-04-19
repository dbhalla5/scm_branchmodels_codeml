#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

##  USAGE 5 arguments to be supplied

# Rscript trim_full_scm_to_gene.R Trophic.Level CYP4A5 \
#  /path/to/SCM/full_tree/Trophic.Level/trees_paml \
#  /path/to/SCM/sub_trees/CYP4A5.phy \
#  /path/to/SCM/CYP4A5/Trophic.Level/trees_paml_trimmed

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage:\n",
      "  Rscript trim_full_scm_to_gene.R <trait> <gene> <full_scm_dir> <alignment.phy> [<out_dir>]\n",
      "Example:\n",
      "  Rscript trim_full_scm_to_gene.R Trophic.Level CYP4A5 \\\n",
      "    /path/to/SCM/full_tree/Trophic.Level/trees_paml \\\n",
      "    /path/to/SCM/sub_trees/CYP4A5.phy \\\n",
      "    /path/to/SCM/CYP4A5/Trophic.Level/trees_paml_trimmed\n",
      sep = "")
  quit(status = 1)
}

TRAIT    <- args[1]
GENE     <- args[2]
SCM_DIR  <- args[3]
PHY_FILE <- args[4]
OUT_DIR  <- ifelse(length(args) >= 5, args[5],
                   file.path(dirname(dirname(SCM_DIR)), GENE, TRAIT, "trees_paml_trimmed"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[INFO] Trait: ", TRAIT, "\n", sep = "")
cat("[INFO] Gene:  ", GENE, "\n", sep = "")
cat("[INFO] Full-tree SCM dir: ", SCM_DIR, "\n", sep = "")
cat("[INFO] Alignment (PHYLIP): ", PHY_FILE, "\n", sep = "")
cat("[INFO] Output dir: ", OUT_DIR, "\n", sep = "")

# ---------- helper functions ----------
norm_name <- function(x) {
  x <- gsub("#.*$", "", x)    # remove PAML #tags
  x <- gsub("\\s+", "_", x)
  toupper(x)
}

# robust PHYLIP taxa reader
read_phy_taxa <- function(path) {
  con <- file(path, "r"); on.exit(close(con))
  header <- readLines(con, n = 1, warn = FALSE)
  m <- regexec("^\\s*([0-9]+)\\s+([0-9]+)", header)
  mm <- regmatches(header, m)[[1]]
  if (length(mm) < 3) stop("Bad PHYLIP header in: ", path)
  nseq <- as.integer(mm[2]); seqlen <- as.integer(mm[3])

  taxa <- character(0)
  i <- 0
  while (i < nseq) {
    line <- ""
    repeat {
      ln <- readLines(con, n = 1, warn = FALSE)
      if (length(ln) == 0) break
      if (nzchar(trimws(ln))) { line <- ln; break }
    }
    if (!nzchar(line)) break
    sp <- strsplit(trimws(line), "\\s+")[[1]][1]
    taxa <- c(taxa, sp)

    # read sequence lines
    got <- 0
    repeat {
      ln <- readLines(con, n = 1, warn = FALSE)
      if (length(ln) == 0) break
      s <- gsub("\\s+", "", ln)
      got <- got + nchar(s)
      if (got >= seqlen) break
    }
    i <- i + 1
  }
  unique(taxa)
}

# ---------- read alignment taxa ----------
phy_taxa_raw <- read_phy_taxa(PHY_FILE)
phy_taxa_norm <- norm_name(phy_taxa_raw)
n_phy <- length(phy_taxa_norm)
cat("[INFO] Alignment taxa: ", n_phy, "\n", sep = "")

# ---------- list SCM trees ----------
trees <- list.files(SCM_DIR, pattern = "_paml\\.nh$", full.names = TRUE)
if (!length(trees)) stop("No *_paml.nh found in ", SCM_DIR)
cat("[INFO] Found ", length(trees), " full-tree SCMs.\n", sep = "")

# ---------- process trees ----------
summary <- data.frame(Tree = basename(trees),
                      Tips_Before = NA_integer_,
                      Tips_Kept = NA_integer_,
                      Duplicates_Removed = NA_integer_,
                      Missing_Vs_Alignment = NA_integer_,
                      stringsAsFactors = FALSE)

for (k in seq_along(trees)) {
  tf <- trees[k]
  tr <- read.tree(tf)

  tips_raw <- tr$tip.label
  tips_norm <- norm_name(tips_raw)

  keep_idx <- which(tips_norm %in% phy_taxa_norm)
  drop_idx <- setdiff(seq_along(tips_norm), keep_idx)

  tr_trim <- if (length(drop_idx)) drop.tip(tr, drop_idx) else tr

  #  NEW: remove duplicated tips
  dup_tips <- which(duplicated(tr_trim$tip.label))
  if (length(dup_tips)) {
    cat("[FIX] ", basename(tf), ": found ", length(dup_tips),
        " duplicated tip(s): ", paste(tr_trim$tip.label[dup_tips], collapse = ", "), "\n", sep = "")
    tr_trim <- drop.tip(tr_trim, dup_tips)
  }

  # NEW: check missing taxa
  missing_in_tree <- setdiff(phy_taxa_norm, norm_name(tr_trim$tip.label))
  if (length(missing_in_tree)) {
    cat("[ERROR] ", basename(tf), ": missing ", length(missing_in_tree),
        " taxa from alignment.\n", sep = "")
    cat("Missing taxa:\n", paste(missing_in_tree, collapse = ", "), "\n")
    stop("Aborting to avoid running codeml with incomplete trees.")
  }

  # clean, bifurcate & unroot for codeml
  tr_trim <- multi2di(tr_trim)
  if (!is.null(tr_trim$edge.length)) {
    z <- which(tr_trim$edge.length <= 1e-12)
    if (length(z)) tr_trim$edge.length[z] <- 1e-6
  }
  if (length(tr_trim$tip.label) >= 3) tr_trim <- unroot(tr_trim)

  n_kept <- length(tr_trim$tip.label)
  out_file <- file.path(OUT_DIR, basename(tf))
  write.tree(tr_trim, file = out_file)

  summary[k, ] <- c(basename(tf),
                    length(tr$tip.label),
                    n_kept,
                    length(dup_tips),
                    0)

  cat(sprintf("[OK] %s  before=%d  kept=%d  dups_removed=%d  missing_vs_alignment=0\n",
              basename(tf), length(tr$tip.label), n_kept, length(dup_tips)))
}

sum_path <- file.path(OUT_DIR, paste0(GENE, "_", TRAIT, "_trim_summary.csv"))
write.csv(summary, sum_path, row.names = FALSE)
cat("[INFO] Summary written to: ", sum_path, "\n", sep = "")

