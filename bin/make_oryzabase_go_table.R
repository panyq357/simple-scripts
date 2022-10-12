#!/usr/bin/env Rscript

OryzabaseGeneListEnURL <- "https://shigen.nig.ac.jp/rice/oryzabase/gene/download?classtag=GENE_EN_LIST"

parse_arguments <- function() {

    raw_args <- commandArgs(trailingOnly = TRUE)

    help_message = sprintf("version: 0.1.0
author: Pan Yongqing

usage: make_oryzabase_go_table.R [-h] [-o OUT_PATH]

Download annotation file from Oryzabase, extract Gene ID and GO annotation (BP),
use GO.db package to iteratively find all ancestor node, combine and export.

Suitable for subsequent GO enrichment analysis using enricher() function of clusterProfiler.

Annotation file URL: \"%s\"

options:
  -h, --help           show this help message and exit
  -o OUT_PATH          output csv filepath (default: \"oryzabase_go_table.csv\")

dependencies:
  - bioconductor-go.db
  - r-readr
  - r-tidyr
", OryzabaseGeneListEnURL)

    args_list <- list()
    option2variable <- c(
        "-o" = "out_path"
    )

    if (length(raw_args) != 0) {
        if (length(raw_args) == 0 || raw_args == "-h" || raw_args == "--help") {
            write(help_message, stderr())
            quit()
        }

        key <- ""
        for (i in 1:length(raw_args)) {
            if (raw_args[i] %in% names(option2variable)) {
                key <- option2variable[raw_args[i]]
                next
            } else {
                if (!is.null(args_list[[key]])) {
                    option_name <- names(option2variable[option2variable == key])
                    stop(sprintf("Option \"%s\" can only accept one argument!", option_name))
                }
                args_list[[key]] <- raw_args[i]
            }
        }
    }

    # Default values
    if (is.null(args_list[["out_path"]]))
        args_list[["out_path"]] <- "oryzabase_go_table.csv"

    for (var in option2variable) {
        if (is.null(args_list[[var]]))
            stop(sprintf("Missing Variable \"%s\"", toupper(var)))
    }

    return(args_list)
}

args_list <- parse_arguments()

cat("Downloading annotation data from Oryzabase ... ")
OryzabaseGeneListEn <- suppressMessages(readr::read_tsv(OryzabaseGeneListEnURL))
cat("DONE\n")

cat("Data wrangling ... ")
# Extract RAPID and GeneOntology from original big table
rapid2go <- setNames(OryzabaseGeneListEn[c("RAP ID", "Gene Ontology")], c("RAPID", "GeneOntology"))
rapid2go <- subset(rapid2go, !is.na(RAPID) & !is.na(GeneOntology))

# Split one row to multiple rows by GO terms
# (commas exist in GO term, so we have to split it using ", GO:")
rapid2go <- tidyr::separate_rows(rapid2go, GeneOntology, sep=", GO:")

# Re-prefix GO ID with "GO:"
mask_without_GO_prefix <- grep("GO:", rapid2go$GeneOntology, invert=T)
rapid2go$GeneOntology[mask_without_GO_prefix] <- paste0(
    "GO:",
    rapid2go$GeneOntology[mask_without_GO_prefix]
)

# Separate GO ID and GO term
rapid2go <- tidyr::separate(rapid2go, GeneOntology, into=c("GOID", "GOTERM"), sep=" - ")

# Add a column to indicate which domain a GO ID is belonging (BP? MF? CC?),
# using information from GO.db package
rapid2go[["ONTOLOGY"]] <- suppressMessages(
    AnnotationDbi::select(GO.db::GO.db, keys=rapid2go[["GOID"]], columns="ONTOLOGY")[["ONTOLOGY"]]
)
cat("DONE\n")

# Interatively get ancestor GO ID
all_gene2go_bp <- list()
rapid2go_splitted_by_gene <- split(rapid2go, rapid2go$RAPID)

# Progress scroll
total_gene_number <- length(rapid2go_splitted_by_gene)
counter <- 0
for (gene in names(rapid2go_splitted_by_gene)) {
    # Progress scroll
    counter <- counter + 1
    cat(sprintf("Finding ancestor GO ID: %s/%s ...\r", counter, total_gene_number))

    # Subset GO ID in BP
    rapid2go_per_gene <- rapid2go_splitted_by_gene[[gene]]
    rapid2go_per_gene_bp <- subset(rapid2go_per_gene, ONTOLOGY == "BP")
    if (nrow(rapid2go_per_gene_bp) == 0) next

    # Get ancestors
    bp_ancestor_go_id <- unlist(as.list(GO.db::GOBPANCESTOR[rapid2go_per_gene_bp[["GOID"]]]))

    # Combine with existing GO ID, and remove duplicates and root node
    all_bp_go_id <- c(rapid2go_per_gene_bp[["GOID"]], bp_ancestor_go_id)
    all_bp_go_id <- all_bp_go_id[!duplicated(all_bp_go_id)]
    all_bp_go_id <- all_bp_go_id[all_bp_go_id != "all"]

    names(all_bp_go_id) <- NULL
    all_gene2go_bp[[gene]] <- data.frame(RAPID=gene, GOID=all_bp_go_id)
}
cat(sprintf("Finding ancestor GO ID: %s/%s ... DNOE\n", counter, total_gene_number))

# Combine, fill-in GO term and export
cat(sprintf("Exporting table to \"%s\" ... ", args_list[["out_path"]]))
all_gene2go_bp <- do.call(rbind, all_gene2go_bp)
all_gene2go_bp[["GOTERM"]] <- suppressMessages(
    AnnotationDbi::select(GO.db::GO.db, keys=all_gene2go_bp[["GOID"]], columns="TERM")[["TERM"]]
)
readr::write_csv(all_gene2go_bp, args_list[["out_path"]])
cat("DONE\n")
