# version: 0.1.0
# author: panyq
# usage: source this file i R programming enviroment.

# Query biomaRt for gene models in a given genomic region, and parse the query
# to return a Gviz compatible data.frame.
get_gviz_df_from_biomaRt <- function(mart, chromosome, start, end) {

    biomart_query_result <- biomaRt::getBM(
        mart=mart,
        filters=c("chromosome_name", "start", "end"),
        values=list(chromosome, start, end),
        attributes=c(
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "transcript_length",
            "exon_chrom_start",
            "exon_chrom_end",
            "5_utr_start",
            "5_utr_end",
            "3_utr_start",
            "3_utr_end",
            "strand"
        )
    )

    strand_number2symbol <- function(number) {
        symbol <- list()
        for (i in number) {
            if (i == "1") symbol <- append(symbol, "+")
            else if (i == "-1") symbol <- append(symbol, "-")
            else stop(sprintf("unknown strand number: %s", number))
        }
        return(unlist(symbol))
    }


    df_for_gviz <- data.frame(matrix(nrow=0, ncol=8))
    names(df_for_gviz) <- c("chromosome", "start", "end", "width", "strand", "feature", "gene", "transcript")

    for (gene_id in unique(biomart_query_result[["ensembl_gene_id"]])) {

        # Subset the query reult, keep each gene's longest transcript's model info.
        per_gene_query_result <- subset(biomart_query_result, ensembl_gene_id == gene_id)
        transcript_and_length <- per_gene_query_result[c("ensembl_transcript_id", "transcript_length")]
        transcript_and_length <- subset(transcript_and_length, !duplicated(transcript_and_length))
        order_of_transcript_by_length <- order(transcript_and_length[["transcript_length"]], decreasing=T)
        longest_transcript <- transcript_and_length[["ensembl_transcript_id"]][order_of_transcript_by_length[1]]
        per_gene_query_result <- subset(per_gene_query_result, ensembl_transcript_id == longest_transcript)

        # Set a empty df for adding features
        per_gene_df_for_gviz <- data.frame(matrix(nrow=0, ncol=6))
        names(per_gene_df_for_gviz) <- c("chromosome", "start", "end", "width", "strand", "feature")

        # there is some rows in query result that contain both exon and utr infomation,
        # utr's range is a subset of exon, so we have to manually split this marginal 
        # exon into utr part and cds part, and add them to per_gene_df_for_gviz

        # add 5'utr containing exon's info
        utr5_rows <- subset(per_gene_query_result, !is.na(`5_utr_start`) & !is.na(`5_utr_end`))
        if (nrow(utr5_rows) > 0) {
            utr5_exon <- IRanges::IRanges(start=utr5_rows[["exon_chrom_start"]],end=utr5_rows[["exon_chrom_end"]])
            utr5_utr <- IRanges::IRanges(start=utr5_rows[["5_utr_start"]],end=utr5_rows[["5_utr_end"]])
            utr5_cds <- IRanges::setdiff(utr5_exon,  utr5_utr)
            per_gene_df_for_gviz <- rbind(
                per_gene_df_for_gviz,
                data.frame(
                    chromosome=chromosome,
                    start=IRanges::start(utr5_utr),
                    end=IRanges::end(utr5_utr),
                    width=IRanges::width(utr5_utr),
                    strand=strand_number2symbol(utr5_rows[["strand"]]),
                    feature="utr5"
                )
            )
            if (length(utr5_cds) > 0) {
                per_gene_df_for_gviz <- rbind(
                    per_gene_df_for_gviz,
                    data.frame(
                        chromosome=chromosome,
                        start=IRanges::start(utr5_cds),
                        end=IRanges::end(utr5_cds),
                        width=IRanges::width(utr5_cds),
                        strand=strand_number2symbol(utr5_rows[["strand"]]),
                        feature="protein_coding"
                    )
                )
            }
        }

        # add 3'utr containing exon's info
        utr3_rows <- subset(per_gene_query_result, !is.na(`3_utr_start`) & !is.na(`3_utr_end`))
        if (nrow(utr3_rows) > 0) {
            utr3_exon <- IRanges::IRanges(start=utr3_rows[["exon_chrom_start"]],end=utr3_rows[["exon_chrom_end"]])
            utr3_utr <- IRanges::IRanges(start=utr3_rows[["3_utr_start"]],end=utr3_rows[["3_utr_end"]])
            utr3_cds <- IRanges::setdiff(utr3_exon,  utr3_utr)
            per_gene_df_for_gviz <- rbind(
                per_gene_df_for_gviz,
                data.frame(
                    chromosome=chromosome,
                    start=IRanges::start(utr3_utr),
                    end=IRanges::end(utr3_utr),
                    width=IRanges::width(utr3_utr),
                    strand=strand_number2symbol(utr3_rows[["strand"]]),
                    feature="utr3"
                )
            )
            if (length(utr3_cds) > 0) {
                per_gene_df_for_gviz <- rbind(
                    per_gene_df_for_gviz,
                    data.frame(
                        chromosome=chromosome,
                        start=IRanges::start(utr3_cds),
                        end=IRanges::end(utr3_cds),
                        width=IRanges::width(utr3_cds),
                        strand=strand_number2symbol(utr3_rows[["strand"]]),
                        feature="protein_coding"
                    )
                )
            }
        }

        # now, add normal exon rows to per_gene_df_for_gviz
        normal_rows <- subset(per_gene_query_result, is.na(`3_utr_start`) & is.na(`5_utr_start`))
        if (nrow(normal_rows) > 0) {
            normal_ir <- IRanges::IRanges(start=normal_rows[["exon_chrom_start"]],end=normal_rows[["exon_chrom_end"]])
            per_gene_df_for_gviz <- rbind(
                per_gene_df_for_gviz,
                data.frame(
                    chromosome=chromosome,
                    start=IRanges::start(normal_ir),
                    end=IRanges::end(normal_ir),
                    width=IRanges::width(normal_ir),
                    strand=strand_number2symbol(normal_rows[["strand"]]),
                    feature="protein_coding"
                )
            )
        }

        # add gene and transcript ID and append to overall df_for_gviz
        per_gene_df_for_gviz[["gene"]] <- gene_id
        per_gene_df_for_gviz[["transcript"]] <- longest_transcript

        df_for_gviz <- rbind(df_for_gviz, per_gene_df_for_gviz)
    }
    
    return(df_for_gviz)
}

