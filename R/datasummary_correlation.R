cor_mat <- function(x, ...) {
  mat <- as.matrix(x)
  n <- ncol(mat)
  out_r <- out_n <- out_p <- matrix(NA, n, n)
  diag(out_r) <- 1
  diag(out_p) <- 0.00
  diag(out_n) <- colSums(!is.na(mat))

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      out_n[i, j] <- out_n[j, i] <- sum(is.finite(mat[, i]) & is.finite(mat[, j]))
      if (out_n[i, j] > 2) {
        tmp <- stats::cor.test(mat[, i], mat[, j], exact = FALSE, ...)
        out_r[i, j] <- out_r[j, i] <- tmp$estimate
        out_p[i, j] <- out_p[j, i] <- tmp$p.value
      }
    }
  }
  colnames(out_r) <- rownames(out_r) <- colnames(mat)
  colnames(out_p) <- rownames(out_p) <- colnames(mat)
  colnames(out_n) <- rownames(out_n) <- colnames(mat)
  list(r = out_r, p = out_p, n = out_n)
}

#' Generate a correlation table for all numeric variables in your dataset.
#'
#' @inheritParams datasummary
#' @param method character or function
#' \itemize{
#'   \item character: "pearson", "kendall", "spearman", or "pearspear"
#'     (Pearson correlations above and Spearman correlations below the diagonal)
#'   \item function: takes a data.frame with numeric columns and returns a
#'   square matrix with unique row.names and colnames.
#' }
#' @param bold numeric value. If not NULL (default), bold correlations with a
#'   two-sided significance below value.
#' @param add_n_note logical. If TRUE, add a note giving the range of the numbers of
#'   observations that were used to calculate correlations. Useful when
#'   the data contains missing values. Deafults to FALSE
#' @export
datasummary_correlation <- function(data,
                                    output = 'default',
                                    fmt = 2,
                                    title = NULL,
                                    notes = NULL,
                                    method = "pearson",
                                    bold = NULL,
                                    add_n_note = FALSE) {

  # sanity checks
  sanity_output(output)

  any_numeric <- any(sapply(data, is.numeric) == TRUE)
  if (any_numeric == FALSE) {
    stop("`datasummary_correlation` can only summarize numeric data columns.")
  }

  checkmate::assert(
    checkmate::check_choice(method, c("pearson", "kendall", "spearman", "pearspear")),
    checkmate::check_function(method)
  )

  if (is.function(method) & (!is.null(bold) | add_n_note)) stop(paste(
    "Parameters 'bold' and 'and_n_note' are not possible when custom method",
    "function is provided."
  ))

  # subset to numeric
  data <- data[, sapply(data, is.numeric), drop = FALSE]

  if (is.character(method)) {
    if (method == "pearspear") {
      cor_list <- cor_mat(data, method = "pearson")
      spear_cor <- cor_mat(data, method = "spearman")
      out <- cor_list$r
      out[lower.tri(out)] <- spear_cor$r[lower.tri(out)]
    } else {
      cor_list <- cor_mat(data, method = method)
      out <- cor_list$r
    }
    if(add_n_note) n_range <- c(min(cor_list$n), max(cor_list$n))
    if(!is.null(bold)) bold_mat <- cor_list$p < bold
  } else {
    out <- method(data)
  }

  if (is.function(method) && !is.matrix(out)) {
    stop("The function supplied to the `method` argument did not return a square matrix with row.names and colnames.")
  }

  out <- data.frame(out)
  out <- cbind(rowname=row.names(out), out)

  clean_r <- function(x) {
    x <- rounding(x, fmt)
    x <- gsub('0\\.', '\\.', x)
    x <- gsub('1\\.0*', '1', x)
    return(x)
  }

  for (i in seq_along(out)) {
    if (is.numeric(out[[i]])) {
      out[[i]] <- clean_r(out[[i]])
    }
  }

  if(is.function(method) || method != "pearspear") {
    for (i in 1:nrow(out)) {
      for (j in 2:ncol(out)) {
        out[i, j] <- ifelse(i + 1 < j, '.', out[i, j])
      }
    }
  }
  colnames(out) <- c(' ', out[[1]])

  align <- paste0('l', strrep('r', ncol(out) - 1))

  if(add_n_note) notes = paste(notes, sprintf(
    "Number of observations for correlations ranges from %d to %d.",
    n_range[1], n_range[2]
  ))

  factory(out,
    align = align,
    hrule = NULL,
    notes = notes,
    output = output,
    title = title)

}
