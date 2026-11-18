# =============================================================================
# IMPORT FUNCTIONS — preserve the biological identity of every measurement
# =============================================================================
# Public functions:
#   read_long(path, kind): a CSV already arranged as one measurement per row.
#   read_blocks(path, sheet, kind): an Excel sheet with stacked table blocks.
# Helpers perform strict numeric/date conversion and parse the blocks.
#
# kind is either "recovery" or "inoculum". The returned count column is named
# recovery_cfu_ml or initial_cfu_ml respectively. Counts must already be CFU/ml;
# these functions do not calculate dilution factors or plated-volume corrections.
# Invalid cells stop import rather than becoming an unnoticed NA.
# =============================================================================
# Parse raw worksheet cells, not a sheet whose first row became column names.
# This preserves each block's own identifiers and avoids Excel name repair.
# Convert text to numbers while retaining genuinely empty cells as missing.
# A cell such as "<5" contains censoring information and must not become 5.
# A numeric-looking cell such as "Inf" is also invalid as a measured count.
strict_numeric <- function(x, label) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  out <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & (!is.finite(out) | is.na(out)))) {
    stop(label, ": nonnumeric or non-finite cell found; inspect the source.")
  }
  out
}

# Dates are accepted in only two unambiguous formats: an Excel serial number
# using the 1900 system, or ISO YYYY-MM-DD. No day/month ordering is guessed.
# Return a Date vector so joins compare actual dates, not spreadsheet formatting.
parse_dates <- function(x) {
  x <- trimws(as.character(x))
  if (anyNA(x) || any(x == "")) stop("Missing experimental date.")
  out <- rep(as.Date(NA), length(x))
# Split serial-number cells from ISO text before conversion. Numeric Excel
  # dates count days from a historical origin; a decimal part represents time
  # of day, which is not retained by this date-level workflow.
  serial <- grepl("^[0-9]+(\\.[0-9]+)?$", x)
  out[serial] <- as.Date(as.numeric(x[serial]), origin = "1899-12-30")
  iso <- !serial & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)
  out[iso] <- as.Date(x[iso], format = "%Y-%m-%d")
  if (anyNA(out)) stop("Dates must be Excel serial dates (1900 system) or ISO YYYY-MM-DD.")
  out
}

# Parse a raw sheet supplied as a matrix or data frame, including all headers.
# Recovery layout: Plate | Date | Repeat | isolate IDs ... | optional Average.
# Inoculum layout: Date | Repeat | sampled isolate IDs ... | optional Average.
# Each later block has its own header and may contain different isolates.
# The output is a long data frame; source_block records the original table.
parse_blocks <- function(raw, kind = c("recovery", "inoculum")) {
  kind <- match.arg(kind)
  x <- as.matrix(raw)
  storage.mode(x) <- "character"
  x[] <- trimws(x)
  x[x == ""] <- NA_character_
# dc/rc are the date/repeat column positions, not experimental metadata.
  # Inoculum has no plate column, so both positions shift one place left.
  dc <- if (kind == "recovery") 2L else 1L
  rc <- dc + 1L
  if (ncol(x) <= rc) stop("Insufficient columns for ", kind, " blocks.")
# Detect actual block headers from their Date and Repeat labels. A blank row
  # alone is not a reliable block boundary because spreadsheets may have several.
  starts <- which(!is.na(x[, dc]) & !is.na(x[, rc]) &
                    tolower(x[, dc]) == "date" & tolower(x[, rc]) == "repeat")
  if (!length(starts)) stop("No Date / Repeat block headers found.")
  # Reject unexpected material before the first header.
  if (starts[1] > 1L && any(!is.na(x[seq_len(starts[1] - 1L), , drop = FALSE]))) {
    stop("Nonblank rows precede the first header; remove titles or notes.")
  }
# A block ends immediately before the next header (or at the end of the sheet).
  # Blank separator rows inside these boundaries are removed later.
  ends <- c(starts[-1] - 1L, nrow(x))
  blocks <- lapply(seq_along(starts), function(b) {
    if (ends[b] <= starts[b]) stop("Header has no measurement rows in block ", b)
# Remove inherited vector names: these cells are header VALUES, not a named
    # selection map. This also avoids tidyselect interpreting IDs as renaming.
    h <- unname(x[starts[b], ])
# All columns after Repeat are potential isolates. Mean/Average are summaries,
    # not isolates; including them would create artificial biological observations.
    candidate <- seq.int(rc + 1L, ncol(x))
    keep <- candidate[!is.na(h[candidate]) &
                        !grepl("^(average|mean)$|^\\.\\.\\.", h[candidate], ignore.case = TRUE)]
    if (!length(keep)) stop("No isolate headers in block ", b)
    ids <- h[keep]
    if (anyDuplicated(ids)) stop("Duplicate isolate header in block ", b)
# Retain nonblank rows, then insist that each is a recognisable repeat.
    # Notes, merged labels and unexpected material must be resolved explicitly.
    rows <- seq.int(starts[b] + 1L, ends[b])
    active <- rows[rowSums(!is.na(x[rows, , drop = FALSE])) > 0L]
    if (!length(active)) stop("No observations in block ", b)
    labels <- x[active, rc]
    if (anyNA(labels) || any(!grepl("^R[0-9]+$", labels))) {
      stop("Nonblank nonmeasurement row in block ", b, "; expected R1, R2, ...")
    }
# A repeat must occur only once within a block. Multiple rows for R1 might
    # represent technical replicates, not extra independent biological repeats.
    if (anyDuplicated(labels)) stop("Duplicate replicate row in block ", b)
    # Counts under blank/unnamed headers must never disappear silently.
    badcols <- candidate[is.na(h[candidate]) | grepl("^\\.\\.\\.", h[candidate])]
    if (length(badcols) && any(!is.na(x[active, badcols, drop = FALSE]))) {
      stop("Data under a missing/unnamed isolate header in block ", b)
    }
    dates <- parse_dates(x[active, dc])
    plate <- if (kind == "recovery") h[1] else NA_character_
    if (kind == "recovery" && is.na(plate)) stop("Missing plate label in block ", b)
# The worksheet matrix is row-wise in biological meaning: one repeat contains
    # many isolates. t() followed by as.vector() preserves that row-wise order.
    # rep(..., each = number of isolates) then gives every count the correct date
    # and repeat; rep(ids, times = number of rows) gives it the correct isolate.
    vals <- strict_numeric(as.vector(t(x[active, keep, drop = FALSE])),
                           paste(kind, "block", b))
    out <- data.frame(
      source_block = b,
      date = rep(dates, each = length(ids)),
      replicate = rep(labels, each = length(ids)),
      isolate_id = rep(ids, times = length(active)),
      stringsAsFactors = FALSE
    )
# Inoculum blocks have no plate assignment in the source. Do not invent one
    # from the block number; later matching uses actual IDs, dates and repeats.
    if (kind == "recovery") out$plate <- plate
    out[[if (kind == "recovery") "recovery_cfu_ml" else "initial_cfu_ml"]] <- vals
    out
  })
  dplyr::bind_rows(blocks)
}

# Read every cell as text with col_names = FALSE. Otherwise the first block
# header becomes R column names while later headers remain data, increasing
# the risk that later counts inherit the wrong isolate IDs. Hidden rows are
# still imported; the user must define the intended input region beforehand.
read_blocks <- function(path, sheet, kind) {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                           col_types = "text", .name_repair = "minimal")
  parse_blocks(raw, kind)
}

# CSV route for an already tidy table. Read IDs as character to preserve
# leading zeros. Keep only the documented columns so unrelated private metadata
# are not inadvertently carried into exported tables.
read_long <- function(path, kind) {
  d <- utils::read.csv(path, colClasses = "character", check.names = FALSE,
                       na.strings = c("", "NA"))
  value <- if (kind == "recovery") "recovery_cfu_ml" else "initial_cfu_ml"
  required <- c("isolate_id", "date", "replicate", value,
                if (kind == "recovery") "plate")
  if (!all(required %in% names(d))) stop("Missing columns in ", path)
  d <- d[, required, drop = FALSE]
  d$date <- parse_dates(d$date)
  d[[value]] <- strict_numeric(d[[value]], value)
  d
}
