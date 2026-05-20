#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr2)
})

project_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
out_dir <- file.path(project_dir, "data", "decode")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/24_decode_token_helper.R list <TOKEN>\n",
    "  Rscript scripts/24_decode_token_helper.R download-readme <TOKEN>\n",
    "  Rscript scripts/24_decode_token_helper.R download-file <TOKEN> <FILE_KEY>\n",
    sep = ""
  )
  quit(status = 1)
}

if (length(args) < 2) usage()

mode <- args[[1]]
token <- args[[2]]

base_url <- "https://download.decode.is"

req_json <- function(path, params = list()) {
  req <- request(base_url) |>
    req_url_path_append(path) |>
    req_url_query(!!!params) |>
    req_timeout(120)
  resp <- req_perform(req)
  resp_body_json(resp, simplifyVector = TRUE)
}

download_binary <- function(path, params, dest) {
  req <- request(base_url) |>
    req_url_path_append(path) |>
    req_url_query(!!!params) |>
    req_timeout(600)
  resp <- req_perform(req)
  writeBin(resp_body_raw(resp), dest)
}

if (mode == "list") {
  x <- req_json("s3/folder", list(token = token))
  saveRDS(x, file.path(out_dir, "decode_folder_listing.rds"))

  files <- as.data.table(x$files)
  if (nrow(files)) {
    fwrite(files, file.path(out_dir, "decode_folder_listing.csv"))
    cat("Saved file listing to:\n")
    cat("  data/decode/decode_folder_listing.csv\n")
    if (!is.null(x$directoryName)) cat("Directory:", x$directoryName, "\n")
    if (!is.null(x$dlTokenValidDays)) cat("Token valid days:", x$dlTokenValidDays, "\n")
    print(files[, .(Key, Size)][1:min(.N, 20)])
  } else {
    cat("No files returned.\n")
    print(x)
  }
  quit(status = 0)
}

if (mode == "download-readme") {
  x <- req_json("s3/folder", list(token = token))
  if (is.null(x$readme$Key) || !nzchar(x$readme$Key)) {
    stop("No readme entry found in folder response.")
  }
  dest <- file.path(out_dir, basename(x$readme$Key))
  download_binary("s3/download", list(token = token, file = x$readme$Key), dest)
  cat("Downloaded readme to:\n")
  cat("  ", dest, "\n", sep = "")
  quit(status = 0)
}

if (mode == "download-file") {
  if (length(args) < 3) usage()
  file_key <- args[[3]]
  dest <- file.path(out_dir, basename(file_key))
  download_binary("s3/download", list(token = token, file = file_key), dest)
  cat("Downloaded file to:\n")
  cat("  ", dest, "\n", sep = "")
  quit(status = 0)
}

usage()

