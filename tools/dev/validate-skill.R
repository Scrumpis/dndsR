# Validates the end-user skill. Base R only.
#   Rscript --vanilla tools/dev/validate-skill.R [repo-root]
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
skill_dir <- file.path(root, "skills", "dndsr")
errors <- character()
fail <- function(...) errors <<- c(errors, paste0(...))

skill_path <- file.path(skill_dir, "SKILL.md")
if (!file.exists(skill_path)) stop("Missing ", skill_path, call. = FALSE)
lines <- readLines(skill_path, warn = FALSE)

# Frontmatter: agentskills.io constraints. Top-level "key: value" lines only.
if (!identical(lines[1], "---")) fail("SKILL.md: frontmatter must start on line 1")
end <- which(lines == "---")[2]
if (is.na(end)) stop("SKILL.md: unterminated frontmatter", call. = FALSE)
front <- lines[2:(end - 1)]
field <- function(key, text = front) {
  hit <- grep(paste0("^\\s*", key, ":\\s*"), text, value = TRUE)
  if (!length(hit)) return(NA_character_)
  gsub("^[\"']|[\"']$", "", sub(paste0("^\\s*", key, ":\\s*"), "", hit[1]))
}
name <- field("name")
description <- field("description")
compatibility <- field("compatibility")
if (is.na(name) || !grepl("^[a-z0-9]+(-[a-z0-9]+)*$", name) || nchar(name) > 64) {
  fail("name: must be 1-64 lowercase letters, digits, or single hyphens")
}
if (!identical(name, basename(skill_dir))) fail("name: must match the directory name")
if (is.na(description) || !nzchar(description) || nchar(description) > 1024) {
  fail("description: must be 1-1024 characters")
}
if (!is.na(compatibility) && nchar(compatibility) > 500) fail("compatibility: over 500 characters")
if (length(lines) > 500) fail("SKILL.md: over the recommended 500 lines")

# The skill must have been reviewed for the package version in this tree.
version <- sub("^Version:\\s*", "", grep("^Version:", readLines(file.path(root, "DESCRIPTION")), value = TRUE))
if (!identical(field("software-version"), version)) {
  fail("metadata.software-version '", field("software-version"), "' differs from DESCRIPTION '", version,
       "': review the skill for this version")
}

# One canonical source: discovery paths must be symlinks to skills/dndsr.
for (link in c(".claude/skills/dndsr", ".agents/skills/dndsr")) {
  path <- file.path(root, link)
  if (!nzchar(Sys.readlink(path)) || normalizePath(path, mustWork = FALSE) != normalizePath(skill_dir)) {
    fail(link, ": must be a symlink to skills/dndsr")
  }
}

docs <- c(skill_path,
          list.files(file.path(skill_dir, "references"), "\\.md$", full.names = TRUE),
          list.files(file.path(root, "docs", "tutorials"), "\\.md$", full.names = TRUE))
sources <- list.files(file.path(root, "R"), "^cli_.*\\.R$", full.names = TRUE)
commands <- sub("^cli_", "", sub(" <- structure.*", "",
  grep("^cli_[a-z_]+ <- structure", unlist(lapply(sources, readLines, warn = FALSE)), value = TRUE)))
if (!length(commands)) fail("no CLI commands found in R/cli_*.R")

digests <- character()
for (doc in docs) {
  rel <- sub(paste0(root, "/"), "", doc, fixed = TRUE)
  text <- readLines(doc, warn = FALSE)

  # Relative Markdown links must resolve.
  links <- unlist(regmatches(text, gregexpr("\\]\\([^)#[:space:]]+", text)))
  for (target in sub("^\\]\\(", "", links)) {
    if (grepl("^[a-z]+:", target)) next
    if (!file.exists(file.path(dirname(doc), target))) fail(rel, ": broken link ", target)
  }

  # Commands at the start of a line must exist in the CLI.
  used <- sub("^dndsr ([a-zA-Z_-]+).*", "\\1", grep("^dndsr [a-zA-Z_]", text, value = TRUE))
  for (cmd in setdiff(unique(used), commands)) fail(rel, ": unknown CLI command '", cmd, "'")

  # Executable instructions must not use a floating container tag.
  if (any(grepl("^[^|#]*(docker|singularity|apptainer)[^|]*:latest", text))) fail(rel, ": uses a floating :latest tag")
  digests <- c(digests, unlist(regmatches(text, gregexpr("sha256:[0-9a-f]+", text))))
}
digests <- unique(digests)
if (any(nchar(digests) != 71)) fail("malformed image digest: ", paste(digests[nchar(digests) != 71], collapse = ", "))
if (length(digests) > 1) fail("skill and tutorials disagree on the image digest: ", paste(digests, collapse = ", "))

if (length(errors)) {
  writeLines(paste("skill validation FAILED:", errors), stderr())
  quit(status = 1, save = "no")
}
cat(sprintf("skill validation OK: %s, version %s, %d CLI commands, %d documents\n",
            name, version, length(commands), length(docs)))
