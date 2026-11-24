#' Reproducible scripts from the Shiny app
#'
#' Each Shiny module describes the panGB calls it ran as a "step": the R code
#' with the user's settings filled in, the settings themselves, the input
#' files it read, and the steps it needs first. A "Reproduce this" button
#' (`mod_repro_ui()` / `mod_repro_server()`) joins a step and everything
#' upstream of it into one script, and offers it as text, as an `.R` file,
#' or as a zip bundle that reruns outside the app with
#' `source("<tab name>.R")`, e.g. `source("marker_validation.R")`.
#'
#' The helpers below are internal; see `repro_recorder()` for how modules
#' register their steps.
#'
#' @name repro
#' @keywords internal
#' @noRd
NULL


# CODE TEMPLATES ------------------------------------------------------------

#' Mark a value as an R object name in generated code
#'
#' `repro_call(x = repro_sym("kasp_data"))` writes `x = kasp_data` rather
#' than the quoted string `x = "kasp_data"`.
#'
#' @param x Object name, or any R expression, as a string.
#' @return The string, with class `repro_sym`.
#' @noRd
repro_sym <- function(x) {
  structure(x, class = "repro_sym")
}


#' Write a value as R code
#'
#' @param x A value: atomic vector (named or not), NULL, a list of those,
#'   or a `repro_sym()`.
#' @return One string of R code that recreates `x`.
#' @noRd
repro_value <- function(x) {
  if (inherits(x, "repro_sym")) {
    return(unclass(x))
  }
  if (is.null(x)) {
    return("NULL")
  }
  if (is.list(x)) {
    items <- vapply(x, repro_value, character(1))
    nms <- names(x)
    if (!is.null(nms)) {
      items <- ifelse(nzchar(nms), paste(repro_name(nms), "=", items), items)
    }
    return(paste0("list(", paste(items, collapse = ", "), ")"))
  }
  # Plain values: drop attributes other than names so the code stays short,
  # and write whole numbers without an L suffix or e+07 notation
  x <- stats::setNames(as.vector(x), names(x))
  if (is.integer(x)) x <- stats::setNames(as.numeric(x), names(x))
  old <- options(scipen = 100)
  on.exit(options(old), add = TRUE)
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}


#' Quote argument names that are not syntactic
#' @noRd
repro_name <- function(x) {
  ok <- make.names(x) == x
  ifelse(ok, x, paste0("`", x, "`"))
}


#' Write a function call as R code, one argument per line
#'
#' @param fn Function name, e.g. `"kasp_color"`.
#' @param args Named list of argument values (see `repro_value()`).
#'   Unnamed elements are passed by position.
#' @param assign Optional object name to assign the result to.
#' @return A character vector of code lines.
#' @examples
#' \dontrun{
#' repro_call("kasp_color",
#'            list(x = repro_sym("kasp_data"), sep = ":"),
#'            assign = "kasp_coded")
#' #> kasp_coded <- kasp_color(
#' #>   x = kasp_data,
#' #>   sep = ":"
#' #> )
#' }
#' @noRd
repro_call <- function(fn, args = list(), assign = NULL) {
  lhs <- if (is.null(assign)) "" else paste0(assign, " <- ")
  if (length(args) == 0) {
    return(paste0(lhs, fn, "()"))
  }

  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  vals <- vapply(args, repro_value, character(1))
  arg_lines <- ifelse(nzchar(nms), paste(repro_name(nms), "=", vals), vals)
  arg_lines <- paste0("  ", arg_lines, c(rep(",", length(arg_lines) - 1), ""))

  c(paste0(lhs, fn, "("), arg_lines, ")")
}


#' One step of a reproducible script
#'
#' @param title Short heading shown as a comment above the code.
#' @param code Character vector of code lines (see `repro_call()`).
#' @param settings Named list of the settings the user chose, written to
#'   `settings.yml`.
#' @param files Named character vector of input files: names are the file
#'   names used in the script (under `data/`), values are local paths to
#'   copy into the bundle.
#' @param needs Keys of the steps that must run first.
#' @return A list with class `repro_step`.
#' @noRd
repro_step <- function(title, code, settings = list(), files = character(),
                       needs = character()) {
  structure(
    list(
      title = title,
      code = as.character(code),
      settings = settings,
      files = files,
      needs = needs
    ),
    class = "repro_step"
  )
}


#' Path of a bundled input file, as written in the script
#'
#' @param name File name, as given in `repro_step(files = )`.
#' @return The relative path `data/<name>`.
#' @noRd
repro_data_path <- function(name) {
  file.path("data", repro_file_name(name))
}


#' Make a file name safe to write into a bundle
#' @noRd
repro_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]", "_", basename(x))
  x[!nzchar(x) | x %in% c(".", "..")] <- "input"
  x
}


# REGISTRY ------------------------------------------------------------------

#' Shared registry of reproducible steps
#'
#' Create one in the app server and pass it to the modules of a workflow.
#' Each module registers a reactive that returns its `repro_step()` (or
#' NULL when it has not run yet). Steps are looked up by key, so a module's
#' "Reproduce this" button can pull in the steps upstream of it.
#'
#' @return A list of functions: `register(key, step)`, `get(key)`,
#'   `keys()` and `collect(key)`.
#' @noRd
repro_recorder <- function() {
  steps <- new.env(parent = emptyenv())

  register <- function(key, step) {
    stopifnot(is.character(key), length(key) == 1, is.function(step))
    assign(key, step, envir = steps)
    invisible(key)
  }

  get_step <- function(key) {
    if (!exists(key, envir = steps, inherits = FALSE)) {
      return(NULL)
    }
    # NULL while the module has not run, or if its inputs are not ready
    if_ready(get(key, envir = steps, inherits = FALSE)())
  }

  # The step and everything upstream of it, in the order they must run.
  # Returns list(steps = <named list>, missing = <keys not ready>)
  collect <- function(key) {
    out <- list()
    missing <- character()
    visit <- function(k, seen) {
      if (k %in% seen) stop("Circular step dependency at '", k, "'.")
      if (k %in% names(out) || k %in% missing) return(invisible())
      st <- get_step(k)
      if (is.null(st)) {
        missing <<- c(missing, k)
        return(invisible())
      }
      for (n in st$needs) visit(n, c(seen, k))
      out[[k]] <<- st
    }
    visit(key, character())
    list(steps = out, missing = missing)
  }

  list(
    register = register,
    get = get_step,
    keys = function() ls(steps),
    collect = collect
  )
}


# SCRIPT AND BUNDLE ---------------------------------------------------------

#' Version string of the installed panGB
#' @noRd
repro_pkg_version <- function() {
  tryCatch(as.character(utils::packageVersion("panGenomeBreedr")),
           error = function(e) "unknown")
}


#' Join steps into one R script
#'
#' @param steps Named list of `repro_step()` objects, in run order.
#' @param title Name of the analysis, used in the header.
#' @param time Time stamp written in the header.
#' @param script_name File name of the script, named after the app tab,
#'   e.g. `"marker_validation.R"`.
#' @param notes Header lines on where the script reads its inputs. The
#'   default is for a bundle with a `data/` folder.
#' @return One string with the whole script.
#' @noRd
repro_script <- function(steps, title = "panGenomeBreedr analysis",
                         time = Sys.time(), script_name = "analysis.R",
                         notes = c(
                           "Inputs are read from data/ and results are written to results/.",
                           "The settings you chose in the app are also listed in settings.yml."
                         )) {
  rule <- strrep("-", 72)
  header <- c(
    paste("#", rule),
    paste("#", title),
    "#",
    paste0("# Exported from the panGenomeBreedr Shiny app on ",
           format(time, "%Y-%m-%d %H:%M %Z"), "."),
    paste0("# panGenomeBreedr ", repro_pkg_version(), ", ",
           R.version$version.string, "."),
    "#",
    "# To rerun: start R in the folder that holds this script and run",
    paste0("#   source(\"", script_name, "\")"),
    paste("#", notes),
    paste("#", rule),
    "",
    "library(panGenomeBreedr)",
    "",
    "dir.create(\"results\", showWarnings = FALSE)"
  )

  body <- unlist(lapply(seq_along(steps), function(i) {
    st <- steps[[i]]
    c("", paste0("# ", i, ". ", st$title, " ",
                 strrep("-", max(3, 68 - nchar(st$title) - nchar(i)))),
      st$code)
  }), use.names = FALSE)

  paste(c(header, body, "", "message(\"Done. Results are in the results/ folder.\")", ""),
        collapse = "\n")
}


#' Settings of all steps, as one list for settings.yml
#' @noRd
repro_settings <- function(steps, title, time = Sys.time()) {
  list(
    analysis = title,
    exported = format(time, "%Y-%m-%d %H:%M:%S %Z"),
    panGenomeBreedr = repro_pkg_version(),
    R = R.version$version.string,
    steps = lapply(steps, function(st) {
      c(list(title = st$title), st$settings,
        if (length(st$files)) list(files = as.list(repro_data_path(names(st$files)))))
    })
  )
}


#' Write a reproducibility bundle
#'
#' The zip holds the script (`script_name`), `settings.yml` (or `settings.json` if the
#' yaml package is not installed), the input files under `data/` and
#' `sessionInfo.txt`.
#'
#' @param file Path of the zip to write.
#' @param steps Named list of `repro_step()` objects, in run order.
#' @param title Name of the analysis.
#' @param script_name File name of the script, named after the app tab.
#' @param max_copy_mb Input files larger than this are not copied; the
#'   README asks the user to put them in `data/` by hand.
#' @return `file`, invisibly.
#' @noRd
repro_bundle <- function(file, steps, title = "panGenomeBreedr analysis",
                         script_name = "analysis.R", max_copy_mb = 200) {
  time <- Sys.time()
  dir <- tempfile("pangb_bundle_")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  dir.create(file.path(dir, "data"), recursive = TRUE)

  writeLines(repro_script(steps, title, time, script_name),
             file.path(dir, script_name))

  settings <- repro_settings(steps, title, time)
  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(settings, file.path(dir, "settings.yml"))
  } else {
    jsonlite::write_json(settings, file.path(dir, "settings.json"),
                         auto_unbox = TRUE, pretty = TRUE, null = "null")
  }

  # Input files
  files <- unlist(unname(lapply(steps, `[[`, "files")))
  files <- files[!duplicated(names(files))]
  not_copied <- character()
  for (nm in names(files)) {
    src <- files[[nm]]
    dest <- file.path(dir, repro_data_path(nm))
    size_mb <- file.size(src) / 1024^2
    if (!is.na(size_mb) && size_mb <= max_copy_mb) {
      file.copy(src, dest)
    } else {
      not_copied <- c(not_copied, repro_file_name(nm))
    }
  }

  readme <- c(
    title,
    "",
    "Rerun this analysis without the app:",
    "  1. Unzip this folder.",
    "  2. Start R in the unzipped folder (or setwd() to it).",
    paste0("  3. Run: source(\"", script_name, "\")"),
    "",
    "Tables and plots are written to results/.",
    "settings.yml lists every setting chosen in the app.",
    "sessionInfo.txt lists the R and package versions used."
  )
  if (length(not_copied)) {
    readme <- c(readme, "",
                "These input files were too large to include. Copy them into data/ before running:",
                paste("  -", not_copied))
  }
  writeLines(readme, file.path(dir, "README.txt"))

  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(dir, "sessionInfo.txt"))

  entries <- list.files(dir, all.files = FALSE, recursive = FALSE)
  if (file.exists(file)) unlink(file)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zipfile = file, files = file.path(dir, entries))
  } else {
    old <- setwd(dir)
    on.exit(setwd(old), add = TRUE, after = FALSE)
    utils::zip(zipfile = file, files = entries, flags = "-r9Xq")
  }
  invisible(file)
}


# BUTTON MODULE -------------------------------------------------------------

#' "Reproduce this" button
#'
#' @param id Module ID.
#' @param label Button text.
#' @param class CSS classes. The default fills the width of a sidebar; use
#'   `"btn-outline-secondary btn-sm"` in a card header.
#' @return A Shiny tag.
#' @noRd
mod_repro_ui <- function(id, label = "Reproduce this",
                         class = "btn-outline-secondary btn-sm w-100 mt-3") {
  ns <- NS(id)
  actionButton(ns("open"), label, icon = icon("code"), class = class)
}


#' Server side of the "Reproduce this" button
#'
#' Opens a dialog with the script for step `key` and the steps upstream of
#' it, with buttons to copy it, download it as `script_name`, or download the
#' zip bundle.
#'
#' @param id Module ID, matching `mod_repro_ui()`.
#' @param recorder A `repro_recorder()`.
#' @param key Key of the step to reproduce.
#' @param title Name of the analysis, used in the header.
#' @param script_name File name of the script, named after the app tab,
#'   e.g. `"marker_validation.R"`. The zip is named after it too.
#' @param labels Named character vector: plain names of the steps, used to
#'   say which step still has to run.
#' @param optional Keys of further steps added to the script when they have
#'   run (with the steps they need), e.g. a plot the user may not have
#'   drawn.
#' @param bundle FALSE to offer the script only, with no zip bundle, e.g.
#'   when the script reads from a database rather than from files.
#' @param notes Header lines of the script on where it reads its inputs
#'   (see `repro_script()`). NULL keeps the default.
#' @noRd
mod_repro_server <- function(id, recorder, key,
                             title = "panGenomeBreedr analysis",
                             script_name = "analysis.R",
                             labels = character(),
                             optional = character(),
                             bundle = TRUE,
                             notes = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    make_script <- function(steps) {
      if (is.null(notes)) {
        repro_script(steps, title, script_name = script_name)
      } else {
        repro_script(steps, title, script_name = script_name, notes = notes)
      }
    }

    collected <- reactive({
      res <- recorder$collect(key)
      for (k in optional) {
        extra <- recorder$collect(k)
        if (length(extra$missing) == 0) {
          new <- setdiff(names(extra$steps), names(res$steps))
          res$steps[new] <- extra$steps[new]
        }
      }
      res
    })

    # Zip name: the tab name plus a time stamp, e.g.
    # marker_validation_20261122_1630.zip
    zip_name <- function() {
      paste0(sub("\\.[Rr]$", "", script_name), "_",
             format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
    }

    observeEvent(input$open, {
      res <- collected()

      if (length(res$missing)) {
        pending <- ifelse(res$missing %in% names(labels),
                          labels[res$missing], res$missing)
        showModal(modalDialog(
          title = "Nothing to reproduce yet",
          easyClose = TRUE,
          footer = modalButton("Close"),
          tags$p("Run these steps first, then try again:"),
          tags$ul(lapply(unname(pending), tags$li))
        ))
        return()
      }

      script <- make_script(res$steps)
      code_id <- ns("code")

      intro <- if (bundle) {
        tags$p(
          "This is the panGenomeBreedr code behind these results, with your settings filled in.",
          "Download the bundle to rerun it later without the app: unzip it, start R in that folder and run",
          tags$code(paste0("source(\"", script_name, "\")")), "."
        )
      } else {
        tags$p(
          "This is the panGenomeBreedr code behind these results, with your settings filled in.",
          "Download the script to rerun it later without the app: start R in the folder you saved it to and run",
          tags$code(paste0("source(\"", script_name, "\")")), ".",
          "It reads from the same database, so it needs access to it."
        )
      }

      showModal(modalDialog(
        title = tagList(icon("code"), " Reproduce this analysis"),
        size = "l",
        easyClose = TRUE,
        intro,
        tags$pre(
          id = code_id,
          style = "max-height: 55vh; overflow: auto; font-size: 0.8rem; white-space: pre;",
          script
        ),
        footer = tagList(
          tags$button(
            type = "button",
            class = "btn btn-outline-secondary",
            onclick = sprintf(repro_copy_js, code_id),
            icon("copy"), " Copy code"
          ),
          downloadButton(ns("download_script"), "Script (.R)",
                         class = if (bundle) "btn-outline-success" else "btn-success"),
          if (bundle) {
            downloadButton(ns("download_bundle"), "Bundle (.zip)",
                           class = "btn-success")
          },
          modalButton("Close")
        )
      ))
    })

    output$download_script <- downloadHandler(
      filename = function() script_name,
      content = function(file) {
        res <- collected()
        req(length(res$missing) == 0)
        writeLines(make_script(res$steps), file)
      }
    )

    output$download_bundle <- downloadHandler(
      filename = zip_name,
      content = function(file) {
        res <- collected()
        req(length(res$missing) == 0)
        repro_bundle(file, res$steps, title, script_name = script_name)
      },
      contentType = "application/zip"
    )
  })
}


#' JavaScript for the "Copy code" button; `%s` is the code block's id
#' @noRd
repro_copy_js <- paste0(
  "(function(btn){",
  "var t=document.getElementById('%s').innerText;",
  "var done=function(){var h=btn.innerHTML;btn.innerHTML='Copied';",
  "setTimeout(function(){btn.innerHTML=h;},1500);};",
  "if(navigator.clipboard&&window.isSecureContext){",
  "navigator.clipboard.writeText(t).then(done);}",
  "else{var a=document.createElement('textarea');a.value=t;",
  "document.body.appendChild(a);a.select();document.execCommand('copy');",
  "document.body.removeChild(a);done();}",
  "})(this)"
)
