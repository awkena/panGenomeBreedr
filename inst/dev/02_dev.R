# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Amend DESCRIPTION with dependencies read from package code parsing
## install.packages('attachment') # if needed.
attachment::att_amend_desc()
attachment::att_from_rscripts() # get list of packages being used.


# import the pangenome package
#usethis::use_dev_package(package = 'panGenomeBreedr',type = 'Imports',  remote = "awkena/panGenomeBreedr")

# # Other packages being used.
# usethis::use_package(package = "grDevices" ,type = 'suggests')
# usethis::use_package(package = "golem" ,type = 'Imports') # necessary
#
#
#
#
#
#
#
# usethis::use_package(package = 'reactable',type = 'suggests')


## Add modules ----
## Create a module infrastructure in R/
golem::add_module('variant_discovery',with_test = TRUE) # variant discovery module
golem::add_module('kasp_marker_design',with_test = TRUE) # Kasp marker design module
golem::add_module('mv_read_kasp_csv',with_test = TRUE) # mv-read_kasp_csv
golem::add_module('mv_nsamples_plate',with_test = TRUE) # mv- nsamples_plate
golem::add_module('mv_get_alleles',with_test = TRUE) # mv- get_alleles
golem::add_module('mv_kasp_color',with_test = TRUE) # mv- kasp_color
golem::add_module('mv_kasp_qc_ggplot',with_test = TRUE) # mv-  kasp_qc_ggplot
golem::add_module('mv_plate_plot',with_test = TRUE) # mv - plate_plot
golem::add_module('mv_pred_sum_stat',with_test = TRUE) #mv - pred_summary, pred_status & plot
golem::add_module('ds_trait_introg_hypothesis_test',with_test = TRUE) # process data for heatmap generation
golem::add_module('ds_marker_ass_bac',with_test = TRUE)  # MABC
golem::add_module('ds_foreground_select',with_test = TRUE)
## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct("helpers", with_test = TRUE) # added already
golem::add_utils("helpers", with_test = TRUE)

## External resources
## Creates .js and .css files at inst/app/www
golem::add_js_file("script")
golem::add_js_handler("handlers")
golem::add_css_file("custom")
golem::add_sass_file("custom")
golem::add_any_file("file.json")

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw(name = "my_dataset", open = FALSE)

## Tests ----
## Add one line by test you want to create
usethis::use_test("app")

# Documentation

## Vignette ----
usethis::use_vignette("PanGB.webapp")
devtools::build_vignettes()

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
##
## (You'll need GitHub there)
usethis::use_github()

# GitHub Actions
usethis::use_github_action()
# Chose one of the three
# See https://usethis.r-lib.org/reference/use_github_action.html
usethis::use_github_action_check_release()
usethis::use_github_action_check_standard()
usethis::use_github_action_check_full()
# Add action for PR
usethis::use_github_action_pr_commands()

# Circle CI
usethis::use_circleci()
usethis::use_circleci_badge()

# Jenkins
usethis::use_jenkins()

# GitLab CI
usethis::use_gitlab_ci()

# Validation check
devtools::check()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")







