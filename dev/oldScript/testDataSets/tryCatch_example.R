tryCatch(
  expr = {
    run_check(counts_inpu, group_list, tempdir())
  },
  error = function(e){
    usethis::ui_oops("check_failed")
  },
  warning = function(w){
    usethis::ui_warn("pleasecheckmore")
  },
  finally = {
    usethis::ui_line("done now")
  }
)
