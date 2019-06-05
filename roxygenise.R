############# build all documentation for the package and update namespace
file.remove(list.files('man/', pattern = '.Rd', full.names = TRUE))
roxygen2::roxygenise()
############# load the pks to access functions
devtools::load_all()
############# run ./tests tests which then redirects to testthat.R. testhat/helper-*.R files are run before tests
devtools::test()
############# test examples
devtools::run_examples()
## path to Rd files - will let you know if there are any errors/warnings
# testthat::test_example("man/parent_base_ext.Rd")
############# manual tests
## test files that cannot be included in standard testthat pipeline
# invisible(lapply(list.files("tests/manual", full.names = TRUE), source))
############# check
# devtools::check()
############# build it
## by default, it is made in parent directoy of the package, we manually set it.
## testhat has changed wd
# setwd('/Users/alabadi/Documents/_Projects/_Personal/someWrappers')

  # binary_dir <- "../__binary"
  # if(!dir.exists(binary_dir)){
  #   dir.create(binary_dir)
  # }
  # pkg <- devtools::build(path = binary_dir)
  # ############# install it in my library
  # install.packages(pkg, repos = NULL)

# CMD + SHift + B
