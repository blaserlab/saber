# package installation ----------------------------------------------------

# # By default, the newest version of all packages available  in the cache
# # at the time of project initiation are linked to the projet library

# # Use this to update the entire project library
# # to the newest versions available in the cache
# blaseRtemplates::get_new_library(newest_or_file = "newest")

# # Use this to update the entire project library
# # to another version of the project library
# blaseRtemplates::get_new_library(newest_or_file = "<path/to/file>")

# # Use this to get or update a package from the cache
# blaseRtemplates::install_one_package("<package name>", how = "link_from_cache")

# # If you need a new package or an update from a repository, try this:
# blaseRtemplates::install_one_package("<package name>", how = "new_or_update")

# # use "bioc::<package name>" for bioconductor packages
# # use "<repo/package name>" for github source packages

# # If you are installing from a "tarball", use this:
# blaseRtemplates::install_one_package("/path/to/tarball.tar.gz")


# load and attach packages ------------------------------------------------

library("blaseRtemplates")
library("conflicted")
library("tidyverse")
library("gert")

# load, install, and/or update the project data -----------------------------

# blaseRtemplates::project_data("<project.datapkg>")

