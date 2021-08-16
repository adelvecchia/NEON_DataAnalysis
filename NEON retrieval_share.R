#### Intro ####

# This script can be used to download a neon product for every site that has
# data available. The data will be saved in folder (i.e. directory):
# data/raw/product_name/site_name



# These packages are needed to run this script
library(neonUtilities)
library(glue)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(feather)
library(tidyverse)

#install.packages("tidyverse", type="source")

# Here, define the product code for the neon product (dpID) you want to download
# and the name you would like it to be called (such as stream_gases,
# stream_chemistry, etc.). This name is only used to save the files in a folder
# so that it is more easily understood than the neon product code

product_code <- 'DP1.20097.001'
product_name <- 'stream_gasses'

#### Helper Functions ####

# These functions are from MacroSHeds. Understanding how these function works is
# not needed to run the script but take a look if you want!

# get_avail_neon_product_sets will retrieve all the sites and components available
# from the NEON website for the given data product. Neon organizes their data into
# year-month components for most data products, where every year-month has a data
# file (or many) with all corresponding meta data files for that month-year.
get_avail_neon_product_sets <- function(prodcode_full){

    #returns: tibble with url, site_name, component columns

    avail_sets = tibble()

    req = httr::GET(paste0("http://data.neonscience.org/api/v0/products/",
                           prodcode_full))
    txt = httr::content(req, as="text")
    neondata = jsonlite::fromJSON(txt, simplifyDataFrame=TRUE, flatten=TRUE)

    urls = unlist(neondata$data$siteCodes$availableDataUrls)

    avail_sets = stringr::str_match(urls,
                                    '(?:.*)/([A-Z]{4})/([0-9]{4}-[0-9]{2})') %>%
        as_tibble(.name_repair='unique') %>%
        rename(url=`...1`, site_name=`...2`, component=`...3`)

    return(avail_sets)
}


# serialize_list_to_dir will take a list of data tables downloaded from
# NEON and save that list as individual feather files for each dataframe in
# the list. Feather files are much like CSVs and can be loaded with
# read_feather(file_path) or saved with write_feather(dataframe, file_path)
serialize_list_to_dir <- function(l, dest){

    #l must be a named list
    #dest is the path to a directory that will be created if it doesn't exist

    #list element classes currently handled: data.frame, character

    elemclasses = lapply(l, class)

    handled = lapply(elemclasses,
                     function(x) any(c('character', 'data.frame') %in% x))

    if(! all(unlist(handled))){
        stop('Unhandled class encountered')
    }

    dir.create(dest, showWarnings=FALSE, recursive=TRUE)

    for(i in 1:length(l)){

        if('data.frame' %in% elemclasses[[i]]){

            fpath = paste0(dest, '/', names(l)[i], '.feather')
            write_feather(l[[i]], fpath)

        } else if('character' %in% x){

            fpath = paste0(dest, '/', names(l)[i], '.txt')
            readr::write_file(l[[i]], fpath)
        }
    }

}


# read_neon_feather will either read files into a list for one site or will
# read all sites and combine them into a list. The structure of this list mimics
# how a list downloaded from neon would look
#    file_path = where the data is stored. ether data/raw/product_name/site or
#                data/raw/product_name
#    by_site = TRUE or FALSE, TRUE to load one site or FALSE to load
#              all sites. Make sure to include the site you want to load in
#              file_path if by_site = TRUE
read_neon_feathers <- function(file_path, by_site){

    if(by_site == TRUE){

        neon_files <- list.files(file_path, full.names = TRUE)

        file_names <- list.files(file_path)

        neon_list <- lapply(neon_files, read_feather)

        names(neon_list) <- file_names

        return(neon_list)

    } else{

        neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
        file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
        file_names <- str_split_fixed(file_names, '[/]', n = Inf)
        inv_file_names <- file_names[,dim(file_names)[2]]
        site_name <- file_names[,dim(file_names)[2]-1]

        inv_file_names <- unique(inv_file_names)
        site_names <- unique(site_name)

        #categoricalCodes
        #readme
        #validation
        #variables

        data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        final_list <- list()
        for(i in 1:length(data_files)){

            all_files <- tibble()
            for(s in 1:length(site_names)){

                path <- glue('{f}/{s}/{i}.feather',
                             f = file_path,
                             i = data_files[i],
                             s = site_names[s])

                one_site <- read_feather(path)

                all_files <- rbind(all_files, one_site)
            }

            all_site_files <- list(all_files)

            names(all_site_files) <- data_files[i]

            final_list <- c(final_list, all_site_files)
        }

        meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        for(i in 1:length(meta_data_files)){

            path <- glue('{f}/{s}/{i}.feather',
                         f = file_path,
                         i = meta_data_files[i],
                         s = site_names[1])

            meta_file <- read_feather(path)

            meta_list <- list(meta_file)
            names(meta_list) <- meta_data_files[i]

            final_list <- c(final_list, meta_list)
        }
        return(final_list)
    }

}


#### Data Download ####

# This section will dowload all data avalible for the given neon data product

# Get all sites available  for the data product
avail_sets <- get_avail_neon_product_sets(product_code)

# Pull out all the sites available and create a vector that can be looped through
avail_sites <- unique(avail_sets$site_name)

# This for loop will perform the code inside of it for each of the sites that is
# in the vector "avail_sites"
for(j in 1:length(avail_sites)){

    # Defines what site to retrieve
    site_name <- avail_sites[j]

    # Download data product for the site
    data_pile <- neonUtilities::loadByProduct(dpID = product_code,
                                              site = site_name,
                                              package='basic',
                                              check.size=FALSE)

    # Create the file path for where the files will be saved. This can be changed
    # if you  want to save the files in a different location
    raw_data_dest <- glue('{wd}/data/raw/{p}/{s}',
                         wd = getwd(),
                         p = product_name,
                         s = site_name)

    # Save the list of all dataframes for the site as individual feather files
    # in the file path raw_data_dest
    serialize_list_to_dir(data_pile, raw_data_dest)

}


#### Load Data ####

# Here there is an example to load data for one site or for all sites

# Load data for one site

site_name <- 'ARIK'

file_path <- glue('{wd}/data/raw/{p}/{s}',
                     wd = getwd(),
                     p = product_name,
                     s = site_name)

arik_stream_gases <- read_neon_feathers(file_path, by_site = TRUE)

# Load data for all sites

file_path <- glue('{wd}/data/raw/{p}',
                  wd = getwd(),
                  p = product_name)

stream_gases <- read_neon_feathers(file_path, by_site = FALSE)


# Write each dataframe in the list of all sites to separate feather files
# (not necessary, only if you want all sites combined into one file)
file_path_all_data <- glue('{wd}/data/{p}',
                  wd = getwd(),
                  p = product_name)

serialize_list_to_dir(stream_gases, file_path_all_data)
list2env(stream_gases, .GlobalEnv)


