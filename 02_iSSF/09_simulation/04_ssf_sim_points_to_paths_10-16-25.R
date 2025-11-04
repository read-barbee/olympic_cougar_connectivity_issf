#### SSF sim paths ####

# Author: Read Barbee

# Date:2025-10-06 

# Purpose:

#1. Import simulated points
#2. Convert points to paths
#3. Subset paths intersecting I5
#4. Output gpkg of paths intersecting I5


################################ Functions #################################
# Function to process one batch of path_ids
process_batch <- function(points, batch) {
  # Filter to current subset of points
  sub <- points %>% filter(path_id %in% batch)
  
  # Convert to sf points
  sub_sf <- sub %>%
    st_as_sf(coords = c("x", "y"), crs = 5070, remove = FALSE)
  
  #identify starting and ending sides
  starting_side <- sub_sf %>% 
    filter(step_id == 1) %>% 
    st_join(ocp_split) %>% 
    as.data.frame() %>% 
    rename(starting_side = side_I5) %>% 
    select(path_id, starting_side)
  
  ending_side <- sub_sf %>% 
    filter(step_id == max(sub_sf$step_id, na.rm=T)) %>% 
    st_join(ocp_split) %>% 
    as.data.frame() %>% 
    rename(ending_side = side_I5) %>% 
    select(path_id, ending_side)
  
  # Build lines for this batch
  paths_sf <- sub_sf %>%
    arrange(path_id, step_id) %>%
    group_by(path_id) %>%
    summarise(geometry = st_combine(geometry), .groups = "drop") %>%
    mutate(geometry = st_cast(geometry, "LINESTRING")) %>%
    st_as_sf() %>% 
    left_join(starting_side, by = "path_id") %>% 
    left_join(ending_side, by = "path_id")
  
  
  # Intersection test
  paths_sf <- paths_sf %>%
    mutate(intersects_i5 = lengths(st_intersects(paths_sf, i5, prepared = TRUE)) > 0)
  
  # Return minimal object
  out <- paths_sf %>% select(path_id, intersects_i5, starting_side, ending_side)
  
  return(out)
}


################################ Libraries #################################
library(tidyverse)
library(sf)
library(data.table)
library(flextable)

################################ Import data #################################

batch_size <- 2500  

#filepath <- "/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity"

filepath <- "/Volumes/Reads_Seagate/ms_thesis/ssf_sims"

files <- list.files(filepath, pattern = ".csv", full.names = T)

i5 <- st_read("results/mapping/I-5_v3.gpkg") %>% 
  st_transform(crs = 5070)

ocp_split <- st_read("data/spatial/ocp_study_area_I5_split_5070.gpkg") %>% st_buffer(5000)

################################ Convert points to lines and calculate intersection with I5 #################################

sim_paths_results <- list()

for(i in 1:length(files)){
  
  points <- fread(files[i])
  
  # Extract unique path IDs and create batches
  all_paths <- unique(points$path_id)
  
  batches <- split(all_paths, ceiling(seq_along(all_paths) / batch_size))
  
  
  # Process all batches: about 2.5 minutes for 10,000 paths
  
  results_list <- list()
  
  for(j in 1:length(batches)){
    results_list[[j]] <- process_batch(points, batches[[j]])
    
    print(paste0("Processed batch ", j, " of ", length(batches)))
  }


#combine results
sim_paths_results[[i]] <- bind_rows(results_list)

print(paste0("Processed file ", i, " of ", length(files)))
}

rm(points)
rm(results_list)

sim_paths_all <- bind_rows(sim_paths_results) %>% 
  mutate(path_id = row_number())

#st_write(sim_paths_all, "pop_sim_lines_200000p_4380s_both_sides_I5_10-20-25.gpkg")
sim_paths_all <- st_read("pop_sim_lines_200000p_4380s_both_sides_I5_10-20-25.gpkg")

#subset paths intersecting I5 (for more manageable file size)
sim_paths_i5 <- sim_paths_all %>% 
  slice(1:200000) %>% 
  filter(intersects_i5==T)

#st_write(sim_paths_i5, "pop_sim_lines_intersecting_I5_10-20-25.gpkg", append=F)

sim_paths_i5 <- st_read("pop_sim_lines_intersecting_I5_10-20-25.gpkg")

#calculate point intersections

sim_i5_crossings <- st_intersection(i5, sim_paths_i5)

sim_i5_crossings_f <- sim_i5_crossings %>% 
  st_collection_extract("POINT") %>% 
  select(path_id:geom) %>% 
  mutate(crossing_id = row_number(), .before = path_id)


#st_write(sim_i5_crossings_f, "pop_sim_paths_I5_intersection_points_10-20-25.gpkg")

##Tally directions of paths (side of I5 for start and end points)
directions <- sim_paths_i5 %>% 
  as.data.table() %>% 
  select(-geom) %>% 
  mutate(direction = paste0(starting_side, "_", ending_side)) %>% 
  group_by(direction) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(percentage = round(((n/200000)*100), 2)) %>% 
  add_row(direction = "Total", n = sum(.$n), percentage = sum(.$percentage)) %>% 
  rename(Direction = direction,
         `Num of paths` = n,
         `% of paths` = percentage) %>% 
  flextable() %>% 
  autofit()

#save_as_docx(directions, path = "ssf_sim_i5_crossing_counts_by_direction_10-20-25.docx")


