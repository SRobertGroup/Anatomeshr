
#' @title prepare .geo file
#'
#'
#' @param section a 2D cell network
#' @param cell_wall_thickness the table or unique value of inner wall for the cell network
#' @param corner_smoothing the table or unique value of smoothing algorithm (smoothr with ksmooth)
#' @keywords GMSH Geo
#' @export
#'
#'



require(sf)
require(smoothr)
require(tidyverse)

prep_geo <- function(section, cell_wall_thickness = 0.2, corner_smoothing=0.5){
  
  if(length(cell_wall_thickness)==1){
    cell_wall_thickness = tibble(type = "default", value = cell_wall_thickness)
  }
  if(length(cell_wall_thickness$value[cell_wall_thickness$type == "outerwall"])==0){
    cell_wall_thickness = rbind(cell_wall_thickness,
                                tibble(type = "outerwall",
                                       value = cell_wall_thickness$value[cell_wall_thickness$type == "default"]))
  }
  
  if(length(corner_smoothing)==1){
    corner_smoothing = tibble(type = "default", value = corner_smoothing)
  }
  
  id_cell_vector = unique(section$id_cell)
  root_cell =NULL
  polygons_list = list()
  k = 0
  for(i in id_cell_vector){
    tmp = section%>% filter(id_cell == i)
    if(nrow(tmp)>2){
      k = k+1
      wall_type = tmp$type[1]
      if(wall_type %!in% cell_wall_thickness$type){wall_type = "default"}
      
      sf_linestring <- sf::st_sfc(st_linestring(as.matrix(rbind(tmp[, c("x", "y")],tmp[1, c("x", "y")]))), crs = 2056)
      my_multilinestring = sf::st_sf(geom = sf::st_sfc(sf_linestring), crs = 2056)

      # Union and polygonize the multilinestring
      poly_geom <- sf::st_union(my_multilinestring) %>% sf::st_polygonize()
      
      # Check if the resulting geometry collection is empty
      if (!sf::st_is_empty(poly_geom)) {
        r_poly <- poly_geom %>% sf::st_collection_extract()
      } else {break}
      
      r_poly_smooth <- smoothr::smooth(r_poly, method = "ksmooth",
                                       smoothness = corner_smoothing$value[corner_smoothing$type == wall_type][1])
      shrunken_polygon <- st_buffer(r_poly_smooth,
                                    -cell_wall_thickness$value[cell_wall_thickness$type == wall_type])
      shrunken_polygon  <- st_simplify(shrunken_polygon, dTolerance = 0.1, preserveTopology = TRUE)
      
      swollen_polygon <- st_buffer(r_poly, cell_wall_thickness$value[cell_wall_thickness$type == "outerwall"])
      # Get the coordinates
      coords <- sf::st_coordinates(shrunken_polygon )[, 1:2]
      
      pol = tibble(x = coords[,1], y = coords[,2],
                   id_cell = i)
      pol = pol %>%
        mutate(x1 = x,
               y1 = y,
               x2 = c(pol$x[-1], pol$x[1]),
               y2 = c(pol$y[-1], pol$y[1]))
      
      root_cell = rbind(root_cell, pol)
      
      polygons_list[[k]] = swollen_polygon
      
    }
  }
  
  # Convert list to sf object
  polygons_sf <- do.call(rbind, lapply(seq_along(polygons_list), function(j) {
    st_sf(id = j, geometry = st_geometry(polygons_list[[j]][1]))
  }))
  
  polygons_sfc <- st_as_sfc(polygons_sf)
  # Perform the union operation
  final_polygon_raw <- st_union(polygons_sfc)
  
  # if(length(corner_smoothing$value[corner_smoothing$type == "outerwall"])==0){
  #   corner_smoothing = rbind(corner_smoothing, tibble(type = "outerwall",
  #                                                     value = corner_smoothing$value[corner_smoothing$type == "default"]))
  # }
  # check = F
  # while(check){
  #   final_polygon <- smoothr::smooth(final_polygon_raw, method = "ksmooth",
  #                                    smoothness = corner_smoothing$value[corner_smoothing$type == "outerwall"])
  #   check = any(is.na((sf::st_coordinates(final_polygon)[,c(1,2)])))
  #   corner_smoothing$value[corner_smoothing$type == "outerwall"] = corner_smoothing$value[corner_smoothing$type == "outerwall"]*1.1
  # }
  # final_polygon <- sf::st_make_valid(final_polygon)
  simplified_polygon <- st_simplify(final_polygon_raw, dTolerance = 0.1, preserveTopology = TRUE)
  # Get the coordinates
  wall_coords <- sf::st_coordinates(simplified_polygon)[, 1:2]
  
  wall = tibble(x = wall_coords[,1], y = wall_coords[,2], id_cell = k+1)
  wall = wall %>%
    mutate(x1 = x,
           y1 = y,
           x2 = c(wall$x[-1], wall$x[1]),
           y2 = c(wall$y[-1], wall$y[1]))
  
  return(rbind (wall, root_cell)%>%mutate(res = 1))
}

`%!in%` <- compose(`!`, `%in%`)

require(RImageJROI)

#Read the cell.roi and convert it to a polygon
get_cell <- function(roi_path){
  x <- RImageJROI::read.ijroi(file = roi_path)
  i = parse_number(x$name)
  poly_sf <- sf::st_polygon(list(rbind(x$coords,x$coords[1,])))
  poly_cell = sf::st_sf(geom = sf::st_sfc(poly_sf), crs = 2056)
  coords <- sf::st_coordinates(poly_cell)[, 1:2]
  pol = tibble(x = coords[,1], y = coords[,2],
               id_cell = i)
  pol = pol %>%
    mutate(x1 = x,
           y1 = y,
           x2 = c(pol$x[-1], pol$x[1]),
           y2 = c(pol$y[-1], pol$y[1]))
  return(pol)
}




write_geo <- function(data, dim = 2, path_geo, Celldomain =F){
  
  date = Sys.time()
  x1 = paste0('// Gmsh project created on ', date,'\nSetFactory("OpenCASCADE");\n//+\n')
  
  if(dim == 2){
    data$z = 0
    data$z1 = 0
    data$z2 = 0
  }
  
  center = c(mean(data$x), mean(data$y))
  central = data %>% filter(id_cell != max(id_cell)) %>% 
    dplyr::group_by(id_cell) %>% 
    dplyr::summarise(mx = mean(x), my = mean(y), .groups = "drop")%>%
    mutate(euc = sqrt((mx-center[1])^2+(my-center[2])^2))
  center_cell_id = central$id_cell[central$euc == min(central$euc)][1] 
  
  k = h = j = 1
  txt = x1
  for(i in sort(unique(data$id_cell))){
    tmp = data%>%filter(id_cell == i)
    i_point = unique_point_id(tmp)
    
    i_point$idx = seq(k,-1+k+nrow(i_point),1)
    i_point$idx2 = c(i_point$idx[-1],i_point$idx[1])
    
    if(i < max(data$id_cell)){
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',1.0};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nSurface(',
                    h,") = {", j,'};\n//+\n', collapse = "")
      if(Celldomain){
        Physical = paste0("Physical Surface(",h,") = {",h,"};\n")
      }else{
        Physical = paste0("//Physical Surface(",h,") = {",h,"};\n")
      }
      
      txt = paste0(txt, dots, Lines, Surf, Physical)
      if(i == center_cell_id){
        fix_curve = paste0('Physical Curve("fix", 2) = {',paste0(i_point$idx, collapse = ", "),'};\n')
        fix_id= i_point$idx
      }
    }else{
      
      all_inner_id = 1:(k-1)
      all_inner_id = all_inner_id[all_inner_id %!in% fix_id]
      Physical_curve = paste0('Physical Curve("inner", 1) = {',paste0( all_inner_id  , collapse = ", "),'};\n')
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',0.5};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nPlane Surface(',
                    h,") = {", paste0(sort(seq(1,j,2), decreasing = T), collapse = ", "),'};\n//+\n', collapse = "")
      Physical = paste0("Physical Surface(0) = {",h,"};\n")
      outer_curve = paste0('Physical Curve("outer", 3) = {',paste0(i_point$idx, collapse = ", "),'};\n')
      txt = paste0(txt, dots, Lines, Surf, Physical, Physical_curve, fix_curve, outer_curve)
    }
    
    h = h+1
    j = j+2
    k = max(i_point$idx)+1
  }
  write(txt, path_geo)
}


#' @title remove duplicated points
#'
#'
#' @param tmp data to filter
#' @keywords GMSH Geo
#' @export
#'
#'

unique_point_id <- function(tmp){
  tmp$id_point = paste0('x:',tmp$x1,'_y:',tmp$y1,'_z:',tmp$z1)
  tmp$id_point2 = paste0('x:',tmp$x2,'_y:',tmp$y2,'_z:',tmp$z2)
  i_point = tibble(id_point = unique(c(tmp$id_point, tmp$id_point2)))
  
  x_coor = unlist(str_split(unlist(str_split(i_point$id_point, pattern = 'x:')),pattern = "_y"))
  y_coor = unlist(str_split(unlist(str_split(i_point$id_point, pattern = 'y:')),pattern = "_z"))
  z_coor = unlist(str_split(i_point$id_point, pattern = 'z:'))
  i_point$x = parse_number(x_coor[seq(2,length(x_coor),3)])
  i_point$y = parse_number(y_coor[seq(2,length(y_coor),3)])
  i_point$z = parse_number(z_coor[seq(2,length(z_coor),2)])
  
  i_point
  return(i_point)
}

