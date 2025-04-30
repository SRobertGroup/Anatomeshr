

#' @title write .geo file
#'
#'
#' @param data a 2D cell network preprocessed with prep_geo()
#' @param dim 
#' @param path_geo file name to store the output 
#' @param Celldomain should the cell domain be meshed? True/False (Default= F)
#' @keywords GMSH Geo
#' @export
#'


write_geo <- function(data, dim = 2, path_geo, Celldomain =F){
  `%!in%` <- compose(`!`, `%in%`)
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