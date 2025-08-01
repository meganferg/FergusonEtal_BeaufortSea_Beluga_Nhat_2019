#Script inla_mesh2sp_sf_fun.R...Megan C. Ferguson...3 July 2024

#MCF adapted this function from Finn Lindgren's function
#posted to https://groups.google.com/g/r-inla-discussion-group/c/z1n1exlZrKM and
#accessed on 3 July 2024

inla.mesh2sp.sf.fun <- function(mesh, crs) {

  #First, create SpatialPolygonsDataFrame for triangles. The SPDF object has the
  #triangles in the same order as in the original mesh, but each triangle loops
  #through the vertices in clockwise order (\code{sp} standard) instead of
  #counterclockwise order (\code{inla.mesh} standard).
  
    tri.SPDF <- SpatialPolygonsDataFrame(
      Sr = SpatialPolygons(lapply(1:nrow(mesh$graph$tv), function(x) {
                                      tv <- mesh$graph$tv[x, , drop = TRUE]
                                      Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                      1:2,
                                      drop = FALSE])),
                                      ID = x)
                                  }),
                            proj4string = crs),
                            data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
                            match.ID = FALSE)
    
  #Convert tri.SPFD to sf object  
    tri.sf <- st_as_sf(tri.SPDF)
    
  #Extract vertices as SpatialPoints object
    vert.SPts <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
    
  #Convert vert.SPts to sf object 
    vert.sf <- st_as_sf(vert.SPts)
    
  #Output all in a list  
    out.list <- list("tri.SPDF"=tri.SPDF, 
                     "tri.sf"=tri.sf, 
                     "vert.SPts"=vert.SPts,
                     "vert.sf"=vert.sf)
    return(out.list)

}
