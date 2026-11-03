extract_coordinates <- function(x, scheme){
  
  sites <- unique(x, by = "eventID")
  
  o.sf  <- sf::st_as_sf(sites,
                        coords = c("x.coord", "y.coord"),
                        crs = 3035)
  
  o.sf2 <- o.sf
  o.sf2 <- sf::st_transform(o.sf, crs = 4326)
  o.sf2 <- sf::st_coordinates(o.sf2)
  o.sf2 <- as.data.frame(o.sf2)
  
  o.sf2$eventID <- sites$eventID
  names(o.sf2)  <- c("x.coord", "y.coord", "eventID")
  
  return(o.sf2)
}