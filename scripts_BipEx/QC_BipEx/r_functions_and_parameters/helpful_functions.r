ellipse <- function(angle, centre_x, centre_y, height, width) {
  k <- 1
  x <- matrix(nrow = 100, ncol=2)
  for(i in seq(0, 2*pi, length.out=100)) {
    x[k,] <- c(centre_x, centre_y) + 
      width * cos(i) * c(cos(angle), sin(angle)) + 
      height * sin(i) * c(-sin(angle), cos(angle))
    k <- k+1
  }
  return(x)
}

in_ellipse <- function(x, y, angle, centre_x, centre_y, height, width) {
	return((((x-centre_x) * cos(angle) + (y-centre_y) * sin(angle)) / width)^2 +
		   (((x-centre_x) * sin(angle) - (y-centre_y) * cos(angle)) / height)^2 < 1)
}

above_below_line <- function(x_points, y_points, m, c, above_or_below) {
	if (above_or_below == "above") {
		return(y_points > (m * x_points + c))
	} else if (above_or_below == "below"){
		return(y_points < (m * x_points + c))
	} else {
		return("Error!")
	}
}