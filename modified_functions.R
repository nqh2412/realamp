#crop
img_crop <- function(img, x_1, x_2, y_1, y_2) {
  c_im1 <- imsub(img, x > width/x_1, y > height/y_1)
  c_im2 <- imsub(c_im1, x < width/x_2, y < height/y_2)
  return(c_im2)
  UseMethod('img_crop')
}
##Modifying build-in function of 'colocr'

#Roi_select
#' @export
roi_select.cimg <- function(img, threshold = 50, shrink = 5, grow = 5, fill = 5,
                            clean = 5, tolerance = .1, n = 1) {
  
  # check valid input
  if(!missing(threshold) & !is.numeric(threshold)) {
    stop('threshold should be a numeric >= 0 and < 100.')
  }
  if(!missing(threshold) & (threshold >= 100 | threshold < 0)) {
    stop('threshold should be a numeric >= 0 and < 100.')
  }
  
  #crop and rotate the image
  
  #img <- img_crop(img, x_1 = 3.5, x_2 = 1.65, y_1 = 2.3, y_2 = 5)
  img <- imrotate(img,270)
  # change image to gray scale
  img.g <- grayscale(img)
  
  # apply threshold
  img.t <- threshold(img.g, paste0(threshold, '%'))
  
  # change to pixset
  px <- as.pixset(1-img.t)
  
  # apply shrink
  px.m <- shrink(px, shrink)
  
  # apply grow
  px.m <- grow(px.m, grow)
  
  # apply fill
  px.m <- fill(px.m, fill)
  
  # apply clean
  px.m <- clean(px.m, clean)
  
  # add labels when n is provided
  labs.px <- .labels_add(px.m, tolerance = tolerance, n = n)
  attr(img, 'label') <- as.numeric(labs.px)
  
  # return object
  return(img)
}

#Roi_show modified to show only original and the pix set images with hightlight

roi_show.cimg <- function(img, ind = c(1,3)) {
  
  # get labels from img
  # transform labels to cimg object
  labels <- attr(img, 'label')
  dims <- dim(grayscale(img))
  a <- array(labels, dim = dims)
  
  px <- cimg(a)
  
  # merge image
  
  plot(img,
       axes = FALSE,
       main = 'Input image')
  # pixset image
  plot(px,
       axes = FALSE,
       main = 'Read image')
  highlight(px)
  
  # return null
  invisible(NULL)
}

#modified function ROI_CHECK

roi_check.cimg <- function(img, ind = c(1:3)) {
  # get pixel intensities
  label <- attr(img, 'label')
  
  #Red channel intensities of ROIS
  r_channel <- channel(img, ind = 1)
  r_rois <- split(r_channel, label)
  str(r_rois)
  red_int <- data.frame(lapply(r_rois,mean))
  
  #Green channel intensities of ROIS
  g_channel <- channel(img, ind = 2)
  g_rois <- split(g_channel, label)
  str(g_rois)
  gre_int <- data.frame(lapply(g_rois,mean))
  
  #Blue channel intensities of ROIS
  b_channel <- channel(img, ind = 3)
  b_rois <- split(b_channel, label)
  str(b_rois)
  blu_int <- data.frame(lapply(b_rois,mean))
  
  img_value <- rbind.data.frame(gre_int/red_int)
  img_value[,1] <- NULL
  return(img_value)
  # return null
  invisible(NULL)
}


roi_check.list <- function(img, ind = c(1,3)) {
  # get the length of the image list
  img_n <- length(img)
  a <- list()
  # repeat argument to match list length
  if(!is.list(ind)) {
    ind <- rep(list(ind), img_n)
  }
  
  # loop over the images of lists and call roi_check
  for(i in 1:img_n) {
    a[[i]] <-  roi_check(img[[i]],
                         ind = ind[[i]])
  }
  
  #name_list_table <- as.data.frame(bind_rows(name_list, .id = "order"))
  a.table <- as.data.frame(bind_rows(a, .id = "Time_minute"))
  a.table[,2] <-  a.table[,2] + 0.06
  a.table$Time_minute <- order(as.numeric(a.table$Time_minute))
  ml1 <- modlist(a.table, 1, 2:8, model = l4)
  plot(ml1, col = rep(1:7, each = 1))
  
  
  ml2 <- modlist(a.table, 1, 3:8, model = l4)
  
  
  # calculating threshold cycles and plot the threshold cycle values, threshold = "cpD2" or numeric.
  calib(ml2, thresh = 0.95, dil = c(10, 100,1000,10000,100000,1000000),
        group = NULL, plot = TRUE, conf = 0.95, B = 200)
  return(a.table)
  # return null
  invisible(NULL)
}

