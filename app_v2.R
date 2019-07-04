library(shiny)
library(shinyBS)
library(imager)
library(colocr)
library(dplyr)
library(qpcR)

#pcrfit from "qpcR"
plot.pcrfit <- function(
  x, 
  type = c("all", "single", "3D", "image"),
  fitted = TRUE, 
  add = FALSE,
  col = NULL, 
  par2D = list(),
  par3D = list(),
  ...) 
{
  type <- match.arg(type)    
  object <- x
  
  print(class(x))
  
  if (class(x)[1] != "modlist") modLIST <- list(object) else modLIST <- object      
  
  ## extract cycles and fluorescence values from all curves
  allCYC <- lapply(modLIST, function(x) x$DATA[, 1])
  allFLUO <- lapply(modLIST, function(x) x$DATA[, 2])
  vecCYC <- do.call(c, allCYC)
  vecFLUO <- do.call(c, allFLUO)
  
  ## make unique cycles  
  CYC <- unique(as.numeric(vecCYC))  
  CYC <- CYC[!is.na(CYC)]    
  
  ## calculate min and max fluo values for defining ylim 
  MIN <- min(vecFLUO, na.rm = TRUE)   
  MAX <- max(vecFLUO, na.rm = TRUE)     
  
  ## length of 'modlist'
  LEN <- length(modLIST)
  ## names of 'modlist'
  NAMES <- sapply(modLIST, function(x) x$names)   
  
  ## define plotting colors
  if (is.null(col)) {
    COL <- rep(1, LEN)    
    if (class(object)[2] == "replist") COL <- rainbow(attr(object, "nlevels"))     
  } else COL <- rep(col, length.out = LEN)   
  
  ## 3D plot empty setup using par3D parameters
  if (type == "3D") {
    do.call(plot3d, modifyList(list(x = CYC, y = 1:LEN, z = MAX, type = "n", axes = FALSE, box = FALSE, xlab = "", 
                                    ylab = "", zlab = "", zlim = c(0, 1.1 * MAX)), par3D))
    do.call(axis3d, modifyList(list('x', at = pretty(CYC), cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Time (minutes)", 'x', line = 2), par3D))     
    do.call(axis3d, modifyList(list('y', at = 1:LEN, label = NAMES, cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Run", 'y', line = 2), par3D))
    do.call(axis3d, modifyList(list('z', cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Green/Red ratio", 'z', line = 2), par3D))
  }   
  
  ## standard 'all' plot empty setup
  if (type == "all" && !add) {   
    tempLIST <- modifyList(list(CYC, rep(MAX, length(CYC)), ylim = c(MIN, MAX), 
                                xlab = "Time (minutes)", ylab = "Green/Red Ratio", las = 1), par2D)
    tempLIST$type <- "n"
    do.call(plot, tempLIST)   
  }
  
  ## plot matrix empty setup
  if (type == "single") {
    DIM <- ceiling(sqrt(LEN))
  } 
  
  ## image plot 
  if (type == "image") {
    RUNS <- 1:length(modLIST)
    nRUNS <- length(RUNS)
    ## unique cycles
    CYCS <- unique(unlist(lapply(modLIST, function(x) x$DATA[, 1])))
    nCYCS <- length(CYCS)
    ## convert list with fluo data to matrix, fll unequal length with NA
    allLIST <- lapply(modLIST, function(x) x$DATA[, 2])
    maxCYCS <- max(sapply(allLIST, length))
    for (i in 1:length(allLIST)) allLIST[[i]] <- c(allLIST[[i]], rep(NA, maxCYCS - length(allLIST[[i]])))
    allDAT <- do.call(cbind, allLIST)
    ## image setup
    allDAT <- allDAT[, ncol(allDAT):1]
    image(allDAT, col = heat.colors(100), axes = FALSE, xlab = "Time (minutes)", ylab = "Runs")
    axis(1, at = seq(0, 1, length.out = nCYCS), labels = CYCS)
    axis(2, at = seq(0, 1, length.out = nRUNS), labels = rev(RUNS))
  }
  
  ## iterate through all curves
  for (i in 1:LEN) {
    DATA <- modLIST[[i]]$DATA    
    DATA <- na.omit(DATA)      
    FITTED <- fitted(modLIST[[i]])       
    m <- match(CYC, DATA[, 1])
    m <- na.omit(m)
    
    ## plot 3D curves
    if (type == "3D") {
      do.call(points3d, modifyList(list(x = DATA[, 1], y = i, z = DATA[, 2], color = COL[i]), par3D))
      if (!is.null(FITTED) && fitted) do.call(lines3d, modifyList(list(x = DATA[m, 1], y = i, z = FITTED[m], color = COL[i]), par3D))      
    }
    
    ## plot 2D curves
    if (type == "all") {
      do.call(points, modifyList(list(DATA[, 1], DATA[, 2], col = COL[i]), par2D))
      if (!is.null(FITTED) && fitted) do.call(lines, modifyList(list(DATA[m, 1], FITTED[m], col = COL[i]), par2D)) 
    } 
    
    ## plot matrix curves
    if (type == "single") {
      NAME <- NAMES[i]
      ## color by failed fit or failed structure
      if (grepl("\\*\\*[[:alnum:]]*", NAME)) colMAIN <- "blue" 
      else if (grepl("\\*[[:alnum:]]*", NAME)) colMAIN <- "red"
      else colMAIN <- "black"
      TRY <- try(do.call(plot, modifyList(list(DATA[, 1], DATA[, 2], main = NAME, cex.main = 0.7, col.main = colMAIN, type = "p", 
                                               xlab = FALSE, ylab = FALSE, xaxt = "n", yaxt = "n", col = COL[i]), par2D)), silent = TRUE)
      if (inherits(TRY, "try-error")) next      
      if (!is.null(FITTED) && fitted) do.call(lines, modifyList(list(DATA[m, 1], FITTED[m], col = COL[i]), par2D))      
    }     
  }     
}  

#Calib from "qpcR"
calib <- function(
  refcurve, 
  predcurve = NULL, 
  thresh = "cpD2", 
  dil = NULL,
  group = NULL,
  plot = TRUE,
  conf = 0.95,
  B = 200
)
{
  if (class(refcurve)[1] != "modlist") stop("'refcurve' is not a 'modlist'!")
  if (!is.null(predcurve) & class(predcurve)[1] != "modlist") stop("'predcurve' is not a 'modlist'!")
  if (thresh != "cpD2" && !is.numeric(thresh)) stop("'thresh' must be either 'cpD2' or numeric!")
  if (is.null(dil)) stop("Please define dilutions!")
  if (!is.null(group) && (length(dil) != length(unique(group)))) stop("Supply as many dilutions as number of PCR groups in 'refcurve'!")
  
  lref <- length(refcurve)
  lpred <- length(predcurve)
  lgroup <- length(unique(group))
  dil <- log10(dil)
  COLref <- rep(rainbow(nlevels(as.factor(dil))), table(as.factor(dil)))
  COLpred <- rep(rainbow(lpred))   
  
  if(is.null(group))  {
    group <- as.factor(1:lref)
    isReps <- FALSE
  } else isReps <- TRUE
  
  LMFCT <- function(dil, ref, pred = NULL, conf) {
    linModY <- lm(ref ~ dil)
    conf.Y <- predict(linModY, interval = "confidence", level = conf)
    eff <- as.numeric(10^(-1/coef(linModY)[2]))
    FOM1 <- AIC(linModY)
    FOM2 <- AICc(linModY)
    FOM3 <- Rsq(linModY)
    FOM4 <- Rsq.ad(linModY)
    
    if (!is.null(pred)) {
      linModX <- lm(dil ~ ref)
      pred.conc <- sapply(as.numeric(pred), function(x) predict(linModX, newdata = data.frame(ref = x), interval = "confidence", level = conf))
    } else pred.conc <- NULL
    
    return(list(linModY = linModY, conf.Y = conf.Y, eff = eff, FOM1 = FOM1, FOM2 = FOM2,
                FOM3 = FOM3, FOM4 = FOM4, pred.conc = pred.conc[1, ], pred.conf = pred.conc[2:3, ]))
  }
  
  print("Calculating threshold time of reference curves...")
  flush.console()
  
  if (thresh == "cpD2") refCt <- sapply(refcurve, function(x) efficiency(x, plot = FALSE)$cpD2)
  else refCt <- as.numeric(sapply(refcurve, function(x) predict(x, newdata = data.frame(Fluo = thresh), which = "x")))   
  
  print("Calculating threshold time of prediction curves...")
  flush.console()
  
  if (!is.null(predcurve)) {
    if (thresh == "cpD2") predCt <- sapply(predcurve, function(x) efficiency(x, plot = FALSE)$cpD2)
    else predCt <- as.numeric(sapply(predcurve, function(x) predict(x, newdata = data.frame(Fluo = thresh), which = "x")))
  } else predCt <- NULL
  
  iterRef <- split(refCt, group)
  
  lmResList <- list()
  iterMat <- matrix(ncol = lgroup, nrow = B)
  
  for (i in 1:B) {
    if (isReps) selRef <- sapply(iterRef, function(x) sample(x, 1))
    else selRef <- unlist(iterRef)
    lmRes <- LMFCT(dil = dil, ref = as.numeric(selRef), pred = predCt, conf = conf)
    lmResList[[i]] <- lmRes
    iterMat[i, ] <- selRef
    
    if (plot) {
      if (i == 1) {
        plot(dil, selRef, col = c(2:7), pch = 16, cex = 1.3, xlab = "Log of DNA copy number", ylab = "Threshold Time (minute)",main = "Threshold value = 0.95", add = FALSE)
      } else {
        points(dil, selRef, col = c(2:7), pch = 16, cex = 1.3)
        abline(lmRes$linModY, lwd = 2)
        #Red line in the threshold graph
        lines(dil, lmRes$conf.Y[, 2], col = 2, lty = 3)
        lines(dil, lmRes$conf.Y[, 3], col = 2, lty = 3)
        
        #show the regression equation
        cf <- round(coef(lmRes$linModY), 4)
        r_sqr <- summary(lmRes$linModY)$r.squared
        eq <- paste0("y = ", cf[1],
                     ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x", ". R.square = ", 
                     round(r_sqr, 4))
        
        ## printing of the equation
        mtext(eq, 3, line=0)
      }
      if (!is.null(predcurve)) {
        points(lmRes$pred.conc, predCt, pch = 15, col = COLpred, cex = 1.5)
        if (is.vector(lmRes$pred.conf)) lmRes$pred.conf <- matrix(lmRes$pred.conf, ncol = 1)
        if (!all(is.na(lmRes$pred.conc))) {
          arrows(lmRes$pred.conf[1, ], predCt, lmRes$pred.conf[2, ], predCt, code = 3, angle = 90, length = 0.1, col = "blue")
        }
      }
    }
  }
  
  summaryList <- list()
  lenRML <- 2:length(lmResList[[1]])
  
  for (i in lenRML) {
    temp <- sapply(lmResList, function(x) x[[i]])
    summaryList[[i - 1]] <- t(temp)
  }
  
  names(summaryList) <- names(lmRes[lenRML])
  
  alpha = 1 - conf
  CONFINT <- function(x, alpha = alpha) quantile(x, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)
  
  CONF.eff <- CONFINT(summaryList$eff, alpha = alpha)
  CONF.AICc <- CONFINT(summaryList$FOM2, alpha = alpha)
  CONF.Rsq.ad <- CONFINT(summaryList$FOM4, alpha = alpha)
  
  if (!is.null(predcurve)) {
    if (nrow(summaryList$pred.conc) == 1) summaryList$pred.conc <- t(summaryList$pred.conc)
    CONF.predconc <- apply(summaryList$pred.conc, 2, function(x) CONFINT(x, alpha = alpha))
    if (!isReps) CONF.predconc <- apply(rbind(lmRes$pred.conf[1, ], lmRes$pred.conf[2, ]) , 2, function(x) CONFINT(x, alpha = alpha))
  } else {
    summaryList$pred.conc <- NULL 
    CONF.predconc <- NULL
  } 
  
  if (plot) {
    boxplot(as.numeric(summaryList$eff), main = "Efficiency", cex = 0.2)
    abline(h = CONF.eff, col = 2, lwd = 2)
    boxplot(as.numeric(summaryList$FOM2), main = "corrected AIC", cex = 0.2)
    abline(h = CONF.AICc, col = 2, lwd = 2)
    #boxplot(as.numeric(summaryList$FOM4), main = "adjusted R-square", cex = 0.2)
    #abline(h = CONF.Rsq.ad, col = 2, lwd = 2)
    if (!is.null(predcurve)) {
      boxplot(summaryList$pred.conc, main = "log(conc) of predicted", cex = 0.2)
      abline(h = CONF.predconc, col = 2, lwd = 2)
    }
  }
  return(list(eff = summaryList$eff, AICc = summaryList$FOM2, Rsq.ad = summaryList$FOM4, predconc = summaryList$pred.conc,
              conf.boot = list(conf.eff = CONF.eff, conf.AICc = CONF.AICc, conf.Rsq.ad = CONF.Rsq.ad, conf.predconc = CONF.predconc)))
}

#crop
img_crop <- function(img, x_1, x_2, y_1, y_2) {
  c_im1 <- imsub(img, x > width/x_1, y > height/y_1)
  c_im2 <- imsub(c_im1, x < width/x_2, y < height/y_2)
  return(c_im2)
  UseMethod('img_crop')
}
##Modifying build-in function of 'colocr'
roi_select <- function(img, threshold, shrink = 5, grow = 5, fill = 5,
                       clean = 5, tolerance = .1, n = 1) {
  UseMethod('roi_select')
}

#' @export
roi_select.default <- function(img, ...) {
  warning(paste("img is of class",
                class(img),
                ". img should be a cimg or a list of cimg objects."))
}

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
  
  # Ideally, I'd like to check type and value of other arguments,
  # however, I currently cannot since these arguments are optional
  #crop and rotate the image
  
  img <- img_crop(img, x_1 = 3.5, x_2 = 1.65, y_1 = 2.3, y_2 = 5)
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

#' @export
roi_select.list <- function(img, threshold, shrink = 5, grow = 5, fill = 5,
                            clean = 5, tolerance = .1, n = 1) {
  # get the length of the image list
  img_n <- length(img)
  
  # repeat arguments to match list length
  inputs <- list(threshold = threshold,
                 shrink = shrink,
                 grow = grow,
                 fill = fill,
                 clean = clean,
                 tolerance = tolerance,
                 n = n)
  
  for(i in seq_along(inputs)) {
    # lenght of argument
    input_n <- length(inputs[[i]])
    
    # use first item and return warning if not a single value or doesn't match
    # length of image list
    if(input_n != img_n & input_n != 1) {
      inputs[[i]] <- inputs[[i]][1]
      warning(paste0("Only first value in ", names(inputs)[[i]], ' will be used.'))
    }
    
    # match length of the arguments to that of the list of images
    if(input_n != img_n) {
      inputs[[i]] <- rep(inputs[[i]], img_n)
    }
  }
  
  # loop over the list of images and call roi_select
  newimgs <- list()
  for(i in 1:img_n) {
    newimgs[[i]] <- roi_select(img[[i]],
                               threshold = inputs$threshold[i],
                               shrink = inputs$shrink[i],
                               grow = inputs$grow[i],
                               fill = inputs$fill[i],
                               clean = inputs$clean[i],
                               tolerance = inputs$tolerance[i],
                               n = inputs$n[i])
  }
  
  # return list of images
  return(newimgs)
}
#Roi_show modified to show only original and the pix set images with hightlight

roi_show <- function(img, ind = c(1,3)) {
  UseMethod('roi_show')
}


roi_show.default <- function(img, ...) {
  warning(paste("img is of class",
                class(img),
                ". img should be a cimg or a list of cimg objects."))
}


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


roi_show.list <- function(img, ind = c(1,3)) {
  
  # get the length of the image list
  img_n <- length(img)
  
  # repeat argument to match list length
  if(!is.list(ind)) {
    ind <- rep(list(ind), img_n)
  }
  
  # loop over the images of lists and call roi_show
  for(i in 1:img_n){
    roi_show(img[[i]],
             ind = ind[[i]])
  }
  
  # return null
  invisible(NULL)
}

#modified function ROI_CHECK
roi_check <- function(img, ind = c(1:3)) {
  UseMethod('roi_check')
}

roi_check.default <- function(img, ...) {
  warning(paste("img is of class",
                class(img),
                ". img should be a cimg or a list of cimg objects."))
}


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



# Define UI for application that draws a histogram
ui <- navbarPage(
  
  # Application page Title
  title = 'ReaLAMP',

  
  ## Main
  tabPanel('Main',
           sidebarLayout(
             sidebarPanel(
               # Input Panel
               
               tags$h3('Input Panel'),
               
               ## Description
               
               ## load Image
               fileInput('image1', 'Upload Images', multiple = TRUE),
               bsTooltip('image1',
                         'Select and upload the chamber images.',
                         'right', options = list(container = "body")),
               numericInput('roi_num', 'Number of Chamber', 7, 1, 50, 1),
               bsTooltip('roi_num',
                         'Select the number of Sample on a chip. Default is 7',
                         'right', options = list(container = "body")),
               tags$hr(),
               #tags$p(' Crop images'),
               ## Selection Parameters
               #numericInput('x1', 'X', 5, min = 1, max = 10),
               #numericInput('x2', 'X 2', 1.3, min = 1, max = 10),
               #numericInput('y1', 'y', 1.9, min = 1, max = 10),
               #numericInput('y2', 'y 2', 5, min = 1, max = 10),
               #tags$hr(),
               tags$p('Adjust parameters'),
               sliderInput('threshold', 'Threshold', 1, 99, 29, 1),
               bsTooltip('threshold',
                         'Choose a threshold for excluding the background.',
                         'right', options = list(container = "body")),
               sliderInput('shrink', 'Shrink', 1, 10, 1, 1),
               bsTooltip('shrink',
                         'Shrink the selected area by eroding the bounderies around it.',
                         'right', options = list(container = "body")),
               sliderInput('grow', 'Grow', 1, 10, 1, 1),
               bsTooltip('grow',
                         'Grow the selected area by dilating the bounderies around it.',
                         'right', options = list(container = "body")),
               sliderInput('fill', 'Fill', 1, 10, 1, 1),
               bsTooltip('fill',
                         'Remove holes in the selected area.',
                         'right', options = list(container = "body")),
               sliderInput('clean', 'Clean', 1, 10, 1, 1),
               bsTooltip('clean',
                         'Remove small isolated parts in the selected area.',
                         'right', options = list(container = "body")),
               sliderInput('tolerance', 'Tolerance', 0, .99, .1, .1),
               bsTooltip('tolerance',
                         'Set value to determine which two neighboring pixels are in same selected area.',
                         'right', options = list(container = "body"))
               
               ),
             mainPanel(
               # Output Views
               #
               # This part of the app contains the different views of the output
               # It is divided into description and four tabs. The four tabs are:
               #    1. Select ROI
               #    2. Pixel Intensities
               #    3. Graph View
               
               fluidRow(
                 
                 ## Description
                 tags$h2('Colorimetric Real time LAMP'),
                 tags$br(),
                 tags$li('Step 1: Upload the images'),
                 tags$li('Step 2: Adjust parameters in the input panel. The auto-selected chambers show in Preview tab with red borders'),
                 tags$li('Step 3: The graphical results show in ReaLAMP analysis tab, the raw data show in Data table tab'),
                 tags$br(),
                 tags$br(),
                 
                 # Tabs
                 tabsetPanel(
                   
                   ## Select ROI
                   tabPanel('Preview',
                            plotOutput("image_plot", height = "6000px")
                           
                   ),
                   
                   ## RGB ratio
                   
                   tabPanel('ReaLAMP analysis',
                  #sliderInput('thresh', 'Threshold value', 0,1.5, 0.95, 0.05),
                            plotOutput('graph_plot', height = "1000")),
                  
            
                   
                  
                   
                   # Graph View
                   tabPanel('Data tables',
                            tableOutput('table'))
                 )
                 )
                 )
             )),
  
  # GitHub
  tabPanel('GitHub',
           "Comments, issues and contributions are welcomed")
  )


# Define server
server <- function(input, output) {
  # intiate interactive values
  values <- reactiveValues()
  
  # load images
  img1 <- reactive({
    image_load(input$image1$datapath)
             
  })
  
  ## calculate the pixset
  px <- reactive({
    roi_select(img1(),
               threshold = input$threshold,
               shrink = input$shrink,
               grow = input$grow,
               fill = input$fill,
               clean = input$clean,
               n = input$roi_num)
  })
  
  ## calculate correlations
  corr <- reactive({
    roi_test(px(), type = 'both')
  })
  b1 <- reactive({
    roi_check(px())  
  })
  # Output Views
  ## Select ROI
  
  # plots
  output$image_plot <- renderPlot({
    req(input$image1)
    
    n <- length(input$image1$name)
    
    par(mfrow=c(n,4), mar = rep(1, 4))
    roi_show(px())
  })
  
  
  ## sigmode plot
  output$graph_plot <- renderPlot({
    req(input$image1)
    par(mfrow=c(4,1))
    b1()
  })
  ## G/R ratio table
  output$table <- renderTable({
    req(input$image1)
    b1()
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
