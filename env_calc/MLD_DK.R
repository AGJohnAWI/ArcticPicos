#credit: https://github.com/Davidatlarge/mixed_layer_depth/blob/master/mixed_layer_depth_DK.R

MLD_DK <- function(depth, # vector for depth, can be negative
                   variable, # temperature or density
                   threshold, # defines the threshold value that marks the end of the MLD
                   print.info = FALSE, # should into be printed in the console als a side effect
                   plot = FALSE, # should the result be plotted as a side effect
                   depth.max = NULL # maximum depth to be plotted, uses absolute value of max depth if no value is supplied
)
{
  data <- data.frame(depth, variable)
  
  data$delta <- abs(unique(data$variable[data$depth==min(abs(data$depth))]) - data$variable) # unique() in case the min(depth) occures more than once
  
  threshold <- abs(threshold)
  
  # find MLD
  if(any(data$delta>=threshold)) #  if the threshold is reached at any depth ...
  {
    MLD <- data$depth[data$delta>=threshold][1] # ... find the first depth where the threshold is exceeded
  } else { # ... otherwise (if water is stratified) set MLD <- NA
    MLD <- NA
  }
  
  # print info
  if(print.info) {
    if(plot&is.null(depth.max)) {
      print("no depth.max supplied, using maximum depth instead")
    }    
    if(is.numeric(MLD)) {
      print(paste("mixed layer depth is", MLD, "[units of depth]", sep = " ")) 
    } else {
      print(paste("The thereshold is not exceeded. The water column is well mixed to the maximum depth of", depth[length(depth)],"[units of depth]",sep=" "))
    }
    
  }
  
  # plot result
  if(plot) {
    if(is.null(depth.max)) {
      depth.max <- depth[length(depth)]
    }
    # mixed
    plot(variable, depth, 
         ylim = c(depth.max, 0), 
         xlim = c(min(variable)-threshold*2, max(variable)+threshold*2)
    )
    abline(h = MLD, col = "blue")
    text(x = min(variable), y = MLD, paste("MLD =", MLD), col = "blue", adj = c(-0.2,-0.5))
    abline(v = c(variable[depth==min(abs(depth))]+threshold, variable[depth==min(abs(depth))]-threshold), col = "red")
    text(x = mean(variable), y = min(abs(depth)), "threshold", col = "red") 
  }
  
  return(MLD)
  
}