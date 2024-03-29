rtopKrige.rtop = function(object, varMatUpdate = FALSE, params = list(), ...) {
  params = getRtopParams(object$params, newPar = params, ...)
  if (!is.null(params$nsim) && params$nsim > 0) return(rtopSim(object, varMatUpdate, params = params))
  observations = object$observations
  
  predictionLocations = object$predictionLocations
  if (!all(c("varMatObs", "varMatPredObs") %in% names(object)) | varMatUpdate) 
    object = varMat(object, varMatUpdate, params = params,  ...)
  
  varMatObs = object$varMatObs
  varMatPredObs = object$varMatPredObs
  
  krigeRes = rtopKrige(object = observations, predictionLocations = predictionLocations, 
                       varMatObs = varMatObs, 
                       varMatPredObs = varMatPredObs, params = params, 
                       formulaString = object$formulaString, ...)
  object$predictions = krigeRes$predictions
  if ("cvInfo" %in% names(krigeRes)) object$cvInfo = krigeRes$cvInfo
  if ("weight" %in% names(krigeRes)) object$weight = krigeRes$weight
  if ("removed" %in% names(krigeRes)) object$removed = krigeRes$removed
  object
}  


rtopKrige.SpatialPolygonsDataFrame = function(object, predictionLocations = NULL,
                                              varMatObs, varMatPredObs, varMat, params = list(), formulaString,  
                                              sel, ...) {
  rtopKrige.default(object, predictionLocations, varMatObs, 
                    varMatPredObs, varMat, params, formulaString,  
                    sel, ...) 
}



rtopKrige.default = function(object, predictionLocations = NULL,
                             varMatObs, varMatPredObs, varMat, params = list(), formulaString,  
                             sel, wret = FALSE, ...) {
  params = getRtopParams(params, ...)
  if (!is.null(params$nsim) && params$nsim > 0) 
    return(rtopSim(object, predictionLocations,
                  varMatObs, varMatPredObs, varMat, params, formulaString, ...))
  #
  cv = params$cv  
  #  else object$params$cv = params$cv = cv
  nmax = params$nmax    
  wlim = params$wlim
  wlimMethod = params$wlimMethod
  maxdist = params$maxdist
  debug.level = params$debug.level
  lambda = params$lambda
  if ("singularSolve" %in% names(params)) singularSolve = params$singularSolve else singularSolve = FALSE
  if (!is.null(lambda)) BLUE = TRUE else BLUE = FALSE
  if (!missing(varMat)) {
    if (missing(varMatObs)) varMatObs = varMat$varMatObs
    if ((missing(varMatPredObs) | is.null(varMatPredObs)) && !cv) varMatPredObs = varMat$varMatPredObs
  }
  depVar = as.character(formulaString[[2]])
  observations = object
  obs0 = observations[[depVar]]
  nobs = dim(object)[1]
  if (inherits(observations, "sf")) {
    obscors = st_coordinates(observations)
    if (cv) newcors = obscors else newcors = suppressWarnings(st_coordinates(st_centroid(predictionLocations)))
  } else {
    obscors = coordinates(observations)
    if (cv) newcors = obscors else newcors = coordinates(predictionLocations)
  }
  
  npred = ifelse(cv,nobs,ifelse(!missing(sel),length(sel),dim(newcors)[1]))
  if (missing(sel)) sel = c(1:npred)
  if (params$unc && "unc" %in% names(observations)) {
    unc0 = observations$unc
  } else unc0 = array(0,nobs)
  #  
  if (inherits(observations, "Spatial")) {
    mdist = sqrt(bbArea(bbox(observations)))
  } else mdist = sqrt(bbArea(st_bbox(observations)))
  if (nobs < nmax && mdist < maxdist && !cv) {
    varMat = rbind(varMatObs,1)
    diag(varMat) = unc0
    varMat = cbind(varMat,1)
    varMat[nobs+1,nobs+1] = 0    
    varInv = solve(varMat)    
    singMat = TRUE 
  } else singMat = FALSE
  #
  if (cv) {
    predictionLocations = observations
    varMatPredObs = varMatObs
  }
  # 
  removed = NULL
  predictions = data.frame(var1.pred = rep(0,npred),var1.var = 0,sumWeights = 0)
  if (cv) {
    predictions = cbind(predictions,observed = observations[[depVar]], residual=0, zscore = 0)
    cvInfo0 = data.frame(inew = c(0), obs = c(1), pred = c(0), var = c(0),
                         residual = c(0), zscore = c(0), neigh = c(0), neighobs = c(0),
                         weight = c(0), sumWeight = c(0), c0arr = c(0))
    cvInfo = list()
  }
  #
  
  if (wret) weight = matrix(0,nrow = npred,ncol = nobs)
  if (interactive() & debug.level == 1 & length(sel) > 1) {
    pb = txtProgressBar(1, length(sel), style = 3)
  }
  
  if (debug.level >= 1) print(paste(ifelse(cv, "cross-validating", "interpolating "), length(sel), "areas"))

    for (inn in 1:length(sel)) {
    inew = sel[inn]
    if (debug.level > 1) print("\n")
    #  for (inew in 1:20) {         
    if (interactive() & debug.level == 1 & length(sel) > 1) setTxtProgressBar(pb, inn)
                                          
    if (cv) {
      if (debug.level > 1) print(paste("Cross-validating location", inew, 
                                       " out of ",npred," observation locations"))
      if (debug.level > 1) print(data.frame(observations)[inew,] )
      #      if (cv == inew && inew > 1) browser()
    } else {
      if (debug.level > 1) print(paste("Predicting location ",inew,
                                       " out of ", npred," prediction locations" ))
      if (debug.level > 1 && is(predictionLocations, "SpatialPolygonsDataFrame")) 
        print(predictionLocations@data[inew,] )
    }
    newcor = newcors[inew,]
    
    ret = rkrige(data.frame(observations), obs0, obscors, newcor, varMatObs, 
                 varMatPredObs[,inew], nmax, inew, cv, 
                 unc0, mdist, maxdist, singMat, varInv, singularSolve, 
                 wlim, debug.level, wlimMethod, BLUE)
    
    predictions$var1.pred[inew] = ret$pred[1]
    predictions$var1.var[inew] = ret$pred[2]
    predictions$sumWeights[inew] = ret$pred[3]
    
    nneigh = ret$nneigh
    lambda = ret$lambda
    neigh = ret$neigh
    if (wret) weight[inew,neigh] = lambda[1:(length(lambda)-1)] 
    obs = ret$obs
    unc = ret$unc
    c0arr = ret$c0arr
    if (debug.level >1) {
      distm = spDistsN1(obscors,newcor)[neigh]
      lobs = observations@data[neigh,]
      lobs = rbind(lobs, mu =  rep(0, (dim(lobs)[2])))
      lobs = cbind(lobs, data.frame(id = c(neigh, 0), 
                                    edist = c(distm, 0), lambda = lambda, c0 = c0arr, 
                                    obs = c(obs, 1), unc = c(unc, 0), 
                                    lambda_times_obs = lambda*c(obs, 0))) 
      print("neighbours")
      print(lobs, 3)
      print("covariance matrix ")
      print(varMat,3)     
    }
    if (cv) {
      predictions$residual[inew] = observations[[depVar]][inew]-predictions$var1.pred[inew]
      predictions$zscore[inew] = predictions$residual[inew]/sqrt(predictions$var1.var[inew])
      cvNew = data.frame(predictions[inew,],c0arr=lambda[nneigh+1],neigh=0)
      for (jnew in 1:nneigh) cvNew = rbind(cvNew,
                                           data.frame(var1.pred = 0,
                                                      var1.var = unc[neigh[jnew]], sumWeights=lambda[jnew],
                                                      observed = data.frame(observations)[[depVar]][neigh[jnew]], residual = 0, zscore = 0,
                                                      c0arr=c0arr[jnew],neigh=neigh[jnew]))
      if (inew == 1) cvInfo = cvNew else cvInfo = rbind(cvInfo,cvNew)
      if (debug.level > 1) {
        print("prediction") 
        print(c(data.frame(predictionLocations)[inew,],predictions[inew,]))
      }
    }
    if (singularSolve & !is.null(ret$removed)) {
      if (!is.null(removed)) removed = rbind(removed, ret$removed) else removed = ret$removed
    }
  }  
  if (interactive() & debug.level == 1 & length(sel) > 1) close(pb)
  if (inherits(predictionLocations, "SpatialPolygons")) {
      if ("data" %in% names(getSlots(class(predictionLocations)))) {
      predictionLocations@data = cbind(predictionLocations@data, predictions)
      predictions = predictionLocations
    } else {
      predictions = addAttrToGeom(predictionLocations, predictions, match.ID = FALSE)
    }
  } else {
    predictions = cbind(predictionLocations, predictions)
  }
  if (cv) predictions$observed = observations[[depVar]]
  ret = list(predictions = predictions)
  if (wret) ret$weight = weight
  if (singularSolve & !is.null(removed)) ret$removed = removed
  if (cv) ret$cvInfo = cvInfo
  ret
}









