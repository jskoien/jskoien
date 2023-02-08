
eventmeans = function(allvar, allobs, dels = 0, plotdir, mfrow = c(2,2), ranges, first = FALSE) {

  if (missing(allobs)) allobs = do.call(rbind, allvar)
  if (missing(plotdir)) plotdir = getwd()
  allparcels = unique(allobs$ID)
  inds = names(allvar)
  
  if (missing(ranges)) {
    ranges = matrix(NA, ncol = 2, nrow = length(inds))
    for (iind in 1:length(inds)) {
      ranges[iind, ] = range(allvar[[iind]]$MEAN_W, na.rm = TRUE)*c(0.95, 1.05)
    }
    sinds = grep("s1_coh6", inds)
    ranges[sinds,1] = min(ranges[sinds,])
    ranges[sinds,2] = max(ranges[sinds,])
  }
  
  if (length(dels) > 1){
    pdf(paste0(plotdir, "\\delrange", country, ".pdf"))
    pzero = dev.cur()
    par(mfrow = mfrow)
  } 
  pdf(paste0(hdir, "\\summaryplots_", country, ".pdf"))
  par(mfrow = mfrow)
  pfirst = dev.cur()
  pdf(paste0(hdir, "\\ttest_", country, ".pdf"))
  par(mfrow = mfrow)
  pthird = dev.cur()
  

  delnums = array(NA, dim = c(10,length(dels), length(weeks)+1))
  
  gmeans = list()
  stds = list()
  tabs = list()
  geventmeans = list()
  agdiffs = list()
  for (ievent in (length(eventslist)+1):1) {
    if (ievent > length(eventslist)) event = "all" else event = eventslist[ievent]
    geventmeans[[ievent]] = list()
    
    mevent = if (ievent > length(eventslist)) paste(country, "all") else 
      paste(country, paste(eventslist[[ievent]], collapse = " - "))
    gdiffs = matrix(NA, nrow = length(allparcels), ncol = 12)
    for (iind in 1:length(inds)) {
      ind = inds[iind]
      dat = allvar[[iind]]
      
      dat$eventnum = NA
      parcels = unique(dat$ID)
      gdd = data.frame(allparcels, diffs = NA)
      for (iparc in 1:length(parcels)) {
        parc = parcels[iparc]
        days = sort(unique(dat$DOY_[dat$ID == parc]))
        for (ii in 1:length(days)) {
          dat$eventnum[dat$ID == parc & dat$DOY_ == days[ii]] = ii
        }
      }
      if (first) dat = dat[dat$eventnum == 1, ]
      
      
      geventmeans[[ievent]][[iind]] = list()
      
      for (iidel in 1:length(dels)) {
        ldat = dat
        idel = dels[iidel]
        delrange = c(-idel, idel)
        if (idel == 0 & !exists("fixrange")) {
          ranges[iind, ] = range(ldat$MEAN_W, na.rm = TRUE)*c(0.95,1.05)
        } else if (idel > 0) {
          for (parc in unique(ldat$ID)) {
            pdat = ldat[ldat$ID == parc, ]
            days = unique(pdat$DOY_)
            for (iday in days) {
              if (any(findInterval(days - iday, c(-1e6,delrange[1],-.5,.5,delrange[2], 1e6)) %in% c(2,4))) {
                ldat$VALID_OBS_W[ldat$ID == parc & ldat$DOY_ == iday] = 0
              }
            }
          }
          ldat = ldat[ldat$VALID_OBS_W == 1, ]
        }  
        
        for (iparc in 1:length(allparcels)) {
          parc = parcels[iparc]
          ids = which(ldat$ID == parc)
          if (length(ids) == 0) next
          lldat = ldat[ids,]
          lgd = try(lldat$MEAN_W[lldat$WEEK == (week0+2)] - lldat$MEAN_W[lldat$WEEK == week0], silent = TRUE)
          if (length(lgd) > 0 & !is(lgd, "try-error")) gdd[iparc,2] = lgd
        }
        
        if (ievent <= length(eventslist)) {
          ldat = ldat[ldat$EVENT_TYPE %in% eventslist[[ievent]], ]
          #      print(paste(country, ievent, dim(ldat[1])))
          if (dim(ldat)[1] < 2*max(weeks)) { 
            if (ind %in% c("ndvi", "s1_coh6")) {
              if (country == "IE" & ind == "ndvi") next
              #          if (country == "IE") stop()
              if (length(dels) > 1) {
                dev.set(pzero)
                #if (country == "IE") stop()
                plot(1,1, ylab = ind, main = mevent)
              }
            }
            next
          }
        }
        gmean = aggregate(ldat$MEAN_W, by = list(ldat$WEEK), mean)
        if (weighted) {
          ldat$wsum = ldat$MEAN_W*ldat$COUNT_W
          gmean = aggregate(ldat$wsum, by = list(ldat$WEEK), sum)
          gweight = aggregate(ldat$COUNT_W, by = list(ldat$WEEK), sum)
          gmean$x = gmean$x/gweight$x
        } else gmean = aggregate(ldat$MEAN_W, by = list(ldat$WEEK), mean)
        
        
        
        if (ievent > length(eventslist)) gmeans[[iind]] = gmean
        gmed = aggregate(ldat$MEAN_W, by = list(ldat$WEEK), median)
        gnum = aggregate(rep(1, length(ldat$MEAN_W)), by = list(ldat$WEEK), sum)
        gstd = aggregate(ldat$MEAN_W, by = list(ldat$WEEK), sd)
        if (ievent > length(eventslist)) stds[[iind]] = gstd
        geventmeans[[ievent]][[iind]][[iidel]] = list(gmean = gmean, gstd = gstd)
        gg = gmean$x
        rang = range(c(ldat$MEAN_W, gmeans[[iind]]$x))
        #       print(paste(country, ievent, dim(ldat[1]), rang[2], max(gmeans[[iind]]$x)))
        if (idel == 0 | TRUE) {
          dev.set(pfirst)
          sweeks = ldat$WEEK %in% weeks
          plot(weeks-week0, gmean$x[weeks], ylim = rang, 
               xlab = "weeks", ylab = inds[[iind]], main =  paste(mevent, idel, "test"))
          lines(weeks - week0, gmean$x[weeks])
          lines(weeks - week0, gmeans[[iind]]$x[weeks], col = "red")
          abline(v = 0)
          if (country == "IE") {
            lines(weeks - week0,allmeans[["CZ"]][[iind]]$x[weeks], col = "green")
          }
          lines(weeks - week0, gmed$x[weeks], col = "blue")
          text(gnum$Group.1, mean(range(gg)), labels = gnum$x)  
          boxplot(MEAN_W ~ WEEK, ldat[sweeks,], main = inds[[iind]])
          
          dev.set(psecond)
          #if (country == "LV" & iind == 2) debug(findQuad)
          #debug(plotConfidence)
          levent = paste(mevent, idel)
          ylim = plotConfidence(gmean$x, stdobs = gstd$x, weeks = gmean$Group.1, 
                                week0 = week0, xlab = "week",
                                ylab = ind, main = levent, ylim = ranges[iind, ] )
          if (idel == 0) ranges[iind, ] = ylim
          abline(v = 0)
          abline(h = 0)
          text(gnum$Group.1-week0, mean(range(gg)), labels = gnum$x)  
          
          tab = array(NA, dim = c(max(weeks), max(weeks),4))
          for (ii in weeks[1:(length(weeks) - 1)]) {
            for (jj in (ii+1): max(weeks)) {
              xdat = ldat$MEAN_W[ldat$WEEK == ii]
              ydat = ldat$MEAN_W[ldat$WEEK == jj]
              if (length(xdat) > 2 & length(ydat) > 2) {
                tt = t.test(xdat, ydat)
                tab[ii,jj,1] = tt$conf.int[1] 
                tab[ii,jj,2] = tt$conf.int[2]
                tab[ii,jj,3] = -diff(tt$estimate)
                tab[ii,jj,4] = tt$p.value
              }
            }
          }
          tabs[[iind]] = tab
          
          
          #        print(inds[iind])
          nweeks = max(weeks)
          ltab = matrix(0,ncol = 4, nrow = nweeks)
          for (ii in c(1:nweeks)) {
            #     ltab[ii,] = c(tt[ii+1, ii,1], tt[ii, ii+1,1], tt[ii,ii+1,2])
            if (ii < week0) {
              ltab[ii,] = c( tab[ii, week0, 1], tab[ii, week0,  2], tab[ii, week0, 3], tab[ii, week0, 4])
            } else if (ii > week0) {
              ltab[ii,] = c(-tab[week0, ii, 1], -tab[week0, ii, 2], -tab[week0, ii, 3], tab[week0, ii, 4])
            }
          }
          df0 = data.frame(var = inds[iind], signif(ltab, 2))
          df0[week0,5] = 1
          if (iind == 1) {
            df = df0
          } else {
            df = rbind(df, df0 )
          }
          dev.set(pthird)
          if (sum(!is.na(df0[,2:4])) <= 4) {
            print(paste("noprint", country, ievent, iind, sum(!is.na(df0[,2:4]))  ))
            #           if (country == "IE" & iind == 5) stop()
            plot(1,1)
          } else {
            if (idel > 0) levent = paste(mevent, idel) else levent = mevent 
            if (length(grep("s1_coh", ind)) > 0) {
              ylim = c(-0.02, 0.1)
              if (max(df0[,3], na.rm = TRUE) > .1) ylim[2] = 0.15 
            } else ylim = range(df0[,2:4], na.rm = TRUE)
            plot(weeks-week0, df0[,4], ylim = ylim, t = "l", 
                 main = levent, xlab = "week", ylab = inds[[iind]])
            lines(weeks - week0, df0[,2], col = "blue")
            lines(weeks - week0, df0[,3], col = "blue")
            lcor = findQuad(df0[,2:3], week0, weeks, yrang = ylim)
            points(weeks - week0, df0[,4], pch = 16, col = gray(sqrt(df0[,5])), cex = 2)
            abline(h = 0)
            legend(lcor[1], lcor[2], legend = c("0.02", "0.05", "0.1", "0.2", "0.5", "0.8"),
                   pch = 16, col = gray(sqrt(c(0.02, 0.05, 0.1, 0.2, 0.5, 0.8))), pt.cex = 2)
          }
          
        }
        print(paste("ccc: ", country, ind, idel, ievent))
        if (ind %in% c("ndvi", "s1_coh6") & !(ind == "ndvi" & country == "IE")) {
          #          if (ind %in% c("ndvi", "s1_coh6"))  {
          print(paste("cc1: ", country, ind, idel, ievent))
          #          if (!(idel %in% c(0,21)) | !(idel == 9 & ievent == 7)) next
          if (idel %in% c(0,21) | (idel == 9 & ievent == 7)) {
            print(paste("cc2: ", country, ind, idel, ievent))
            if (idel == 21 & ievent == 7) next
            print(paste("cc3: ", country, ind, idel, ievent))
            
            if (length(dels) > 1) {
              dev.set(pzero)
              plotConfidence(gmean$x, stdobs = gstd$x, weeks = gmean$Group.1, 
                             week0 = week0, xlab = "week",
                             ylab = ind, main = paste(mevent, idel, sep = " - "), ylim = ranges[iind, ])
              abline(v = 0)
              text(gnum$Group.1-week0, mean(range(gg)), labels = gnum$x)  
              if (ind == "s1_coh6") delnums[[country]][ievent, iidel, 1:length(weeks)] = gnum$x
              if (ind == "s1_coh6") delnums[[country]][ievent, iidel, length(weeks) + 1] = mean(gnum$x)
            }
          }
        }
      }
      # 800303303
      gdiffs[,iind] = gdd[,2]
    }
    agdiffs[[ievent]] = gdiffs
  }
  #  plot(1:100,1:100, pch = 16, col = gray(sqrt(seq(0.01, 1, .01))))
  
  
  #  grid.table(df[1:24,])
  #  grid.newpage()
  #  grid.table(df[25:48,])
  dev.off(pfirst)
  dev.off(psecond)
  dev.off(pthird)
  if (!names(pzero) == "null device") dev.off(pzero)
  allranges[[country]] = ranges
  cgeventmeans[[country]] = geventmeans
  allmeans[[country]] = gmeans
  cagdiffs[[country]] = agdiffs
}
