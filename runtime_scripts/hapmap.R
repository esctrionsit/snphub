hm_main <- function(co, ro, mis, lon, lat, mer, scal){
  withProgress(message = 'Drawing', detail = "  Drawing HapMap...", value = 5, {
    library(ggplot2)
    library(ggmap)

    co <- paste(strsplit(co, split = " ")[[1]], collapse="")
    fra_hm_detail <<- data.frame()
    if(is.na(as.numeric(scal)) || as.numeric(scal) < 0) { return(hm_error_message(6)) }

    text_hm_NAwarning <<- ""

    hm_chr <- pub_check_chr_value(ro)

    if(grepl("-", ro)){
        hm_beg <- pub_check_pos_begin(ro)
        hm_end <- pub_check_pos_end(ro)
    }else{
        hm_pos <- hm_check_pos_value(ro)
        hm_beg <- hm_pos
        hm_end <- hm_pos
    }    

    if(sum(is.na(hm_chr), is.na(hm_beg), is.na(hm_end)) > 0 || length(hm_chr) == 0) { return(hm_error_message(2)) }
    if(hm_range_is_too_long(hm_beg, hm_end)) { return(hm_error_message(3)) }

    hm_sam <- pub_check_sample_name(c(co))

    if(sum(is.na(hm_sam)) > 0 || length(hm_sam) == 0) { return(hm_error_message(4)) }

    hm_sam <- pub_unique_sample(hm_sam)

    fet_code <- hm_fetch_data(hm_chr, hm_beg, hm_end, hm_sam)

    if(fet_code != 0) { return(hm_error_message(fet_code + 100)) }

    hm_fra_extra <- fra_hm_orivcf[, 3:4]
    fra_hm_orivcf <<- fra_hm_orivcf[, c(-3,-4)]

    tran_code <- hm_trans_data(hm_sam, mis, scal, hm_fra_extra)

    int_hm_warning <<- tran_code

    NAname_lis <- hm_recheck_data()

    if(length(NAname_lis)!=0){
        text_hm_NAwarning <<- paste("Warning: Sample ", paste(NAname_lis, collapse=","), " is/are ignored for don't having location information.")
    }

    if(nrow(fra_hm_detail) == 0){ return(hm_error_message(102)) }
      
    p <- hm_draw_plot(lon, lat, mer, co, ro)
    p
  })
}

hm_trans_data <- function(sam, mis, scal, hm_fra_extra){
    stat <- 0
    fac <- as.numeric(scal)
    text_hm_drawsite <<- paste(paste(as.character(fra_hm_orivcf[1,1]), as.character(fra_hm_orivcf[1,2]), sep=":"), ", Ref: ", as.character(hm_fra_extra[1,1]), ", Alt: ", as.character(hm_fra_extra[1,2]), sep="")
    fra_hm_orivcf <<- fra_hm_orivcf[,-c(1,2)]
    res <- data.frame(sam, rep(0, time=length(sam)), rep(0, time=length(sam)), rep(0, time=length(sam)))
    names(res) <- c("Acession", "REF", "ALT", "MISSING")
    if(nrow(fra_hm_orivcf)>1){
        stat <- 1
    }
    for(i in 1:ncol(fra_hm_orivcf)){
        tmp <- as.character(fra_hm_orivcf[1, i])
        s <- strsplit(tmp, split = "/")[[1]]
        for(s_tmp in s){
            if(mis == "Yes" && s_tmp == "."){
                res[i, 4] = res[i, 4] + fac
            }else if(s_tmp == "0"){
                res[i, 2] = res[i, 2] + fac
            }else if(s_tmp != "."){
                res[i, 3] = res[i, 3] + fac
            }
        }
    }

    fra_hm_drawdata <<- res
    #return
    stat
}

hm_recheck_data <- function(){
    NAlis <- c()
    NAname_lis <- c()
    for(i in 1:nrow(fra_hm_drawdata)){
        tmp <- fra_glo_samloc[which(fra_glo_samloc$Acession==fra_hm_drawdata[i,1]),]
        if(nrow(tmp)!=1 || is.na(tmp[1,3]) || is.na(tmp[1,4])){
            NAlis <- c(NAlis, i)
            NAname_lis <- c(NAname_lis, as.character(fra_hm_drawdata[i,1]))
        }
    }
    if(length(NAlis)>0){
        fra_hm_drawdata <<- fra_hm_drawdata[-NAlis,]
    }
    NAname_lis
}

hm_draw_plot <- function(ilon, ilat, imer, co, ro){
    p <- tryCatch({
        dff <- merge(fra_glo_samloc, fra_hm_drawdata, by = 'Acession')

        fra_hm_detail <<- dff
        var_type_lis <- c()
        for(i in 1:nrow(fra_hm_detail)){
            if(as.numeric(fra_hm_detail[i,5])==0 && as.numeric(fra_hm_detail[i,6])==0){
                var_type_lis <- c(var_type_lis, "MISSING")
            }else if(as.numeric(fra_hm_detail[i,5])==0 && as.numeric(fra_hm_detail[i,7])==0){
                var_type_lis <- c(var_type_lis, "HOMO")
            }else if(as.numeric(fra_hm_detail[i,6])==0 && as.numeric(fra_hm_detail[i,7])==0){
                var_type_lis <- c(var_type_lis, "HOMO")
            }else if(as.numeric(fra_hm_detail[i,7])==0){
                var_type_lis <- c(var_type_lis, "HETE")
            }
        }
        fra_hm_detail <<- cbind(fra_hm_detail[,c(1,2)], var_type_lis)
        names(fra_hm_detail) <<- c("Accession", "Location", "Mutation_Type")

        # disthr should change with size of map to avoid points overlaping
        dff <- hm_draw_mergeloc(dff, disthr = imer)
        dff$radius <- log(apply(dff[,-c(1:4)], 1, function(x) sum(x)) +1)/3
        pie <- dff[,-c(1,2)]

        # ggplot map data
        world = map_data("world", resolution=0)

        t = read.table(pipe("date +%Y%m%d%H%M%S"))
        d = as.character(t[1,1])
        filenametag <- paste("SnpHub_HapMap_", d, sep="")
        parameter <- paste("Parameter: ", co, "; ", ro,"; ", sep="")

        p <- ggplot(data = world, aes(x=long, y=lat, group=group)) + 
            geom_polygon(fill = "white", color = "#A9D3EB") + 
            coord_quickmap(xlim = ilon, ylim = ilat) +
            labs(title = "Haplotype map", caption = paste0(filenametag, "\n", parameter)) +
            ylab("Latitude") + 
            xlab("Longitude") + 
            theme(
                panel.background = element_rect(fill = "#A9D3EB"),
                panel.grid.minor = element_blank(), 
                panel.grid.major = element_line(colour = "grey90", size = 0.5), 
                legend.position = "top",
                legend.direction = "horizontal",
                legend.title = element_blank(),
                legend.spacing.x = unit(0.01, "npc"),
                legend.justification = "right",
                plot.title = element_text(size = 25,hjust=0.5),
                axis.title = element_text(size=20),
                strip.text = element_text(size=15),
                plot.margin = grid::unit(c(0.05,0.05,0.05,0.05), "npc"))

        pie.list <- pie %>% 
            tidyr::gather(type, value, -lon, -lat, -radius) %>%
            tidyr::nest(type, value) %>%
            
            # make a pie chart from each row, & convert to grob
            mutate(pie.grob = purrr::map(data,
                function(d) ggplotGrob(ggplot(d, 
                    aes(x = 1, y = value, fill = type)) +
                    geom_col(color = "black",
                        show.legend = FALSE) +
                    coord_polar(theta = "y") +
                    theme_void()))) %>%
            
            # convert each grob to an annotation_custom layer. I've also adjusted the radius
            # value to a reasonable size (based on my screen resolutions).
            rowwise() %>%
            mutate(radius = radius * 4) %>%
            mutate(subgrob = list(annotation_custom(grob = pie.grob,
                xmin = lon - radius, xmax = lon + radius,
                ymin = lat - radius, ymax = lat + radius)))

        p <- p + 
            
            # Optional. this hides some tiles of the corresponding color scale BEHIND the
            # pie charts, in order to create a legend for them
            geom_tile(data = pie %>% tidyr::gather(type, value, -lon, -lat, -radius),
                aes(x = lon, y = lat, fill = type), 
                color = "black", width = 0.01, height = 0.01, 
                inherit.aes = FALSE) +
            pie.list$subgrob

    }, error = function(e) {
        paste("202: Unknown draw time errror." , e, collapse = "\n")
    })
    p
}

hm_draw_mergeloc <- function(locdf, disthr = 0.1){

  # merge same location, ggmap use location as group
  if(disthr < 0.1){disthr = 0.1}
  
  newlocdf <- locdf[FALSE,]
  for (i in seq_len(nrow(locdf))){
    if (i == 1){
      newlocdf <- rbind(newlocdf, locdf[i,])
    }else{
      locnum <- nrow(newlocdf)
      n = 0
      for (j in seq_len(nrow(newlocdf))){
        if (as.numeric(dist(rbind(locdf[i,c(3,4)], newlocdf[j,c(3,4)]))) < disthr){
          newlocdf[j,-c(1:4)] <- newlocdf[j,-c(1:4)] + locdf[i,-c(1:4)]
          break
        }else{
          n=n+1
        }
      }
      if (n == locnum){
        newlocdf <- rbind(newlocdf, locdf[i,])
      }
    }
  }
  return(newlocdf)
}

hm_range_is_too_long <- function(hm_beg, hm_end){ F }

hm_error_message <- function(code){
    code
}

hm_check_pos_value <- function(ro) {
    inputs <- ro
    res <- c()
    for(i in inputs){
        if (grepl(":", i)){
            tmp <- strsplit(i, split = ":")[[1]]
            if(length(tmp)>1){
                if(!is.na(as.numeric(tmp[2]))) {
                    res <- c(res, tmp[2])
                }else {
                    res <- c(res, NA)
                }
            }else{
                res <- c(res, NA)
            }
        }
    }
    #return
    res
}
