shinyServer(function(input, output, session){

	reaobj = reactiveValues()
	reaobj$snp_stat <- "System Info: Ready"
	reaobj$hn_stat <- "System Info: Ready"
	reaobj$hp_stat <- "System Info: Ready"
	reaobj$dt_stat <- "System Info: Ready"
	reaobj$hm_stat <- "System Info: Ready"
	reaobj$lp_stat <- "System Info: Ready"
	reaobj$ns_stat <- "System Info: Ready"
	reaobj$snp_ui_on <- F
	reaobj$hp_ui_on <- F
	reaobj$hp_err_on <- F
	reaobj$hn_ui_on <- F
	reaobj$hn_err_on <- F
	reaobj$dt_ui_on <- F
	reaobj$dt_err_on <- F
	reaobj$hm_ui_on <- F
	reaobj$hm_err_on <- F
	reaobj$lp_ui_on <- F
	reaobj$lp_err_on <- F
	reaobj$int_hp_plot_width <- 1000
	reaobj$int_hp_plot_height <- 500

	reaobj$fra_snp_res <- data.frame()
	reaobj$plot_hp_res <- "ggplot"
	reaobj$plot_hm_res <- "ggplot"
	reaobj$text_hm_warnt <- ""
	reaobj$fra_ns_res <- data.frame()

	reaobj$text_snp_para <- "Parameter: "
	reaobj$text_hp_para <- "Parameter: "
	reaobj$text_hn_para <- "Parameter: "
	reaobj$text_hm_para <- "Parameter: "
	reaobj$text_dt_para <- "Parameter: "
	reaobj$text_lp_para <- "Parameter: "
	reaobj$text_ns_para <- "Parameter: "

	output$sam_table <- DT::renderDataTable(
		fra_glo_samshow, filter = 'top', server = TRUE, rownames = FALSE,
		options = list(autoWidth = TRUE, pageLength = 10)
	)

	output$gro_table <- DT::renderDataTable(
		fra_glo_groshow, filter = 'top', server = TRUE, rownames = FALSE,
		options = list(autoWidth = TRUE, pageLength = 10)
	)

	output$chr_table <- DT::renderDataTable(
		fra_glo_chrshow, filter = 'top', server = TRUE, rownames = FALSE,
		options = list(autoWidth = TRUE, pageLength = 10)
	)

	output$sys_table <- DT::renderDataTable(
		fra_glo_sysdata, filter = 'top', server = TRUE, rownames = FALSE,
		options = list(autoWidth = TRUE, pageLength = 10)
	)

	output$snp_ui_1 <- renderUI({
		if(reaobj$snp_ui_on){
			hr()
		}
	})

	output$snp_ui_2 <- renderUI({
		if(reaobj$snp_ui_on){
			textInput("snp_maf", "Minimum allele frequency (MAF)", "0")
		}
	})

	output$snp_ui_3 <- renderUI({
		if(reaobj$snp_ui_on){
			radioButtons("snp_bso", "Biallelic sites only",
                                choices = list("Yes","No"),
                                selected = "No",
                                inline = T)
		}
	})

	output$snp_ui_4 <- renderUI({
		if(reaobj$snp_ui_on){
			textInput("snp_mlr", "Maximum missing rate", "1")
		}
	})

	output$snp_status <- renderText({
		reaobj$snp_stat
	})

	output$snp_restable <- renderTable({
        if(input$snp_run){
        	isolate({
        		reaobj$fra_snp_res <- snp_main(input$snp_ty, input$snp_oi, input$snp_co_t, input$snp_co_f, input$snp_co_e, input$snp_ro, input$snp_ro_ext, input$snp_maf, input$snp_bso, input$snp_mlr)
        		reaobj$text_snp_para <- text_snp_currpara
        		if(names(reaobj$fra_snp_res) == c("Info")){
        			reaobj$snp_stat <- as.character(reaobj$fra_snp_res[1,1])
        			reaobj$fra_snp_res <- data.frame()
        			reaobj$text_snp_para <- "Parameter: "
        		}else{
        			reaobj$snp_stat <- "System Info: Done"
        		}
        		reaobj$fra_snp_res
        	})
        }else{
            data.frame()
        }
    })

	output$snp_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_VarTable", d, ".csv", sep = "")
       },
       content = function(file){
            write.csv(isolate({reaobj$fra_snp_res}), file, row.names = F)
       }
    )

    output$snp_down_tran <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_VarTable", d, ".nucleotide.csv", sep = "")
       },
       content = function(file){
       		tmp_fra <- reaobj$fra_snp_res
       		for(i in 1:nrow(tmp_fra)){
       			ref <- tmp_fra[i, 3]
       			alt <- strsplit(tmp_fra[i, 4], split = ",")[[1]]
       			tran_lis <- c(ref, alt)
       			for(j in 5:ncol(tmp_fra)){
       				ori_s <- tmp_fra[i,j]
       				ori_s <- strsplit(ori_s, split = ":")[[1]]
       				s <- strsplit(ori_s[1], split = "/")[[1]]
       				s1 <- as.numeric(s[1])
       				if(is.na(s1)){s1 <- "N"}else{s1 <- tran_lis[s1+1]}
       				s2 <- as.numeric(s[2])
       				if(is.na(s2)){s2 <- "N"}else{s2 <- tran_lis[s2+1]}
       				s <- paste(s1, s2, sep = "/")
       				ori_s[1] <- s
       				tmp_fra[i,j] <- paste(ori_s, collapse = ":")
       			}
       		}
            write.csv(tmp_fra, file, row.names = F)
       }
    )

	output$hp_plot <- renderUI({
	    
	    output$plot1 <- renderPlot({
			isolate({reaobj$plot_hp_res})
	    })
	    plotOutput("plot1", height = reaobj$int_hp_plot_height+reaobj$int_hp_plot_width*0, brush = brushOpts("snp_brush", delay = 500, delayType ="debounce", resetOnNew = T))
	})

	output$hp_ui_0 <- renderUI({
		if(reaobj$hp_ui_on){
			hr()
		}
	})

	output$hp_ui_1 <- renderUI({
		if(reaobj$hp_ui_on){
			radioButtons("hp_dty", "Figure fomat",
                                choices = list("pdf","png"),
                                selected = "pdf",
                                inline = T)
		}
	})

	output$hp_ui_2 <- renderUI({
		if(reaobj$hp_ui_on){
			if(input$hp_dty == "pdf"){
				textInput("hp_dwi", "Figure width (inch)", "5")
			}else{
				textInput("hp_dwi", "Figure width (100px)", "5")
			}
		}
	})

	output$hp_ui_3 <- renderUI({
		if(reaobj$hp_ui_on){	
			if(input$hp_dty == "pdf"){
				textInput("hp_dth", "Figure height (inch)", "5")
			}else{
				textInput("hp_dth", "Figure height (100px)", "5")
			}
		}
	})

	output$hp_ui_4 <- renderUI({
		if(reaobj$hp_ui_on){
			downloadButton("hp_down", "Download")
		}
	})

	output$hp_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("Snphub_Heatmap", d, ".", input$hp_dty, sep = "")
       },
       content = function(file){
            ggsave(reaobj$plot_hp_res, filename = file, height = as.numeric(input$hp_dth), width = as.numeric(input$hp_dwi), limitsize = FALSE)
       }
    )

	output$hp_status <- renderText({
		reaobj$hp_stat
	})

	output$hn_status <- renderText({
		reaobj$hn_stat
	})

	output$hn_ui_0 <- renderUI({
		if(reaobj$hn_ui_on){
			hr()
		}
	})

	output$hn_ui_1 <- renderUI({
		if(reaobj$hn_ui_on){
			radioButtons("hn_dty", "Figure fomat",
                                choices = list("pdf","png"),
                                selected = "pdf",
                                inline = T)
		}
	})

	output$hn_ui_2 <- renderUI({
		if(reaobj$hn_ui_on){
			if(input$hn_dty == "pdf"){
				textInput("hn_dwi", "Figure width (inch)", "5")
			}else{
				textInput("hn_dwi", "Figure width (100px)", "5")
			}
		}
	})

	output$hn_ui_3 <- renderUI({
		if(reaobj$hn_ui_on){
			if(input$hn_dty == "pdf"){
				textInput("hn_dth", "Figure height (inch)", "5")
			}else{
				textInput("hn_dth", "Figure height (100px)", "5")
			}
		}
	})

	output$hn_ui_4 <- renderUI({
		if(reaobj$hn_ui_on){
			downloadButton("hn_down", "Download")
		}
	})


	output$hn_plot <- renderPlot({
		if(input$hn_run){
			reaobj$hn_err_on <- F
			isolate({plot_hn_res <- hn_main(input$hn_co, input$hn_ro, input$hn_ro_ext)})
			reaobj$text_hn_para <- text_hn_currpara
			if(length(plot_hn_res) != 1){
		    	reaobj$hn_stat <- "System Info: Done"
		    }else{
		    	reaobj$text_hn_para <- "Parameter: "
		    	code <- as.numeric(plot_hn_res)
		    	plot_hn_res <- NA
		    	if(code == 1){
			        reaobj$hn_stat <- "Error 0001:Flanking length must be integer."
			    }else if(code == 2){
			        reaobj$hn_stat <- "Error 0002:Invalid region."
			    }else if(code == 3){
			        reaobj$hn_stat <- "Error 0003:Region is too long."
			    }else if(code == 4){
			        reaobj$hn_stat <- "Error 0004:Invalid accession detected."
			    }else if(code == 101){
			        reaobj$hn_stat <- "Error 0101:No variation found in current regions."
			    }else if(code == 6){
			        reaobj$hn_stat <- "Error 0006:Group not found."
			    }
			    reaobj$hn_err_on <- T
		    }
		    plot_hn_res
        }
	})

	output$hn_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_HapNet", d, ".", input$hn_dty, sep = "")
       },
       content = function(file){
       		if(input$hn_dty == "pdf"){
       			pdf(file, width = as.numeric(input$hn_dwi), height = as.numeric(input$hn_dth))
       		}else if(input$hn_dty == "png"){
       			png(file, width = as.numeric(input$hn_dwi)*100, height = as.numeric(input$hn_dth)*100)
       		}
       		
       		isolate({hn_main(input$hn_co, input$hn_ro, input$hn_ro_ext)})
       		dev.off()
       }
    )

	output$dt_plot <- renderPlot({
		if(input$dt_run){
			reaobj$dt_err_on <- F
			bool_dt_warning <<- F
			isolate({plot_dt_res <- dt_main(input$dt_ty, input$dt_tty, input$dt_dbr, input$dt_co, input$dt_ro, input$dt_ro_ext)})
	    	reaobj$text_dt_para <- text_dt_currpara
	    	if(length(plot_dt_res) != 1){
	    		if(bool_dt_warning){
	    				reaobj$dt_stat <- "Warning:Has NA! Using njs."
	    			}else{
	    				reaobj$dt_stat <- "System Info: Done"
	    			}
		    }else{
		    	reaobj$text_dt_para <- "Parameter: "
		    	code <- as.numeric(plot_dt_res)
		    	#reaobj$plot_dt_res <- NA
		    	if(is.na(code)){
		    		reaobj$dt_stat <- plot_dt_res
		    		reaobj$dt_err_on <- T
		    	}else{
		    		if(code == 1){
				        reaobj$dt_stat <- "Error 0001:Flanking length must be integer."
				    }else if(code == 2){
				        reaobj$dt_stat <- "Error 0002:Invalid region."
				    }else if(code == 3){
				        reaobj$dt_stat <- "Error 0003:Region is too long."
				    }else if(code == 4){
				        reaobj$dt_stat <- "Error 0004:Invalid accession detected."
				    }else if(code == 101){
				        reaobj$dt_stat <- "Error 0101:No variation found in current regions."
				    }else if(code == 201){
				        reaobj$dt_stat <- "Error 0201:Missing genotype presents in given region."
				    }else if(code == 6){
				        reaobj$dt_stat <- "Error 0006:Group not found."
				    }
				    reaobj$dt_err_on <- T
		    	}
		    }
		    plot_dt_res
		}
	})

	output$dt_status <- renderText({
		reaobj$dt_stat
	})

	output$dt_ui_ty_tree <- renderUI({
		if(input$dt_ty == "NJ-tree"){
			selectInput("dt_tty", "Tree layout", 
                                choices = c("phylogram", "cladogram", "fan", "unrooted", "radial"),
                                selected = "unrooted")
		}
	})

	output$dt_ui_ty_tree_byreal <- renderUI({
		if(input$dt_ty == "NJ-tree"){
			radioButtons("dt_dbr", "Real branch length",
                                choices = list("Yes","No"),
                                selected = "Yes",
                                inline = T)
		}
	})

	output$dt_ui_0 <- renderUI({
		if(reaobj$dt_ui_on){
			hr()
		}
	})

	output$dt_ui_1 <- renderUI({
		if(reaobj$dt_ui_on){
			radioButtons("dt_dty", "Figure fomat",
                                choices = list("pdf","png"),
                                selected = "pdf",
                                inline = T)
		}
	})

	output$dt_ui_2 <- renderUI({
		if(reaobj$dt_ui_on){
			if(input$dt_dty == "pdf"){
				textInput("dt_dwi", "Figure width (inch)", "5")
			}else{
				textInput("dt_dwi", "Figure width (100px)", "5")
			}
		}
	})

	output$dt_ui_3 <- renderUI({
		if(reaobj$dt_ui_on){
			if(input$dt_dty == "pdf"){
				textInput("dt_dth", "Figure height (inch)", "5")
			}else{
				textInput("dt_dth", "Figure height (100px)", "5")
			}
		}
	})

	output$dt_ui_4 <- renderUI({
		if(reaobj$dt_ui_on){
			downloadButton("dt_down", "Download")
		}
	})

	output$dt_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_PhyloTree", d, ".", input$dt_dty, sep = "")
       },
       content = function(file){
			if(input$dt_dty == "pdf"){
       			pdf(file, width = as.numeric(input$dt_dwi), height = as.numeric(input$dt_dth))
       		}else if(input$dt_dty == "png"){
       			png(file, width = as.numeric(input$dt_dwi)*100, height = as.numeric(input$dt_dth)*100)
       		}
       		
       		isolate({dt_main(input$dt_ty, input$dt_tty, input$dt_dbr, input$dt_co, input$dt_ro, input$dt_ro_ext)})
       		dev.off()
       	}
    )

    output$hm_plot <- renderPlot({
		reaobj$plot_hm_res
	})

	output$hm_status <- renderText({
		reaobj$hm_stat
	})

	output$hm_ui_0 <- renderUI({
		if(reaobj$hm_ui_on){
			hr()
		}
	})

	output$hm_ui_1 <- renderUI({
		if(reaobj$hm_ui_on){
			radioButtons("hm_dty", "Figure fomat",
                                choices = list("pdf","png"),
                                selected = "pdf",
                                inline = T)
		}
	})

	output$hm_ui_2 <- renderUI({
		if(reaobj$hm_ui_on){
			if(input$hm_dty == "pdf"){
				textInput("hm_dwi", "Figure width (inch)", "5")
			}else{
				textInput("hm_dwi", "Figure width (100px)", "5")
			}
		}
	})

	output$hm_ui_3 <- renderUI({
		if(reaobj$hm_ui_on){
			if(input$hm_dty == "pdf"){
				textInput("hm_dth", "Figure height (inch)", "5")
			}else{
				textInput("hm_dth", "Figure height (100px)", "5")
			}
		}
	})

	output$hm_ui_4 <- renderUI({
		if(reaobj$hm_ui_on){
			downloadButton("hm_down", "Download")
		}
	})

	output$hm_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_HapMap", d, ".", input$hm_dty, sep = "")
       },
       content = function(file){
            ggsave(reaobj$plot_hm_res, filename = file, height = as.numeric(input$hm_dth), width = as.numeric(input$hm_dwi), limitsize = FALSE)
       }
    )

    output$hm_warntext <- renderText({
    	reaobj$text_hm_warnt
    })

    output$lp_status <- renderText({
		reaobj$lp_stat
	})

	output$lp_ui_0 <- renderUI({
		if(reaobj$lp_ui_on){
			hr()
		}
	})

	output$lp_ui_1 <- renderUI({
		if(reaobj$lp_ui_on){
			radioButtons("lp_dty", "Figure fomat",
                                choices = list("pdf","png"),
                                selected = "pdf",
                                inline = T)
		}
	})

	output$lp_ui_2 <- renderUI({
		if(reaobj$lp_ui_on){
			if(input$lp_dty == "pdf"){
				textInput("lp_dwi", "Figure width (inch)", "5")
			}else{
				textInput("lp_dwi", "Figure width (100px)", "5")
			}
		}
	})

	output$lp_ui_3 <- renderUI({
		if(reaobj$lp_ui_on){
			if(input$lp_dty == "pdf"){
				textInput("lp_dth", "Figure height (inch)", "5")
			}else{
				textInput("lp_dth", "Figure height (100px)", "5")
			}
		}
	})

	output$lp_ui_4 <- renderUI({
		if(reaobj$lp_ui_on){
			downloadButton("lp_down", "Download")
		}
	})


	output$lp_plot <- renderPlot({
		if(input$lp_run){
			reaobj$lp_err_on <- F
			isolate({plot_lp_res <- lp_main(input$lp_co, input$lp_ro, input$lp_mft, input$lp_ro_ext)})
			reaobj$text_lp_para <- text_lp_currpara
			if(length(plot_lp_res) != 1){
		    	reaobj$lp_stat <- "System Info: Done"
		    }else{
		    	reaobj$text_lp_para <- "Parameter: "
		    	code <- as.numeric(plot_lp_res)
		    	plot_lp_res <- NA
		    	if(code == 1){
			        reaobj$lp_stat <- "Error 0001:Flanking length must be integer."
			    }else if(code == 2){
			        reaobj$lp_stat <- "Error 0002:Invalid region."
			    }else if(code == 3){
			        reaobj$lp_stat <- "Error 0003:Region is too long."
			    }else if(code == 4){
			        reaobj$lp_stat <- "Error 0004:Invalid accession detected."
			    }else if(code == 5){
			        reaobj$lp_stat <- "Error 0005:Maxium feature tracks must be integer."
			    }else if(code == 101){
			        reaobj$lp_stat <- "Error 0101:No variation found in current regions."
			    }
			    reaobj$lp_err_on <- T
		    }
		    plot_lp_res
        }
	})

	output$lp_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_SnpFreq", d, ".", input$lp_dty, sep = "")
       },
       content = function(file){
       		if(input$lp_dty == "pdf"){
       			pdf(file, width = as.numeric(input$lp_dwi), height = as.numeric(input$lp_dth))
       		}else if(input$lp_dty == "png"){
       			png(file, width = as.numeric(input$lp_dwi)*100, height = as.numeric(input$lp_dth)*100)
       		}
       		
       		isolate({lp_main(input$lp_co, input$lp_ro, input$lp_mft, input$lp_ro_ext)})
       		dev.off()
       }
    )

    output$ns_status <- renderText({
		reaobj$ns_stat
	})

	output$ns_restable <- renderTable({
        if(input$ns_run){
        	isolate({
        		reaobj$fra_ns_res <- ns_main(input$ns_ty1, input$ns_ty2, input$ns_co, input$ns_ro)
        		if(names(reaobj$fra_ns_res) == c("Info")){
        			reaobj$ns_stat <- as.character(reaobj$fra_ns_res[1,1])
        			reaobj$text_ns_para <- paste("Parameter: ", sep="")
        			reaobj$fra_ns_res <- data.frame()
        		}else{
        			reaobj$ns_stat <- "System Info: Done"
        			reaobj$text_ns_para <- paste("Parameter: ", names(reaobj$fra_ns_res), sep="")
        			names(reaobj$fra_ns_res) <- ""
        		}
        		reaobj$fra_ns_res
        	})
        }else{
            data.frame()
        }
    })

    output$ns_down <- downloadHandler(
       filename = function(){
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            paste("SnpHub_SeqMaker", d, ".fasta", sep = "")
       },
       content = function(file){
            write.table(isolate({reaobj$fra_ns_res}), file, row.names = F, col.names = F, quote = F, sep = "")
       }
    )

    output$hp_text <- renderText({
		if(reaobj$hp_err_on){
			"Error Occured"
		}else{
			""
		}
	})

	output$hn_text <- renderText({
		if(reaobj$hn_err_on){
			"Error Occured"
		}else{
			""
		}
	})

	output$dt_text <- renderText({
		if(reaobj$dt_err_on){
			"Error Occured"
		}else{
			""
		}
	})

	output$hm_text <- renderText({
		if(reaobj$hm_err_on){
			"Error Occured"
		}else{
			""
		}
	})

	output$lp_text <- renderText({
		if(reaobj$lp_err_on){
			"Error Occured"
		}else{
			""
		}
	})

	output$snp_res_para <- renderText({
    	reaobj$text_snp_para
    })

    output$hp_res_para <- renderText({
    	reaobj$text_hp_para
    })

    output$hn_res_para <- renderText({
    	reaobj$text_hn_para
    })

    output$dt_res_para <- renderText({
    	reaobj$text_dt_para
    })

    output$hm_res_para <- renderText({
    	reaobj$text_hm_para
    })

    output$lp_res_para <- renderText({
    	reaobj$text_lp_para
    })

    output$ns_para_title <- renderText({
    	reaobj$text_ns_para
    })

    observeEvent(input$snp_adv, {
    	reaobj$snp_ui_on <- !(reaobj$snp_ui_on)
    })

    observeEvent(input$hp_do, {
    	reaobj$hp_ui_on <- !(reaobj$hp_ui_on)
    })

    observeEvent(input$hn_do, {
    	reaobj$hn_ui_on <- !(reaobj$hn_ui_on)
    })

    observeEvent(input$dt_do, {
    	reaobj$dt_ui_on <- !(reaobj$dt_ui_on)
    })

    observeEvent(input$hm_do, {
    	reaobj$hm_ui_on <- !(reaobj$hm_ui_on)
    })

    observeEvent(input$lp_do, {
    	reaobj$lp_ui_on <- !(reaobj$lp_ui_on)
    })

    observeEvent(input$hp_run, {
    	############################
    	## Haplotype Plot
    	############################
    	reaobj$hp_err_on <- F
    	isolate({reaobj$plot_hp_res <- hp_main(input$hp_co, input$hp_ro, input$hp_ro_ext)})
    	reaobj$text_hp_para <- text_hp_currpara
    	if(length(reaobj$plot_hp_res) != 1){
	    	reaobj$int_hp_plot_width <- nrow(fra_hp_orivcf)*15 + 160
	    	reaobj$int_hp_plot_height <- (ncol(fra_hp_orivcf)-1)*20 + 100
	    	reaobj$hp_stat <- "System Info: Done"
	    }else{
	    	reaobj$text_hp_para <- "Parameter: "
	    	code <- as.numeric(reaobj$plot_hp_res)
	    	reaobj$plot_hp_res <- NA
	    	if(code == 1){
		        reaobj$hp_stat <- "Error 0001:Flanking length must be integer."
		    }else if(code == 2){
		        reaobj$hp_stat <- "Error 0002:Invalid region."
		    }else if(code == 3){
				reaobj$hp_stat <- "Error 0003:Region is too long."
		    }else if(code == 4){
		        reaobj$hp_stat <- "Error 0004:Invalid accession detected."
		    }else if(code == 5){
		        reaobj$hp_stat <- "Error 0005:No variation found in current regions."
		    }else if(code == 6){
		        reaobj$hp_stat <- "Error 0006:Group not found."
    		}
    		reaobj$hp_err_on <- T
    		reaobj$int_hp_plot_width <- 500
    		reaobj$int_hp_plot_height <- 500
	    }
    })

    observeEvent(input$hm_run, {
    	############################
    	## HapMap
    	############################
    	reaobj$hm_err_on <- F
    	isolate({reaobj$plot_hm_res <- hm_main(input$hm_co, input$hm_ro, input$hm_lon, input$hm_lat, input$hm_mer)})
    	reaobj$text_hm_warnt <- text_hm_NAwarning
    	if(length(reaobj$plot_hm_res) != 1){
    		reaobj$text_hm_para <- paste("Parameter: ", "samples ", input$hm_co, " ; site ", input$hm_ro, sep="")
    		if(int_hm_warning == 0){
	    		reaobj$hm_stat <- "System Info: Done"
    		}else if(int_hm_warning == 1){
    			reaobj$hm_stat <- paste("Warning:Only the first site (", text_hm_drawsite, ") is used in the inputed region.", sep="")
    		}
	    }else{
	    	reaobj$text_hm_para <- "Parameter: "
	    	code <- as.numeric(reaobj$plot_hm_res)
	    	if(is.na(code)){
	    		reaobj$hm_stat <- as.character(reaobj$plot_hm_res)
	    		reaobj$hm_err_on <- T
	    	}else{
		    	if(code == 1){
			        reaobj$hm_stat <- "Error 0001:Flanking length must be integer."
			    }else if(code == 2){
			        reaobj$hm_stat <- "Error 0002:Invalid region. Did you noticed that just ONE site is needed here?"
			    }else if(code == 3){
					reaobj$hm_stat <- "Error 0003:Region is too long."
			    }else if(code == 4){
			        reaobj$hm_stat <- "Error 0004:Invalid accession detected."
			    }else if(code == 101){
			        reaobj$hm_stat <- "Error 0101:No variation found in current regions."
	    		}
	    		reaobj$hm_err_on <- T
	    	}
	    	reaobj$plot_hm_res <- NA
	    }
    })
})
