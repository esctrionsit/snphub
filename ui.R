source("./config.R")
source("./scripts/public.R")
source("./scripts/haplo.R")
source("./scripts/hapnet.R")
source("./scripts/distree.R")
source("./scripts/newstring.R")
source("./scripts/select.R")
source("./scripts/hapmap.R")
source("./scripts/lollipop.R")

tagList(
    includeCSS("style.css"),
    navbarPage("SnpHub",
        tabPanel("SampleInfo",
            h2("Avaliable samples"),
            shiny::tags$table (width = "100%",
                shiny::tags$tr(
                    shiny::tags$th(width = "100%", style="word-break:break-all;",
                        DT::dataTableOutput("sam_table")
                    )
                )
            ),

            hr(),

            h2("Avaliable chromosomes"),
            shiny::tags$table (width = "100%",
                shiny::tags$tr(
                    shiny::tags$th(width = "100%", style="word-break:break-all;",
                        DT::dataTableOutput("chr_table")
                    )
                )
            ),

            hr(),

            h2("Pre-defined sample groups"),
            shiny::tags$table (width = "100%",
                shiny::tags$tr(
                    shiny::tags$th(width = "100%", style="word-break:break-all;",
                        DT::dataTableOutput("gro_table")
                    )
                )
            ),
            h2("Other system informations"),
            shiny::tags$table (width = "100%",
                shiny::tags$tr(
                    shiny::tags$th(width = "100%", style="word-break:break-all;",
                        DT::dataTableOutput("sys_table")
                    )
                )
            )
        ),

        tabPanel("VarTable",
            sidebarLayout(
                sidebarPanel(
                    radioButtons("snp_ty", "Variation type",
                        choices = list("snp","indel", "snp+indel"),
                        selected = "snp+indel",
                        inline = T),

                    textInput("snp_co_t", "List of samples, must have variant", json_glo_UIsetting$VarTable$Samples1),

                    textInput("snp_co_f", "List of samples, must NOT have variant", json_glo_UIsetting$VarTable$Samples2),

                    textInput("snp_co_e", "List of samples, independent of having variant", json_glo_UIsetting$VarTable$Samples3),

                    textInput("snp_ro", "Region", json_glo_UIsetting$VarTable$Region),

                    textInput("snp_ro_ext", "Flanking region length", json_glo_UIsetting$VarTable$Flanking),

                    actionButton("snp_adv", "Advanced options"),

                    uiOutput("snp_ui_1"),

                    uiOutput("snp_ui_2"),

                    uiOutput("snp_ui_3"),

                    uiOutput("snp_ui_4"),

                    hr(),

                    actionButton("snp_run", "Run"),

                    downloadButton("snp_down", "Download raw results as CSV"),

                    downloadButton("snp_down_tran", "Download genotype in nucleotide as CSV")
                ),
                mainPanel(
                    verbatimTextOutput("snp_status", placeholder = TRUE),
                    verbatimTextOutput("snp_res_para", placeholder = TRUE),
                    h4("Results"),
                    tableOutput("snp_restable")
                )
            )
        ),

        tabPanel("Heatmap",
            sidebarLayout(
                sidebarPanel(
                    textInput("hp_co", "Groups", json_glo_UIsetting$Heatmap$Groups),

                    helpText("Ex: Group1{Sample1,Sample2},Group2{Sample3,Sample4}"),

                    textInput("hp_ro", "Region", json_glo_UIsetting$Heatmap$Region),

                    textInput("hp_ro_ext", "Flanking region length", json_glo_UIsetting$Heatmap$Flanking),

                    actionButton("hp_run", "Draw"),

                    actionButton("hp_do", "Download options"),

                    uiOutput("hp_ui_0"),

                    uiOutput("hp_ui_1"),

                    uiOutput("hp_ui_2"),

                    uiOutput("hp_ui_3"),

                    uiOutput("hp_ui_4")
                ),
                mainPanel(
                    verbatimTextOutput("hp_status", placeholder = TRUE),
                    verbatimTextOutput("hp_res_para", placeholder = TRUE),
                    h3(textOutput("hp_text")),
                    uiOutput("hp_plot")
                )
            )
        ),

        tabPanel("HapNet",
            sidebarLayout(
                sidebarPanel(
                    textInput("hn_co", "Groups", json_glo_UIsetting$HapNet$Groups),

                    helpText("Ex: Group1{Sample1,Sample2},Group2{Sample3,Sample4}"),

                    textInput("hn_ro", "Region", json_glo_UIsetting$HapNet$Region),

                    textInput("hn_ro_ext", "Flanking region length", json_glo_UIsetting$HapNet$Flanking),

                    actionButton("hn_run", "Draw"),

                    actionButton("hn_do", "Download options"),

                    uiOutput("hn_ui_0"),

                    uiOutput("hn_ui_1"),

                    uiOutput("hn_ui_2"),

                    uiOutput("hn_ui_3"),

                    uiOutput("hn_ui_4")
                ),
                mainPanel(
                    verbatimTextOutput("hn_status", placeholder = TRUE),
                    verbatimTextOutput("hn_res_para", placeholder = TRUE),
                    h3(textOutput("hn_text")),
                    plotOutput("hn_plot")
                )
            )
        ),

        tabPanel("PhyloTree",
            sidebarLayout(
                sidebarPanel(
                    radioButtons("dt_ty", "Type",
                                    choices = list("NJ-tree","MDS"),
                                    selected = "NJ-tree",
                                    inline = T),

                    uiOutput("dt_ui_ty_tree"),

                    uiOutput("dt_ui_ty_tree_byreal"),

                    textInput("dt_co", "Groups", json_glo_UIsetting$PhyloTree$Groups),

                    helpText("Ex: Group1{Sample1,Sample2},Group2{Sample3,Sample4}"),

                    textInput("dt_ro", "Region", json_glo_UIsetting$PhyloTree$Region),

                    textInput("dt_ro_ext", "Flanking region length", json_glo_UIsetting$PhyloTree$Flanking),

                    actionButton("dt_run", "Draw"),

                    actionButton("dt_do", "Download options"),

                    uiOutput("dt_ui_0"),

                    uiOutput("dt_ui_1"),

                    uiOutput("dt_ui_2"),

                    uiOutput("dt_ui_3"),

                    uiOutput("dt_ui_4")
                ),
                mainPanel(
                    verbatimTextOutput("dt_status", placeholder = TRUE),
                    verbatimTextOutput("dt_res_para", placeholder = TRUE),
                    h3(textOutput("dt_text")),
                    plotOutput("dt_plot", height = "500px")
                )
            )
        ),

        tabPanel("HapMap",
            sidebarLayout(
                sidebarPanel(
                    textInput("hm_co", "Samples", json_glo_UIsetting$HapMap$Samples),

                    textInput("hm_ro", "Site", json_glo_UIsetting$HapMap$Region),

                    helpText("Input should be a single site (Ex: chr1A:220037)"),

                    sliderInput(inputId = "hm_lon",
                                        label = "Longitude range",
                                        min = -155,
                                        max = 140,
                                        value = c(-155,140)
                                    ),

                    sliderInput(inputId = "hm_lat",
                                        label = "Latitude range",
                                        min = -50,
                                        max = 71,
                                        value = c(-50, 71)
                                    ),

                    sliderInput(inputId = "hm_mer",
                                        label = "Diatance merge range",
                                        min = 0.1,
                                        max = 10,
                                        value = 1
                                    ),

                    actionButton("hm_run", "Draw"),

                    actionButton("hm_do", "Download Options"),

                    uiOutput("hm_ui_0"),

                    uiOutput("hm_ui_1"),

                    uiOutput("hm_ui_2"),

                    uiOutput("hm_ui_3"),

                    uiOutput("hm_ui_4")
                ),
                mainPanel(
                    verbatimTextOutput("hm_status", placeholder = TRUE),
                    verbatimTextOutput("hm_res_para", placeholder = TRUE),
                    h3(textOutput("hm_text")),
                    plotOutput("hm_plot"),
                    textOutput("hm_warntext")
                )
            )
        ),

        tabPanel("SnpFreq",
            sidebarLayout(
                sidebarPanel(
                    textInput("lp_co", "Samples", json_glo_UIsetting$SnpFreq$Samples),

                    textInput("lp_ro", "Region", json_glo_UIsetting$SnpFreq$Region),

                    textInput("lp_ro_ext", "Flanking region length", json_glo_UIsetting$SnpFreq$Flanking),

                    textInput("lp_mft", "Maxium feature tracks", json_glo_UIsetting$SnpFreq$tracks),

                    actionButton("lp_run", "Draw"),

                    actionButton("lp_do", "Download Options"),

                    uiOutput("lp_ui_0"),

                    uiOutput("lp_ui_1"),

                    uiOutput("lp_ui_2"),

                    uiOutput("lp_ui_3"),

                    uiOutput("lp_ui_4")
                ),
                mainPanel(
                    verbatimTextOutput("lp_status", placeholder = TRUE),
                    verbatimTextOutput("lp_res_para", placeholder = TRUE),
                    h3(textOutput("lp_text")),
                    plotOutput("lp_plot", height = "600px")
                )
            )
        ),

        tabPanel("SeqMaker",
            sidebarLayout(
                sidebarPanel(
                    radioButtons("ns_ty1", "Variation type to be replaced",
                                choices = list("snp","indel", "both"),
                                selected = "both",
                                inline = T),

                    radioButtons("ns_ty2", "Replace homozygous/heterozygous?",
                                choices = list("homo_only", "all"),
                                selected = "all",
                                inline = T),
                    
                    textInput("ns_co", "Samples", json_glo_UIsetting$SeqMaker$Samples),

                    textInput("ns_ro", "Region", json_glo_UIsetting$SeqMaker$Region),
                    
                    actionButton("ns_run", "Run"),

                    downloadButton("ns_down", "Download as fasta")
                ),
                mainPanel(
                    verbatimTextOutput("ns_status", placeholder = TRUE),
                    verbatimTextOutput("ns_para_title", placeholder = TRUE),
                    h4("Consensus sequence(s)"),
                    tableOutput("ns_restable")
                )
            )
        )
    )
)