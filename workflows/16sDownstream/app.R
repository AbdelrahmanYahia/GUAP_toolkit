library("shiny")
library("shinythemes")
library("shinyjs")
suppressMessages(library("phyloseq", quietly = T))
suppressMessages(library("ggplot2", quietly = T))
suppressMessages(library("microbiome", quietly = T))

options(browser="google-chrome")

alphabetaplot <- function(rare_data, bray_curtis, condition="condition"){
  # plot alpha and beta diversity of rarefied data
  p1 <- plot_richness(rare_data , measures = "Shannon", color = condition)+
    facet_wrap(condition, scales = "free_x", nrow = 1)+
    ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)
  p2 <- plot_richness(rare_data, x=condition, measures=c("Observed", "Shannon", "Chao1")) + geom_boxplot()
  p3 <- plot_ordination(rare_data, bray_curtis, color = condition)+stat_ellipse()+
    ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)
  outer <- list("alpha_all" = p1, "alpha_cond" = p2, "beta_all" = p3)
  return(outer)
}
betaSplot <- function(rare_data, bray_curtis, condition="condition",p_name="p_name" ){
  p1 <- plot_ordination(rare_data, bray_curtis, color = p_name, shape = condition)+
    ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)+
    geom_line(aes(group=p_name))
  return(p1)
}
allgenuspercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Genus") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  outer <- list("plot" = plot.mpg, "legend" = legend)
  return(outer)
}
allorderpercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Order") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  outer <- list("plot" = plot.mpg, "legend" = legend)
  return(outer)
}
allfamilypercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Family") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  outer <- list("plot" = plot.mpg, "legend" = legend)
  return(outer)
}
allphylumpercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Phylum") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='bottom')
  return(plot.mpg)
}
bar_compare <- function(ps.top20, condition="condition", comp_name){
  p1 <- plot_bar(ps.top20, "condition", fill="condition", facet_grid=~Genus)+ 
    geom_bar(aes(color=condition, fill=condition), stat="identity", position="stack") + theme(legend.position="bottom")
  p2 <- plot_bar(ps.top20, "Genus", "Abundance", "Genus", 
                 facet_grid="condition") + 
    geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
  outer <- list("plot1" = p1, "plot2" = p2)
  return(outer)
}

ui = fluidPage(
  theme = shinytheme("slate"),
  titlePanel("16s Visualizer"),
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      helpText("Opens phyloseq R object and 
               renders basic 16s plots."),
      fileInput("ps",
                "Phyloseq object:"
      ),
      hidden(p("First, please upload RDS file",
               id = "p_nofiles",
               style = "font-weight:bold;color:red;")),
      sliderInput("raredata", 
                  label = "Rare depth:",
                  min = 500, max = 100000, value = 25000,step = 10),
      textInput("my_condition", h3("Condition: "), 
                value = "condition"),
      checkboxGroupInput("target", h3("Target Organisms"),
                         choices = list("Lactobacillus" = "Lactobacillus", 
                                        "Streptococcus" = "Streptococcus",
                                        "Staphylococcus" = "Staphylococcus",
                                        "Escherichia-Shigella" = "Escherichia-Shigella",
                                        "Klebsiella" = "Klebsiella",
                                        "Bacillus" = "Bacillus",
                                        "Corynebacterium" = "Corynebacterium",
                                        "anthracis" = "anthracis",
                                        "Clostridium" = "Clostridium",
                                        "Enterococcus" = "Enterococcus"),
                         selected = c("Staphylococcus", "Clostridium")),
      actionButton("do", "Run", class = "btn-lg btn-success")
      ),
    mainPanel(
      hidden(p("Processing...",
               id = "p2_nofiles",
               style = "font-weight:bold;")),
      htmlOutput("alpha_1t"),
      plotOutput("alpha_1"),
      htmlOutput("alpha_2t"),
      plotOutput("alpha_2"),
      htmlOutput("alpha_3t"),
      plotOutput("alpha_3"),
      htmlOutput("alpha_4t"),
      plotOutput("alpha_4"),
      htmlOutput("a1t"),
      plotOutput("a1"),
      # plotOutput("a12"),
      htmlOutput("a2t"),
      plotOutput("a2"),
      # plotOutput("a22"),
      htmlOutput("a3t"),
      plotOutput("a3"),
      # plotOutput("a32"),
      htmlOutput("a4t"),
      plotOutput("a4"),
      htmlOutput("a5t"),
      plotOutput("a5"),
      htmlOutput("alpha_5t"),
      plotOutput("alpha_5")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$do,{
    if (is.null(input$ps)) {
      shinyjs::show("p_nofiles")
      shinyjs::hide("p2_nofiles")

    } else {
    shinyjs::show("p2_nofiles")
    shinyjs::hide("p_nofiles")
    ps <- readRDS(input$ps$datapath)
    raredata <- input$raredata
    my_condition=input$my_condition
    rare_data <- rarefy_even_depth(ps, rngseed =1, sample.size = raredata)
    target <<- strsplit(input$target, split=" ")
    ex3 <- subset_taxa(microbiome::transform(rare_data, "compositional"), Genus %in% target)
    pseq2 <- aggregate_taxa(rare_data, "Genus") 
    pseq.rel <- microbiome::transform(rare_data, "compositional") # raltive counts 
    # bar plot of top species per condition
    top20 <- names(sort(taxa_sums(pseq2), decreasing=TRUE))[1:20]
    ps.top20 <- transform_sample_counts(pseq2, function(OTU) OTU/sum(OTU))
    ps.top20 <- prune_taxa(top20, ps.top20)
    bray_curtis <- ordinate(rare_data, method = "PCoA")
    alpha_1 <<- alphabetaplot(rare_data, bray_curtis)
    gen <<- allgenuspercondition(pseq.rel,condition = my_condition)
    comp2 <<- bar_compare(ex3,my_condition,"target_genus") 
    fam <<- allfamilypercondition(pseq.rel,condition = my_condition)
    comp1 <<- bar_compare(ps.top20,my_condition,"top20")
    ord <<- allorderpercondition(pseq.rel,condition = my_condition)
    
    output$alpha_1 <- renderPlot({
      output$alpha_1t <- renderText({
        paste("<h3>Alpha Diversityof all samples: ", "</h3>")
      })
      alpha_1$alpha_all
    })
    
    output$alpha_2 <- renderPlot({
      output$alpha_2t <- renderText({
        paste("<h3>Alpha Diversity ber Condition: ", "</h3>")
      })
      alpha_1$alpha_cond
    })
    
    output$alpha_3 <- renderPlot({
      output$alpha_3t <- renderText({
        paste("<h3>Beta Diversity PCoA: ", "</h3>")
      })
      alpha_1$beta_all
    })
    
    output$alpha_4 <- renderPlot({
      output$alpha_4t <- renderText({
        paste("<h3>Beta Diversity sample Co: ", "</h3>")
      })
      betaSplot(rare_data,bray_curtis)
    })
    
    output$alpha_5 <- renderPlot({
      output$alpha_5t <- renderText({
        paste("<h3>Bar Plot of selected Genus: ", "</h3>")
      })
      comp2$plot1
    })

    output$a1 <- renderPlot({
      output$a1t <- renderText({
        paste("<h3>Bar Plot of Genus: ", "</h3>")
      })
      gen$plot
    })
    # output$a12 <- renderPlot(gen$legend,res = 150, width = 600, height = 1200)
    output$a2 <- renderPlot({
      output$a2t <- renderText({
        paste("<h3>Bar Plot of Family: ", "</h3>")
      })
      fam$plot
    })
    # output$a22 <- renderPlot(fam$legend)
    output$a3 <- renderPlot({
      output$a3t <- renderText({
        paste("<h3>Bar Plot of Order: ", "</h3>")
      })
      ord$plot
    })
    # output$a32 <- renderPlot(ord$legend)
    output$a4 <- renderPlot({
      output$a4t <- renderText({
        paste("<h3>Bar Plot Phylum: ", "</h3>")
      })
      allphylumpercondition(pseq.rel,condition = my_condition)
    })

    output$a5 <- renderPlot({
      output$a5t <- renderText({
        paste("<h3>Bar Plot of top 20 Genus: ", "</h3>")
      })
      comp1$plot2
    })

    showModal(modalDialog(fade = 50,title = NULL,footer = NULL,size = "s",easyClose = T, "Figures are renderring, Please be patient"))
  }
  shinyjs::hide("p2_nofiles")
})}

shinyApp(ui = ui, server = server)