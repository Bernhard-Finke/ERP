list.of.packages <- c("shiny", "dplyr", "tidyr", "ggplot2", "RSQLite", "UpSetR", "devtools", "igraph", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.dev.packages <- c( "synaptome.db", "BiocManager", "WGCNA", "clusterCons", "AnNet")
new.dev.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if ("synaptome.db" %in% new.dev.packages){
    devtools::install_github('lptolik/synaptome.db',ref = "nonBioconductorLFS")
}
if ("BiocManager" %in% new.dev.packages){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
}
if ("org.Hs.eg.db" %in% new.dev.packages){
    BiocManager::install("org.Hs.eg.db")
}
if ("WGCNA" %in% new.dev.packages){
    BiocManager::install("WGCNA")
}
if ("clusterCons" %in% new.dev.packages){
    devtools::install_github("biomedicalinformaticsgroup/clusterCons",ref = "main")
}
library(synaptome.db)
library(WGCNA)
library(org.Hs.eg.db)
library(clusterCons)

if ("AnNet" %in% new.dev.packages){
    devtools::install_github("lptolik/AnNet",ref = "develop")
}


library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RSQLite)
library(UpSetR)
library(AnNet)
library(igraph)
library(ggrepel)

pre_cluster <- read.graph("pre_cluster.gml", format=c("gml"))
post_cluster <- read.graph("post_cluster.gml", format=c("gml"))

gg <- calcCentrality(post_cluster)

wt_post <- makeConsensusMatrix(gg, N=100, mask=10, alg="wt", type=2)
fc_post <- makeConsensusMatrix(gg, N=100, mask=10, alg="fc", type=2)
infomap_post <- makeConsensusMatrix(gg, N=100, mask=10, alg="infomap", type=2)
louvain_post <- makeConsensusMatrix(gg, N=100, mask=10, alg="louvain", type=2)
# next
sgG1_pre <- makeConsensusMatrix(gg, N=100, mask=10, alg="sgG1", type=2)
sgG2_pre <- makeConsensusMatrix(gg, N=100, mask=10, alg="sgG2", type=2)
sgG5_pre <- makeConsensusMatrix(gg, N=100, mask=10, alg="sgG5", type=2)

saveRDS(wt_post, file="wt_post_matrix.rds")


con <- dbConnect(drv=RSQLite::SQLite(), dbname="synaptic.proteome_SR_20210408.db.sqlite")


GeneDFraw <- dbGetQuery(conn=con, statement="SELECT Localisation, BrainRegion, MethodID,
                                          GeneID, Paper, Year FROM FullGeneFullPaperFullRegion")

Methods <- dbGetQuery(conn=con, statement="SELECT ID, Name FROM Method")

GeneDF <- as.data.frame(GeneDFraw) %>% merge(Methods, by.x="MethodID", by.y="ID", all=TRUE) %>% select(-MethodID) %>%
    rename('Method' = 'Name')

dbDisconnect(con)

scale <- function(x, VALUE=NULL){
    
    x = as.numeric(as.vector(x))
    
    xmin <- min(x,na.rm=T)
    xmax <- max(x,na.rm=T)
    
    if( is.null(VALUE) ){
        
        x  <- x-xmin
        x  <- ifelse(!is.na(x), x/(xmax-xmin), NA) 
        
        return(x)
    }
    
    value = as.numeric(as.vector(value)[1])
    value = value-xmin
    value = ifelse(!is.na(value), value/(xmax-xmin), NA) 
    return(value)
}

GenerateSubset <- function(loc, reg, met){
    subset <- data.frame(GeneDF)
    if (loc != "All"){
        subset <- subset %>% filter(Localisation == loc)
    }
    if (reg != "All"){
        subset <- subset %>% filter(BrainRegion == reg)
    }
    if (met != "All"){
        subset <- subset %>% filter(Method == met)
    }
    return(subset)
}

GetNumPapersSubset <- function(loc, reg, met){
    subset <- GenerateSubset(loc, reg, met)
    PaperDF <- subset %>% group_by(GeneID) %>% summarise(NumPapers = length(unique(Paper)))
    return(merge(subset, PaperDF, by.x="GeneID", by.y="GeneID", all=TRUE))
}

SubsetNumPapers <- function(loc, reg, met, numpap){
    subset <- GetNumPapersSubset(loc, reg, met)
    return(subset[subset$NumPapers >= as.numeric(numpap), ])
}

text <- "No proteins in current selection"
EmptyPlot <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = text) +  theme_void()


title1a1b <- "Fig 1a and 1b in original paper: Discovery rate of new proteins across all studies, \n where the number of proteins is plotted against the frequency of identification"
title1c1d <- "Fig 1c and 1d in original paper:  Contribution of each study \n to the total number of genes in each category"
title1e1f <- "Fig 1e and 1f in original paper: Accumulation of new genes (black) \n compared to the total datasets (blue) over years"
title1g <- "Fig 1g in original paper: Non-linear fit predicting the total size of consensus genes"
title1h <- "Fig 1h in original paper: Overlap of three synaptic datasets: presynaptic, \n postsynaptic and synaptosomal. Bars correspond to the number of unique genes in each \n compartment and their intersections"

GenerateFig1a1b <- function(loc, reg, met){
    df <- GetNumPapersSubset(loc, reg, met)
    first_occurrence <- df[match(unique(df$GeneID), df$GeneID),]
    dfOut <- as.data.frame(table(first_occurrence$NumPapers)) 
    if (nrow(dfOut) != 0){
        dfOut <- dfOut %>% rename('NumGenes' = 'Var1')
    }
    return(dfOut)
}

GenerateFig1c1d <- function(loc, reg, met, numpap){
    df <- SubsetNumPapers(loc, reg, met, numpap)
    total = length(unique(df$GeneID))
    df <- df %>% group_by(Paper, Year) %>% summarise(NumGenes = length(unique(GeneID)), .groups="keep") %>% 
        mutate(NumGenes = NumGenes / total)
    df$stat = "Found"
    dfNot = df %>% select("Paper", "Year", "NumGenes")
    dfNot = dfNot %>% mutate(NumGenes = 1 - NumGenes)
    dfNot$stat = "NotFound"
    dffull <- rbind(dfNot, df)
    if (nrow(dffull) != 0){
        dffull$stat <- relevel(as.factor(dffull$stat), "NotFound")    
    }
    return(dffull)
}

GenerateFig1e1f <- function(loc, reg, met, numpap){
    df <- SubsetNumPapers(loc, reg, met, numpap) %>% arrange(Year)
    if (nrow(df) != 0){
        df <- df[match(unique(df$GeneID), df$GeneID),] %>% group_by(Year) %>% summarise(GenesInYear = length(GeneID)) %>%
            complete(Year = 2000:2020, fill = list(GenesInYear = 0))
        df$TotalGenes <- cumsum(df$GenesInYear)
    }
    return(df)
}


GenerateFig1h <- function(reg, met, numpap){
    df <- SubsetNumPapers("All", reg, met, numpap) %>% select(GeneID, Localisation)
    df <- unique(df[,c('GeneID','Localisation')])
    if (nrow(df) != 0){
        setPost <- df %>% filter(Localisation == "Postsynaptic") %>% select(-Localisation)
        setPre <- df %>% filter(Localisation == "Presynaptic") %>% select(-Localisation)
        setSynap <- df %>% filter(Localisation == "Synaptosome") %>% select(-Localisation)
        setAll <- intersect(intersect(setPost, setSynap), setPre)
        setonlyPre <- setdiff(setdiff(setPre, setPost), setSynap)
        setonlyPost <- setdiff(setdiff(setPost, setPre), setSynap)
        setonlySynap <- setdiff(setdiff(setSynap, setPre), setPost)
        setSynapPost <- setdiff(intersect(setPost, setSynap), setPre)
        setPostPre <- setdiff(intersect(setPre, setPost), setSynap)
        setSynapPre <- setdiff(intersect(setSynap, setPre), setPost)
        output <- c("syn" = nrow(setonlySynap), "psd" = nrow(setonlyPost), "pres" = nrow(setonlyPre), "syn&psd" = nrow(setSynapPost), "syn&pres" = nrow(setSynapPre), "psd&pres" = nrow(setPostPre), "psd&pres&syn" = nrow(setAll))
        check <- 0
        for (x in output){
            if (x != 0){
                check <- check + 1
            }
            if (check > 1){
                break
            }
        }
        if (check == 1){
            output = c(1, 1)
        }
    }
    else {
        output = 1
    }
    return(output)
}

GenerateFig2c <- function(loc){
    gg <- calcCentrality(presynaptic)
    conmat<- makeConsensusMatrix(gg, N=100, mask=10, alg="louvain", type=2)
    bm<-getBridgeness(gg,'louvain',conmat)
    return(conmat, bm)
}


selectLocalisation <- function(id){
    # selector for localisation
    selectInput(
        inputId = paste("Localisation", id, sep=""),
        label = "Sub-cellular Localisation",
        choices =  c("Postsynaptic", "Presynaptic", "Synaptosome", "All"),
        selected = "Postsynaptic",
        multiple = FALSE
    )
}

selectRegion <- function(id){
    # selector for brain region
    selectInput(
        inputId = paste("BrainRegion", id, sep=""),
        label = "Brain Region",
        choices = c("All", "Brain", "Forebrain", "Midbrain", "Cerebellum", "Hypothalamus",
                    "Hippocampus", "Striatum", "Cerebral cortex", "Frontal lobe", "Occipital lobe",
                    "Temporal lobe", "Parietal Lobe", "Telencephalon", "Prefrontal cortex", 
                    "Motor cortex", "Visual cortex"
        ), 
        selected = "All",
        multiple = FALSE
    )
}

selectMethodName <- function(id){
    # selector for method
    selectInput(
        inputId = paste("Method", id, sep=""),
        label = "Method",
        choices = c("Target", "Shotgun", "All"
        ), 
        selected = "All",
        multiple = FALSE
    )
}

selectNumPapers <- function(id, num=1){
    # selector for method
    selectInput(
        inputId = paste("NumPapers", id, sep=""),
        label = "Evidence (in >= n papers)",
        choices = c(1, 2, 3
        ), 
        selected = num,
        multiple = FALSE
    )
}




ui <- shinyUI(fluidPage(
    
    titlePanel("Executable Research Paper: A unified resource and configurable model of the synapse proteome and its role in disease"),
    
    fluidRow(
        
        column(5,
            selectInput(
                inputId = "fig1",
                label = "Download dataset as csv",
                choices = c("1a,1b", "1c,1d", "1e,1f", "1g", "1h"),
                selected = "1a,1b",
                multiple = FALSE
            ),
            downloadButton("downloaddata"),
            h5("             ")
        ),
        column(5,
            selectInput(
                inputId = "fig2",
                label = "Download figure as png",
                choices = c("1a,1b", "1c,1d", "1e,1f", "1g"),
                selected = "1a,1b",
                multiple = FALSE
               ),
            downloadButton("downloadfig")
            
        )

    ),
    
    sidebarLayout(
        
        sidebarPanel(
            # selector for localisation
            selectLocalisation(1),
            # selector for brain region
            selectRegion(1),
            # selector for method
            selectMethodName(1)
        ),
        
        mainPanel(
            h4(title1a1b),
            # fig1a1b goes here
            plotOutput("fig1a1b"),
            h5('Selection "Localisation = Postsynaptic, BrainRegion = All, Method = All" corresponds to plot 1a in the original paper. 
               \n Similarly, selection "Localisation = Presynaptic, BrainRegion = All, Method = All" corresponds to plot 1b in the original paper.')
        )
        
    ),
    
    sidebarLayout(
        
        sidebarPanel(
            # selector for localisation
            selectLocalisation(2),
            # selector for brain region
            selectRegion(2),
            # selector for method
            selectMethodName(2),
            # selector for number of papers
            selectNumPapers(2)
        ),
        
        mainPanel(
            h4(title1c1d),
            # fig1c1d goes here
            plotOutput("fig1c1d"),
            h5('Selection "Localisation = Postsynaptic, BrainRegion = All, Method = All, Evidence = 1" corresponds to plot 1c in the original paper.
               \n Similarly, selection "Localisation = Presynaptic, BrainRegion = All, Method = All, Evidence = 1" corresponds to plot 1d in the original paper.')
        )
        
    ),
    
    sidebarLayout(
        
        sidebarPanel(
            # selector for localisation
            selectLocalisation(3),
            # selector for brain region
            selectRegion(3),
            # selector for method
            selectMethodName(3),
            # selector for number of papers
            selectNumPapers(3)
        ),
        
        mainPanel(
            h4(title1e1f),
            # fig1e1f goes here
            plotOutput("fig1e1f"),
            h5('Selection "Localisation = Postsynaptic, BrainRegion = All, Method = All, Evidence = 1" corresponds to plot 1e in the original paper.
               \n Similarly, selection "Localisation = Presynaptic, BrainRegion = All, Method = All, Evidence = 1" corresponds to plot 1f in the original paper.')
        )
        
    ),
    
    sidebarLayout(
        
        sidebarPanel(
            # selector for localisation
            selectLocalisation(4),
            # selector for brain region
            selectRegion(4),
            # selector for method
            selectMethodName(4),
            # selector for number of papers
            selectNumPapers(4, 2),
            # selector for fit type
            selectInput(
                inputId = "fitType",
                label = "Type of Fit",
                choices = c("Logistic", "Linear"
                ), 
                selected = "Logistic",
                multiple = FALSE
            )
        ),
        
        mainPanel(
            h4(title1g),
            # fig1g goes here
            plotOutput("fig1g"),
            h5('Selection "Localisation = Postsynaptic, BrainRegion = All, Method = All, Evidence = 2, Type of Fit = Logistic" corresponds to plot 1g in the original paper.')
        )
        
    ),
    
    sidebarLayout(
        
        sidebarPanel(
            # selector for brain region
            selectRegion(5),
            # selector for method
            selectMethodName(5),
            # selector for number of papers
            selectNumPapers(5)
        ),
        
        mainPanel(
            h4(title1h),
            # fig1h goes here
            plotOutput("fig1h"),
            h5('Selection "BrainRegion = All, Method = All, Evidence = 1" corresponds to plot 1h in the original paper.')
        )
        
    ),
    
    sidebarLayout(
        
        sidebarPanel(
            selectRegion(6)
        ),
        
        mainPanel(
            h4("bleep bloop"),
            plotOutput("fig2c")
        )
    )
))

server <- shinyServer(function(input, output) {
    
    fig1a1b <- reactive({
        data1a1b  <- GenerateFig1a1b(input$Localisation1, input$BrainRegion1, input$Method1)
        if (nrow(data1a1b) != 0) {
            ggplot(data1a1b, aes(x=NumGenes, y=Freq)) + geom_bar(stat="identity") + xlab("Number of Protein Identifications") + 
                ylab("Number of Proteins")
        }
        else {
            EmptyPlot
        }
    })
    
    output$fig1a1b <- renderPlot({fig1a1b()})
    
    fig1c1d <- reactive({
        data1c1d  <- GenerateFig1c1d(input$Localisation2, input$BrainRegion2, input$Method2, input$NumPapers2)
        if (nrow(data1c1d) != 0) {
            ggplot(data1c1d, aes(x=reorder(Paper, Year), y=NumGenes, fill= stat)) + geom_bar(stat="identity")  + ylab("count") +
                theme(axis.text.x=element_text(angle=90)) + xlab("") + scale_fill_manual(values=c("darkorchid4", "yellow"))
        }
        else {
            EmptyPlot
        }
    })
    
    output$fig1c1d <- renderPlot({fig1c1d()})
    
    fig1e1f <- reactive({
        data1e1f  <- GenerateFig1e1f(input$Localisation3, input$BrainRegion3, input$Method3, input$NumPapers3)
        if (nrow(data1e1f) != 0) {
            ggplot(data1e1f, aes(x=Year, y=TotalGenes)) + geom_area(fill="#56B4E9") + geom_area(aes(x=Year, y=GenesInYear)) + 
                xlab("years") + ylab("total")
        }
        else {
            EmptyPlot
        }
    })
    
    output$fig1e1f <- renderPlot({fig1e1f()})
    
    fig1g <- reactive({
        data1g <- GenerateFig1e1f(input$Localisation4, input$BrainRegion4, input$Method4, input$NumPapers4)
        if (nrow(data1g) != 0) {
            if (input$fitType == "Linear"){
                fit <- geom_smooth(method='lm', formula= y~x, fullrange = TRUE, se = FALSE)
                error <- round(summary(lm(Year~TotalGenes, data=data1g))$r.squared, 2)
                ggplot(data1g, aes(x=Year, y=TotalGenes)) + geom_point() + fit + xlab("years") + 
                    ylab("total") + geom_text(label = paste("Standard error:", error), x = 2015, y = 150, size = 5)
            } 
            else if (input$fitType == "Logistic"){
                check <- 0
                for (x in data1g$GenesInYear){
                    if (x != 0){
                        check <- check + 1
                    }
                    if (check > 1){
                        break
                    }
                }
                if (check == 1){
                    text <- "Only one paper in current selection: \n no logistic regression can be fitted"
                    ggplot() + annotate("text", x = 4, y = 25, size=8, label = text) +  theme_void()
                }
                else {
                    fit <- geom_smooth(method='nls', formula= y ~ SSlogis(x, Asym, xmid, scal), fullrange = TRUE, se = FALSE)
                    asymptote <- round(coef(nls(TotalGenes ~ SSlogis(Year, Asym, xmid, scal), data = data1g))[1])
                    ggplot(data1g, aes(x=Year, y=TotalGenes)) + geom_point() + fit + 
                        xlab("years") + ylab("total") + xlim(2000, 2050) + geom_hline(yintercept = asymptote) +
                        geom_text(label = paste("Total predicted \n genes:", asymptote), x = 2040, y = asymptote / 2, size = 5)
                }
            }
            
        }
        else {
            EmptyPlot 
        }
    })
    
    output$fig1g <- renderPlot({fig1g()})
    
    fig1h <- reactive({
        data1h <- GenerateFig1h(input$BrainRegion5, input$Method5, input$NumPapers5)
        if (length(data1h) == 7) {
            upset(fromExpression(data1h), text.scale=1.5, mainbar.y.label = "Number of genes", sets.x.label = "Compartments", 
                  sets.bar.color = c("purple4", "#56B4E9", "orangered3"))
        }
        else if (length(data1h) == 1) {
            EmptyPlot
        }
        else{
            text <- "No overlap: current selection contains \n proteins uniquely found in one compartment"
            ggplot() + annotate("text", x = 4, y = 25, size=8, label = text) +  theme_void()
        }
    })
    
    output$fig1h <- renderPlot({fig1h()})
    
    
    fig2c <- reactive({
        
        bm<-getBridgeness(gg,'louvain',conmat)
        VIPs=c('8495','22999','8927','8573','26059','8497','27445','8499')
        indx <- match(V(gg)$name,VIPs)
        group <- ifelse(is.na(indx), 0,1)
        MainDivSize <- 0.8
        Xlab <- "Semilocal Centrality (SL)"
        Ylab <- "Bridgeness"
        X <- as.numeric(igraph::get.vertex.attribute(gg,"SL",V(gg)))
        X <- scale(X)
        Y <- as.numeric(as.vector(bm[,dim(bm)[2]]))
        lbls <- ifelse(!is.na(indx),V(gg)$GeneName,"")
        dt<-data.frame(X=X,Y=Y,vips=group,entres=V(gg)$name,name=V(gg)$GeneName)
        dt_vips <-dt[dt$vips==1,]
        dt_res <- dt[dt$vips==0,]
        baseColor="royalblue2"

        ggplot(dt,aes(x=X,y=Y,label=name))+#geom_point()+
            geom_point(data=dt_vips,
                       aes(x=X,y=Y),colour=baseColor,size = 7,shape=15,show.legend=F)+
            geom_point(data=dt_res,
                       aes(x=X,y=Y, alpha=(X*Y)), size = 3,shape=16,show.legend=F)+
            geom_label_repel(aes(label=as.vector(lbls)),
                             fontface='bold',color='black',fill='white',box.padding=0.1,
                             point.padding=NA,label.padding=0.15,segment.color='black',
                             force=1,size=rel(3.8),show.legend=F,max.overlaps=200)+
            labs(x=Xlab,y=Ylab,title=sprintf("%s","louvain"))+
            scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
            theme(
                axis.title.x=element_text(face="bold",size=rel(2.5)),
                axis.title.y=element_text(face="bold",size=rel(2.5)),
                legend.title=element_text(face="bold",size=rel(1.5)),
                legend.text=element_text(face="bold",size=rel(1.5)),
                legend.key=element_blank())+
            theme(panel.grid.major = element_line(colour="grey40",size=0.2),
                  panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill="white"),
                  panel.border = element_rect(linetype="solid",fill=NA))+
            geom_vline(xintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
            geom_hline(yintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)
    })
    
    output$fig2c <- renderPlot({fig2c()})
    
    
    datasetInput <- reactive({
        switch(input$fig1,
        "1a,1b" = GenerateFig1a1b(input$Localisation1, input$BrainRegion1, input$Method1),
        "1c,1d" = GenerateFig1c1d(input$Localisation2, input$BrainRegion2, input$Method2, input$NumPapers2),
        "1e,1f" = GenerateFig1e1f(input$Localisation3, input$BrainRegion3, input$Method3, input$NumPapers3),
        "1g" = GenerateFig1e1f(input$Localisation4, input$BrainRegion4, input$Method4, input$NumPapers4),
        "1h" = GenerateFig1h(input$BrainRegion5, input$Method5, input$NumPapers5)
        )
    })
    
    FigInput <- reactive({
        switch(input$fig2,
               "1a,1b" = fig1a1b(),
               "1c,1d" = fig1c1d(),
               "1e,1f" = fig1e1f(),
               "1g" = fig1g()
        )
    })
    
    
    output$downloaddata <- downloadHandler(
        filename = function(){
            paste("data", input$fig1, ".csv", sep="")
        },
        content = function(file){
            if (input$fig1 == "1h"){
                keep.names = TRUE
            }
            else {
                keep.names = FALSE
            }
            write.csv(datasetInput(), file, row.names = keep.names)
        }
    )
    
    output$downloadfig <- downloadHandler(
        filename = function(){
            paste("fig", input$fig2, ".png", sep="")
        },
        content = function(file){
            ggsave(file, plot = FigInput(), device = 'png')
        }
    )
    
    
})

shinyApp(ui = ui, server = server)