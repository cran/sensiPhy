
### decide which plot should be printed (for plot_samp):
which_plot <- function(param = "estimate", graphs = "all",
                       s1 = s1, s2 = s2, opt = opt, s4 = s4,
                       i1 = i1, i2 = i2, i4 = i4, model = model){
    
### ERROR checking:-------------------------------------------------------------
if(param != "estimate" & param != "intercept")
    stop("param = ", param, " is not valid. Valid options are: `estimate` or",
    "`intecept`")
if(graphs != 1 & graphs != 2 & graphs != 3 & graphs != 4 & graphs != "all")
    stop("graphs = ", graphs, " is not valid. Valid options are: `all`, 1, 2, 3, 4.")
    
    
### if model == "BM" or "trend" | estimate: ----------------------------------------
if (model == "BM" | model == "trend" & param == "estimate" & graphs == "all")
    print(multiplot(s1, "No `optpar` plot available for model == `BM` or `trend`", 
                     s2, s4, cols=2))
if (model == "BM" | model == "trend" & param == "estimate" & graphs == 1)
    print(s1)
if (model == "BM" | model == "trend" & param == "estimate" & graphs == 2)
    print(s2)
if (model == "BM" | model == "trend" & param == "estimate" & graphs == 3)
    stop("No `optpar` plot available for model == `BM` or `trend`")
if (model == "BM" | model == "trend" & param == "estimate" & graphs == 4)
    print(s4)
    
### if model == "BM" or "trend" | INTERCEPT: -----------------------------------
if (model == "BM" | model == "trend" & param == "intercept" & graphs == "all")
    print(multiplot(i1, "No `optpar` plot available for model == `BM`", 
                     i2, i4, cols=2))
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 1)
    print(i1)
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 2)
    print(i2)
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 3)
    stop("No `optpar` plot available for model == `BM`")
if (model == "BM" | model == "trend" & param == "intercept" & graphs == 4)
    print(i4)

### if model != "BM" | ESTIMATE: ----------------------------------------------
if (param == "estimate" & graphs=="all")
    suppressMessages(print(multiplot(s1,s4, s2,opt, cols=2)))
if (param == "estimate" & graphs==1)
    suppressMessages(print(s1))
if (param == "estimate" & graphs==2)
    suppressMessages(print(s2))
if (param == "estimate" & graphs==3)
    suppressMessages(print(s4))
if (param == "estimate" & graphs==4)
    suppressMessages(print(opt))
if (param == "intercept" & graphs=="all")
    suppressMessages(print(multiplot(i1,i4,i2,opt,cols=2)))
if (param == "intercept" & graphs==1)
    suppressMessages(print(i1))
if (param == "intercept" & graphs==2)
    suppressMessages(print(i2))
if (param == "intercept" & graphs==3)
    suppressMessages(print(i4))
if (param == "intercept" & graphs==4)
    suppressMessages(print(opt))
}

### Function to plot multiple ggplo2 graphs:------------------------------------
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
        }
    }
}
