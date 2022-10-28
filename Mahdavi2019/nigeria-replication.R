# Nigerian oil network: 
#   R code for generating tables and graphics for
#   "Scraping public co-occurrences for statistical network analysis 
#    of political elites"
# 
# 
# Created by Paasha Mahdavi
# This version: 9 July 2017
#


#-----------------------------------------------------------------------------#
#    Loading packages and data                                                #
#-----------------------------------------------------------------------------#

rm(list = ls())

# Install required packages
pkg <- c("TeachingDemos","statnet","stargazer","miscTools","stringr","ergm.count","texreg",
         "ggplot2","gridExtra","dplyr")
inst <- pkg %in% installed.packages()
if (length(pkg[!inst]) > 0) install.packages(pkg[!inst])
rm(pkg, inst)

# Load required packages 
library("TeachingDemos")
library("statnet")
library("stargazer")
library("miscTools")
library("stringr")
library("ergm.count")
library("texreg")
library("ggplot2")
library("gridExtra")
library("dplyr")



#-----------------------------------------------------------------------------#
#    Setting up code for generating log file                                  #
#-----------------------------------------------------------------------------#

txtStart(file = "nigeria-replication.log", commands = TRUE, results = TRUE)



#-----------------------------------------------------------------------------#
#   2012 appointments                                                         #
#-----------------------------------------------------------------------------#

# Read in scraped dataset

nigeria <- read.csv("nigeria2012.csv",header=FALSE)


# Need to rearrange since hits are in row after each pair

nigeria$degree <- NA

for(i in seq(2,2756,2)){
  nigeria$degree[i-1] <- as.numeric(as.character(nigeria$V2[i]))
}

nigeria <- nigeria[nigeria$V1!=unique(nigeria$V1)[2],]


# Rbind to make the resulting adjacency matrix symmetric: 
#   1326 (52-choose-2) rows of a-b pairs
#   52 rows of a-a pairs
#   1326 (52-choose-2) rows of b-a pairs

names(nigeria) <-c("lastname1","lastname2","degree")
nigeria.b <- nigeria[1:1326,c("lastname2","lastname1","degree")]
names(nigeria.b) <- c("lastname1","lastname2","degree")
nigeria2 <- data.frame(rbind(nigeria,nigeria.b))
nigeria2 <- nigeria2 %>% arrange(lastname2)
nigeria2 <- nigeria2 %>% arrange(lastname1)


# Making into network matrix (adjacency, aka sociomatrix)

nwk.s <- reshape(nigeria2, v.names="degree", timevar="lastname2", 
                 direction="wide", idvar="lastname1")
nwk.s <- as.matrix(nwk.s[1:nrow(nwk.s),2:ncol(nwk.s)])


# Loading in list to get name labels for rows and columns of adj. matrix

cent <- read.csv("nigerialist2012.csv",header=TRUE)
cent <- cent %>% arrange(labels)
rownames(nwk.s) <- cent$labels
colnames(nwk.s) <- cent$labels

# Get shortened names 

shortnames <- function(x,y=cent$labels){
  str_split(y," ")[[x]][length(str_split(y," ")[[x]])]
}
short <- sapply(X = 1:52,FUN = shortnames)


# Making a valued-edge network using statnet commands
#=======================================================================

array_2012 = rep(0, ncol(nwk.s))

for (i in 1 : ncol(nwk.s)){
  for(j in 1 : nrow(nwk.s)){
    if(nwk.s[i,j] > 0){
      print("yes")
      array_2012[i] = array_2012[i] +1
    }
  }
}



hist(array_2012, xlab ='number of connected edges', main = 'Histogram of 2012 data')


array_2015 = rep(0, ncol(nwk2015.s))

for (i in 1 : ncol(nwk2015.s)){
  for(j in 1 : nrow(nwk2015.s)){
    if(nwk2015[i,j] > 0){
      print("yes")
      array_2015[i] = array_2015[i] +1
    }
  }
}



hist(array_2015, xlab ='number of connected edges', main = 'Histogram of 2015 data')


#=======================================================================







nwk <- network(x = nwk.s, 
               directed=FALSE, 
               matrix.type="a", 
               ignore.eval=FALSE, 
               names.eval="hits")

nigeria3 <- as.matrix(nwk, attrname="hits", matrix.type="edgelist")

nwk1 <- nwk

delete.edges(nwk1, seq_along(nwk1$mel))

nwk1[nigeria3[,1:2], names.eval="hits", add.edges=TRUE] <- nigeria3[,3]



# Plotting network
# First getting network characteristics stored as objects

k <- max(log(as.matrix(nwk1, attrname="hits")+1))
nwk.col <- matrix(gray(1-(log(as.matrix(nwk1, attrname="hits")+1)/k)), 
                  nrow=network.size(nwk1))

nwk.width <- log(as.matrix(nwk1, attrname="hits")+1)
nwk.type <- as.numeric(cent$type)+1
nwk.board <- as.numeric(cent$board)+1
nwk.board[nwk.board==3] <- 0
nwk.size <- cent$selfhits <- log(nigeria[1327:1378,3]+1)



#-----------------------------------------------------------------------------#
#    Figure 1 - left panel                                                    #
#-----------------------------------------------------------------------------#

# Network plot

# Loading coordinates
#leftplotcoord <- read.csv(file="leftplotcoord.csv")
plot(nwk1)

#pdf(file = "Figure1-left.pdf", width = 7, height = 7)
plot(nwk1, edge.col = nwk.col, usecurve = TRUE,
     displayisolates=TRUE,
     coord = leftplotcoord,
     edge.curve = .01, edge.lwd = nwk.width, 
     vertex.col = nwk.board+2, vertex.lty = 0)
segments(x0 = 21.65, y0 = -24, x1 = 21, y1 = -32, lwd = 2, lty = 2, col = 2)
mtext(text = "Jonathan", side = 1, col = 2)
#dev.off()


plot(nwk1, edge.label=round(E(nwk1)$weight, 3))

nwk1_matrix = as.matrix(nwk1)

graph_1 = graph.adjacency(nwk1_matrix)


# Computing centrality measures for later modeling

nwk.s <- log(nwk.s+1) # taking logs first because of high skew

cent$net.bw <- betweenness(graph_1)
cent$net.close <- closeness(graph_1)
cent$net.ev <- evcent(graph_1E)
cent$net.degree <- degree(graph_1)

closeness_g1 = closeness(graph_1)
print(paste("mean closeness of 2012 is", mean(closeness_g1, na.rm = TRUE)))
#calculate reciprocity of BA mode
reciprocity_g1 = reciprocity(graph_1)
print(paste("reciprocity of 2012 is", reciprocity_g1))
#calculate transitivity of BA mode
transitivity_g1 = transitivity(graph_1)
print(paste("transitivity of 2012 is", transitivity_g1))
#calculate degree assortativity of BA mode
assortativity_g1 = assortativity_degree(graph_1)
print(paste("assortativity of 2012 is", assortativity_g1))

View(nwk.s)


# 
# 
# 
# # t-test of diff in means
# 
# t.test(net.degree ~ board, data = cent[cent$type!="President",])
# 
# 
# # Modeling using ERGM
# # first adding covariates into network dataframe
# 
# nwk1 %v% "board" <- cent$board
# nwk1 %v% "selfhits" <- nwk.size
# nwk1 %v% "type" <- nwk.type
# nwk1 %v% "region" <- as.numeric(cent$region)
# nwk1 %v% "presregion" <- cent$presregion
# nwk1 %v% "south" <- cent$south
# 
# 
# # Modeling using ERGM (no demographic controls)
# 
# fit.ergm <- ergm(nwk1 ~ nodecov("board") + edges + 
#                    nodematch("type",diff=TRUE,keep=c(1)))
# 
# 
# # Modeling using ERGM (all controls)
# 
# fit.ergm2 <- ergm(nwk1 ~ nodecov("board") + edges + 
#                     nodematch("region",diff=FALSE) +
#                     nodecov("presregion") + 
#                     nodecov("selfhits") + nodematch("type",diff=TRUE,keep=c(1)))
# 
# 
# # Modeling using ERGM-counts (all controls)
# #
# # Getting initial values, intercept init = logged ratio of hits to num. dyads
# #
# m <- sum(nwk1 %e% "hits")/network.dyadcount(nwk1)
# nwk1.sum.init <- log(m)
# #
# # WARNING: Takes a minute, converges after 7 iterations
# #
# fit.ergmc <- ergm(nwk1 ~ sum + nodecov("board") +
#                   nodecov("selfhits") +
#                   nodematch("type",diff=TRUE,keep=c(1)) +
#                   nodematch("region",diff=FALSE) +
#                   nodecov("presregion"),
#                   response = "hits", reference = ~Poisson,
#                   control=control.ergm(init=c(nwk1.sum.init,0,0,0,0,0),
#                                        MCMLE.maxit = 200,
#                                        seed = 1553659388))
# 
# 
# # Modeling using conventional logistic regression
# 
# fit.2012.logit.1 <- glm(board ~ net.ev + presregion + south + selfhits, 
#                         data = cent, 
#                         family = binomial(link = "logit"))
# 
# fit.2012.logit.2 <- glm(board ~ net.degree + presregion + south + selfhits, 
#                         data = cent, 
#                         family = binomial(link = "logit"))
# 



#-----------------------------------------------------------------------------#
#   2015 appointments                                                         #
#-----------------------------------------------------------------------------#

nigeria2015 <- read.csv("nigeria2015.csv",header=FALSE)
nigeria2015$degree <- NA

for(i in seq(2,2756,2)){
  nigeria2015$degree[i-1] <- as.numeric(as.character(nigeria2015$V2[i]))
}

nigeria2015 <- nigeria2015[nigeria2015$V1!=unique(nigeria2015$V1)[2],]
names(nigeria2015) <-c("lastname1","lastname2","degree")

nigeria2015.b <- nigeria2015[1:1326,c("lastname2","lastname1","degree")]
names(nigeria2015.b) <- c("lastname1","lastname2","degree")

nigeria2015.2 <- data.frame(rbind(nigeria2015,nigeria2015.b))
nigeria2015.2 <- nigeria2015.2 %>% arrange(lastname2)
nigeria2015.2 <- nigeria2015.2 %>% arrange(lastname1)


# Making into network matrix (adjacency or sociomatrix)

nwk2015.s <- reshape(nigeria2015.2, 
                     v.names="degree", 
                     timevar="lastname2", 
                     direction="wide", 
                     idvar="lastname1")

nwk2015.s <- as.matrix(nwk2015.s[1:nrow(nwk2015.s),2:ncol(nwk2015.s)])


# Creating centrality list

cent2015 <- read.csv("nigerialist2015.csv",header=TRUE)
cent2015 <- cent2015 %>% arrange(labels)
rownames(nwk2015.s) <- cent2015$labels
colnames(nwk2015.s) <- cent2015$labels
shortnames <- function(x,y=cent2015$labels){
  str_split(y," ")[[x]][length(str_split(y," ")[[x]])]
}
short <- sapply(X = 1:52,FUN = shortnames)


# Making a valued-edge network using statnet commands

nwk2015 <- network(nwk2015.s, 
                   directed=FALSE, 
                   matrix.type="a", 
                   ignore.eval=FALSE, 
                   names.eval="hits")

nigeria2015.3 <- as.matrix(nwk2015, attrname="hits",matrix.type="edgelist")
nwk2015.1 <- nwk2015
delete.edges(nwk2015.1, seq_along(nwk2015.1$mel))
nwk2015.1[nigeria2015.3[,1:2], names.eval="hits", add.edges=TRUE] <- nigeria2015.3[,3]


# Plotting network


nwk2_matrix = as.matrix(nwk2015)

graph_2 = graph.adjacency(nwk2_matrix)



k <- max(log(as.matrix(nwk2015.1, attrname="hits")+1))
nwk2015.col <- matrix(gray(1-(log(as.matrix(nwk2015.1, attrname="hits")+1)/k)), nrow=network.size(nwk2015.1))
nwk2015.width <- log(as.matrix(nwk2015.1, attrname="hits")+1)
nwk2015.type <- as.numeric(cent2015$type)+1
nwk2015.board <- as.numeric(cent2015$board)+1
nwk2015.size <- cent2015$selfhits <- log(nigeria2015[1327:1378,3]+1)



#-----------------------------------------------------------------------------#
#    Figure 1 - right panel                                                   #
#-----------------------------------------------------------------------------#

# Network Plot

# Loading coordinates
rightplotcoord <- read.csv(file="rightplotcoord.csv")

#pdf(file = "Figure1-right.pdf", width = 7, height = 7)
plot(nwk2015.1, edge.col = nwk2015.col, usecurve = TRUE,
     displayisolates=TRUE, 
     coord = rightplotcoord,
     edge.curve = .01, edge.lwd = nwk2015.width, 
     vertex.col = nwk2015.board+2, vertex.lty = 0)
segments(x0 = 16.3, y0 = -8.8, x1 = 20, y1 = -20, lwd = 2, lty = 2, col = 2)
mtext(text = "Buhari", side = 1, col = 2, at = 20)
#dev.off()


#==================================================================================
# Computing centrality measures

nwk2015.s <- log(nwk2015.s+1) # taking logs first because of high skew
cent2015$net.bw <- betweenness(nwk2015.s, cmode="undirected", ignore.eval=FALSE)
cent2015$net.ev <- evcent(nwk2015.s, ignore.eval=FALSE)
cent2015$net.degree <- degree(nwk2015.s, ignore.eval=FALSE)



cent$net.bw <- betweenness(graph_2)
cent$net.close <- closeness(graph_2)
cent$net.ev <- evcent(graph_2)
cent$net.degree <- degree(graph_2)

closeness_g2 = closeness(graph_2)
print(paste("mean closeness of 2015 is", mean(closeness_g2, na.rm = TRUE)))
#calculate reciprocity of BA mode
reciprocity_g2 = reciprocity(graph_2)
print(paste("reciprocity of 2015 is", reciprocity_g2))
#calculate transitivity of BA mode
transitivity_g2 = transitivity(graph_2)
print(paste("transitivity of 2015 is", transitivity_g2))
#calculate degree assortativity of BA mode
assortativity_g2 = assortativity_degree(graph_2)
print(paste("assortativity of 2015 is", assortativity_g2))




# 
# # t-test of diff in means
# t.test(net.degree ~ board, data = cent2015[cent2015$type!="President",])
# 
# 
# # Modeling using ERGM
# 
# nwk2015.1 %v% "board" <- cent2015$board
# nwk2015.1 %v% "selfhits" <- nwk2015.size
# nwk2015.1 %v% "type" <- nwk2015.type
# 
# 
# # Modeling using ERGM (no demographic controls)
# 
# fit2015.ergm <- ergm(nwk2015.1 ~ nodecov("board") + edges + 
#                    nodematch("type",diff=TRUE,keep=c(1)))
# 
# 
# # Modeling using ERGM (all controls)
# 
# fit2015.ergm2 <- ergm(nwk2015.1 ~ nodecov("board") + edges + 
#                     nodecov("selfhits") + nodematch("type",diff=TRUE,keep=c(1)))
# 
# 
# # Modeling using ERGM-counts (all controls)
# # WARNING: Takes a few minutes, converges after 176 iterations
# 
# m <- sum(nwk2015.1 %e% "hits")/network.dyadcount(nwk2015.1)
# nwk2015.1.sum.init <- log(m)
# 
# fit2015.ergmc <- ergm(nwk2015.1 ~ sum + nodecov("board") + 
#                       nodecov("selfhits") +
#                       nodematch("type",diff=TRUE,keep=c(1)), 
#                       response = "hits", reference = ~Poisson, 
#                       control=control.ergm(init=c(nwk2015.1.sum.init,0,0,0), 
#                                            MCMLE.maxit = 200,
#                                            seed = 1553659388))
# #mcmc.diagnostics(fit2015.ergmc)
# #summary(fit2015.ergmc)
# 
# 
# 
# #-----------------------------------------------------------------------------#
# #    Appendix Table 2                                                         #
# #-----------------------------------------------------------------------------#
# 
# # Combined table (2012 and 2015)
# clabs <- c("Board appointee","Edges (density)","Cabinet homophily",
#            "Region homophily","Pres co-ethnic","Google self-hits",
#            "Sum (density)","Board appointee.c","Google self-hits.c",
#            "Cabinet homophily.c","Region homophily.c","Pres co-ethnic.c")
# 
# cat("### --- Appendix Table 2 --- ###")
# 
# texreg(list(fit.ergm, fit.ergm2, fit.ergmc,
#             fit2015.ergm, fit2015.ergm2, fit2015.ergmc), 
#        custom.coef.names = clabs, 
#        dcolumn = TRUE)
# 
# 
# 
# 
# # Modeling using Logit
# 
# fit.2015.logit.1 <- glm(board ~ net.ev + selfhits, 
#                   data = cent2015, 
#                   family = binomial(link = "logit"))
# 
# fit.2015.logit.2 <- glm(board ~ net.degree + selfhits, 
#                   data = cent2015, 
#                   family = binomial(link = "logit"))
# 
# 
# 
# 
# #-----------------------------------------------------------------------------#
# #    Appendix Table 3                                                         #
# #-----------------------------------------------------------------------------#
# 
# # Full logit table for appendix
# 
# cat("### --- Appendix Table 3 --- ###")
# 
# texreg(list(fit.2012.logit.1,fit.2012.logit.2,
#             fit.2015.logit.1,fit.2015.logit.2), 
#        reorder.coef = c(1,2,6,3,4,5),
#        custom.coef.names = c("Intercept", 
#                              "EV centrality", 
#                              "Pres co-ethnic",
#                              "South dummy (region)",
#                              "Google self-hits", 
#                              "Degree centrality"),
#        dcolumn = TRUE)
# 
# 
# 
# #-----------------------------------------------------------------------------#
# #    Stop generating log file                                                 #
# #-----------------------------------------------------------------------------#
# 
# txtStop()
