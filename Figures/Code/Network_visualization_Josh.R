7


#library
library(tidyverse)
library(readr)
library(igraph)



##input
dat.file <- "Tables/field_lab_networkstatsTop20.csv"



##read in data
dat <- read.csv(dat.file, header = T, row.names = 1) %>%
  select(VarA, VarB, rho.clr, fdr.clr)

#replace periods in names with underscores so names remane the same throughout

dat$VarA <- str_replace_all(dat$VarA,"\\.","-")
dat$VarB <- str_replace_all(dat$VarB,"\\.","-")

###define signficance factor and only keep 

#significance value
sig.val <- 0.05

  
### remove relationships below sig.val and correlations of a variable with itself (where VarA = VarB)
dat.sig <- dat %>%
  filter(fdr.clr < sig.val) %>%
  filter(!VarA == VarB)

####select for only positive interactions
dat.sig.pos <- dat.sig %>%
  filter(rho.clr > 0)



#define nodes of network
nodes.A <- dat.sig.pos %>% 
  select(VarA) %>%
  unique() %>%
  rename("node" = VarA)

nodes.B <- dat.sig.pos %>% 
  select(VarB) %>%
  unique() %>%
  rename("node" = VarB)


nodes.all1 <- full_join(nodes.A, nodes.B, by=NULL)


#define edges of network
edges <- dat.sig.pos %>%
  select(VarA, VarB, rho.clr)
  


#create network object
net1 <- graph_from_data_frame(d=edges, vertices=nodes.all1, directed = F)

#Make edge withs proportional to strength of rho.clr correlation
E(net1)$width <- (E(net1)$rho.clr^3)*5

#plot network
plot(net1, vertex.size = 4, vertex.label.dist = 1, layout = layout_with_fr(net))

#save plot
# library(ggplot2)
# cairo_pdf("Figures/Preliminary/Draft_MS2/Figure_S_network_all.pdf", height = 10, width = 10)
# plot(net1, vertex.size = 4, vertex.label.dist = 1, layout = layout_with_fr(net))
# dev.off()




#try making pretty with ggraph

cairo_pdf("Figures/Preliminary/Draft_MS2/Figure_S_network_all.pdf", height = 5, width = 5)
ggraph(net1, layout= "fr") + 
  geom_edge_link(color="gray85",  aes(width = width), show.legend = FALSE) + 
  scale_edge_width(range = c(0.2, 2))+
  geom_node_point(color="orange", size=3) + 
  geom_node_text(aes(label = nodes.all1$node), size=2, color="black", repel=T) + 
  theme_void()
dev.off()



###Play around with labels slightly 
microbes <- c("Chaetoceros_neogracilis.2", 
                   "Candidatus.Pelagibacter.6",
                   "Fragilariopsis_sublineata.23",
                   "Polaribacter.sp..L3A8.58",
             "Phaeocystis_antarctica.4")

net.2 <- net

V(net.2)$label <- V(net.2)
V(net.2)$name <- ifelse(V(net.2)$name %in% microbes, V(net.2)$name, NA)

V(net.2)$name

#plot network
plot(net.2, vertex.size = 4,
     vertex.label = V(net.2)$name,
     vertex.label.dist = 1, 
     layout = layout_with_fr(net.2))







#-------------------------------------------------------------------------------------------------------------------
#Redo with field samples only


##input
dat.file <- "Tables/field_networkstatsTop20.csv"



##read in data
dat <- read.csv(dat.file, header = T, row.names = 1) %>%
  select(VarA, VarB, rho.clr, fdr.clr)


#replace periods in names with underscores so names remane the same throughout

dat$VarA <- str_replace_all(dat$VarA,"\\.","-")
dat$VarB <- str_replace_all(dat$VarB,"\\.","-")

###define signficance factor and only keep 

#significance value
sig.val <- 0.05


### remove relationships below sig.val and correlations of a variable with itself (where VarA = VarB)
dat.sig <- dat %>%
  filter(fdr.clr < sig.val) %>%
  filter(!VarA == VarB)

####select for only positive interactions
dat.sig.pos <- dat.sig %>%
  filter(rho.clr > 0)



#define nodes of network, this is point where rhizo is getting renamed from .12 to .14 and .15
nodes.A <- dat.sig.pos %>% 
  select(VarA) %>%
  unique() %>%
  rename("node" = VarA)

nodes.B <- dat.sig.pos %>% 
  select(VarB) %>%
  unique() %>%
  rename("node" = VarB)


nodes.all2 <- full_join(nodes.A, nodes.B)

#define edges of network
edges <- dat.sig.pos %>%
  select(VarA, VarB, rho.clr)



#create network object
net2 <- graph_from_data_frame(d=edges, vertices=nodes.all2, directed = F)

#Make edge withs proportional to strength of rho.clr correlation
E(net2)$width <- (E(net2)$rho.clr^3)*5

#plot network
plot(net2, vertex.size = 4, vertex.label.dist = 1, layout = layout_with_fr(net))

#save plot
# library(ggplot2)
# cairo_pdf("Figures/Preliminary/Draft_MS2/Figure_S_network_field.pdf", height = 10, width = 10)
# plot(net2, vertex.size = 4, vertex.label.dist = 1, layout = layout_with_fr(net), edge.color="gray85",vertex.color="orange",
#      vertex.label=nodes.all2$node, vertex.label.font=1, vertex.label.color="black", vertex.label.cex=1)
# dev.off()


#try making pretty with ggraph
install.packages("ggraph")
library(ggraph)

cairo_pdf("Figures/Preliminary/Draft_MS2/Figure_S_network_field.pdf", height = 5, width = 5)
ggraph(net2, layout= "fr") + 
  geom_edge_link(color="gray85",  aes(width = width), show.legend = FALSE) + 
  scale_edge_width(range = c(0.2, 2))+
  geom_node_point(color="orange", size=3) + 
  geom_node_text(aes(label = nodes.all2$node), size=2, color="black", repel=T) + 
  theme_void()
dev.off()


#try plotting together
all <- ggraph(net1, layout= "fr") + 
  geom_edge_link(color="gray85",  aes(width = width), show.legend = FALSE) + 
  scale_edge_width(range = c(0.2, 2))+
  geom_node_point(color="orange", size=3) + 
  geom_node_text(aes(label = nodes.all1$node), size=2, color="black", repel=T) + 
  theme_void()

field <- ggraph(net2, layout= "fr") + 
  geom_edge_link(color="gray85",  aes(width = width), show.legend = FALSE) + 
  scale_edge_width(range = c(0.2, 2))+
  geom_node_point(color="orange", size=3) + 
  geom_node_text(aes(label = nodes.all2$node), size=2, color="black", repel=T) + 
  theme_void()
field

combo <- plot_grid(all, field, labels=c("C", "D"), rel_widths = c(1, 1),rel_heights =  c(1, 1), ncol=2, align = "hv", axis = "l", scale=0.9)
combo

save_plot("Figures/Preliminary/Draft_MS2/Figure_S_networks.pdf", combo, base_height = 8, base_width = 10, units="in")


