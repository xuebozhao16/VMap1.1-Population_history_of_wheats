library(admixturegraph)
#####测试代码1
data(bears)
bears
leaves <- c("BLK", "PB",
            "Bar", "Chi1", "Chi2", "Adm1", "Adm2",
            "Denali", "Kenai", "Sweden") 
inner_nodes <- c("R", "PBBB",
                 "Adm", "Chi", "BC", "ABC",
                 "x", "y", "z",
                 "pb_a1", "pb_a2", "pb_a3", "pb_a4",
                 "bc_a1", "abc_a2", "x_a3", "y_a4")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "pb_a1"),
                        edge("pb_a1", "pb_a2"),
                        edge("pb_a2", "pb_a3"),
                        edge("pb_a3", "pb_a4"),
                        edge("pb_a4", "PBBB"),
                        
                        edge("Chi1", "Chi"),
                        edge("Chi2", "Chi"),
                        edge("Chi", "BC"),
                        edge("Bar", "BC"),
                        edge("BC", "bc_a1"),
                        
                        edge("Adm1", "Adm"),
                        edge("Adm2", "Adm"),
                        
                        admixture_edge("bc_a1", "pb_a1", "ABC", "a"),
                        edge("Adm", "ABC"),
                        
                        edge("ABC", "abc_a2"),
                        admixture_edge("abc_a2", "pb_a2", "x", "b"),
                        
                        edge("Denali", "x"),
                        edge("x", "x_a3"),
                        admixture_edge("x_a3", "pb_a3", "y", "c"),
                        
                        edge("Kenai", "y"),
                        edge("y", "y_a4"),                        
                        admixture_edge("y_a4", "pb_a4", "z", "d"),
                        
                        edge("Sweden", "z"),
                        
                        edge("z", "PBBB"),
                        edge("PBBB", "R")))

bears_graph <- agraph(leaves, inner_nodes, edges)
plot(bears_graph, show_admixture_labels = TRUE)
#####测试代码2
fit <- fit_graph(bears, bears_graph)
fit
summary(fit)
plot(fit)
#####测试代码3
data(bears)
leaves <- c("BLK", "PB", "AK", "ABC_A", "ABC_BC", "YB", "EBB") 
inner_nodes <- c("R", "a", "b", "c", "d", "e", "f", "g", "h",
                 "abc_a", "G", "E")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "h"),
                        edge("AK", "h"),
                        edge("ABC_A", "abc_a"),
                        admixture_edge("abc_a", "f", "g"),
                        edge("ABC_BC", "g"), 
                        edge("YB", "e"),
                        edge("EBB", "c"),
                        edge("h", "f"),
                        edge("f", "d"),
                        edge("g", "G"),
                        admixture_edge("G", "d", "e"),
                        edge("d", "b"),
                        edge("e", "E"),
                        admixture_edge("E", "b", "c"),
                        edge("b", "a"),
                        edge("c", "a"),
                        edge("a", "R")))


admixtures <- admixture_proportions(c(admix_props("abc_a", "f", "g", "alpha"),
                                      admix_props("G", "d", "e", "beta"),
                                      admix_props("E", "b", "c", "gamma")))

bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
plot(bears_graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)


########开始做自己的数据
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/qpgraph/EA_land_f3")
qpAsian = read.delim("qpgraph_out4.txt")
qpAsian$Z.value = abs(qpAsian$Z.value)
leaves <- c("out", "out1",
            "L1", "L7", "L8", "L3", "L6",
            "L5", "L4", "L2") 
inner_nodes <- c("R", "out1BB",
                 "Adm", "Chi", "BC", "ABC",
                 "x", "y", "z",
                 "pb_a1", "pb_a2", "pb_a3", "pb_a4",
                 "bc_a1", "abc_a2", "x_a3", "y_a4")

edges <- parent_edges(c(edge("out", "R"),
                        edge("out1", "pb_a1"),
                        edge("pb_a1", "pb_a2"),
                        edge("pb_a2", "pb_a3"),
                        edge("pb_a3", "pb_a4"),
                        edge("pb_a4", "out1BB"),
                        
                        edge("L7", "Chi"),
                        edge("L8", "Chi"),
                        edge("Chi", "BC"),
                        edge("L1", "BC"),
                        edge("BC", "bc_a1"),
                        
                        edge("L3", "Adm"),
                        edge("L6", "Adm"),
                        
                        admixture_edge("bc_a1", "pb_a1", "ABC", "a"),
                        edge("Adm", "ABC"),
                        
                        edge("ABC", "abc_a2"),
                        admixture_edge("abc_a2", "pb_a2", "x", "b"),
                        
                        edge("L5", "x"),
                        edge("x", "x_a3"),
                        admixture_edge("x_a3", "pb_a3", "y", "c"),
                        
                        edge("L4", "y"),
                        edge("y", "y_a4"),                        
                        admixture_edge("y_a4", "pb_a4", "z", "d"),
                        
                        edge("L2", "z"),
                        
                        edge("z", "out1BB"),
                        edge("out1BB", "R")))

qpAsian_graph <- agraph(leaves, inner_nodes, edges)
plot(qpAsian_graph, show_admixture_labels = TRUE)
#####测试代码2
fit <- fit_graph(qpAsian, qpAsian_graph)
fit
summary(fit)
plot(fit)














