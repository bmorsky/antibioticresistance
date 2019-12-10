load("/Users/brycemorsky/Documents/Notre Dame/ND paper 2/Data/best_protocol.Rda")

library(ggplot2)
library(cowplot) 
theme_set(theme_cowplot())
library(viridis)

q1 <- ggplot(best_protocol, aes(x=On,y=Off,size=Success,colour=Competition)) + geom_point() + scale_color_viridis()
q2 <- ggplot(best_protocol, aes(x=Competition,y=Success,colour=On)) + geom_point()

save_plot(plot_grid(q1,q2,ncol=2),filename="~/Desktop/result.png",base_height = 4,base_width = 8)
