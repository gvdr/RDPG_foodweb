library(ggplot2)
library(gtable)
library(RColorBrewer)

p_in_12 <- ggplot(data=traits_in_df, aes(x=X1, y=X2, label=rownames(traits_in_df)))
p_in_12 <- p_in_12 + geom_point(aes(alpha=.45)) + geom_text(angle=30,size=3) + theme_bw()
p_in_23 <- ggplot(data=traits_in_df, aes(x=X1, y=X3, label=rownames(traits_in_df)))
p_in_23 <- p_in_23 + geom_point(aes(alpha=.45)) + theme_bw()


p1 <- p_in_12   + guides(alpha=FALSE)
p2 <- p_in_23   + guides(alpha=FALSE)

gt1 <- ggplot_gtable(ggplot_build(p1))
gt2 <- ggplot_gtable(ggplot_build(p2))

gt <- gtable(widths = unit(c(6, 1), "null"), height = unit(c(1, 1), "null"))
gt <- gtable_add_grob(gt, gt1[,-5], 1, 1)
gt <- gtable_add_grob(gt, gt2[,-5], 2, 1)

grid.newpage()
grid.draw(gt)