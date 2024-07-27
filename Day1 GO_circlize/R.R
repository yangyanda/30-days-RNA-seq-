rm(list = ls())


####----load R Package----####
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
source("combine_enrichment_function.R")

####----load Data----###
# 筛选上调和下调都有的GO Term, 并且个数是前20个的数据
plot_data <- combine_enrichment_files(files = c("GO_result_Up.csv", "GO_result_Down.csv"),
                                      n = 20)


####----plot----####

# 首先绘制整体布局
plot_data_1 <- plot_data %>% dplyr::select(ID, start, end)

# set color
ONTOLOGY_Count <- table(plot_data$ONTOLOGY) %>% as.numeric()
GO_color <- c(rep("#9e9ac8",ONTOLOGY_Count[1]), rep("#d9f0a3", ONTOLOGY_Count[2]), rep("#7bccc4", ONTOLOGY_Count[3]))

pdf("GO_enrichment_circlize_version.pdf", width = 15, height = 8)
circos.par("start.degree" = 90)
# 第一圈 GO Term
circos.initializeWithIdeogram(plot_data_1, plotType = NULL)
circos.track(
  ylim = c(0, 1), 
  track.height = 0.05,  # 轨道高度
  bg.border = NA,  # 不要边框
  bg.col = GO_color, # 添加颜色
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data("ycenter")  
    xlim = get.cell.meta.data("xcenter")
    sector.name = get.cell.meta.data("sector.index")  # 提取 GO Term 
    circos.axis(h = "top", labels.cex = 0.7, major.tick.percentage = 0.4, labels.niceFacing = FALSE)  # 刻度线
    circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  # 添加 GO Term 
  } )

# 第二圈 画 富集到的基因个数
plot_data_2 <- plot_data %>% dplyr::select(ID, start, total_number)
# label data
plot_data_2_label <- plot_data_2 %>% dplyr::select(total_number)

circos.genomicTrackPlotRegion(
  plot_data_2,
  track.height = 0.08, 
  bg.border = "#000000", 
  stack = TRUE,  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "#fa9fb5", border = NA, ...)  
    ylim = get.cell.meta.data("ycenter")  
    xlim = plot_data_2_label[get.cell.meta.data("sector.index"),1] / 2
    sector.name = plot_data_2_label[get.cell.meta.data("sector.index"),1]
    circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  
  } )

# 第三圈 上调和下调的比例
# Up
plot_data_Up <- plot_data %>% 
  dplyr::select(ID, start, Up_percent, end) %>%
  dplyr::mutate(end2 = Up_percent * end) %>%
  dplyr::select(ID, start, end2) %>%
  dplyr::mutate(type = 1) %>%
  dplyr::rename(end = end2)

plot_data_Down <- plot_data %>% 
  dplyr::select(ID, start, Up_percent, end) %>%
  dplyr::mutate(end2 = Up_percent * end) %>%
  dplyr::select(ID, end2, end) %>%
  dplyr::mutate(type = 2) %>%
  dplyr::rename(start = end2)

plot_data_3 <- rbind(plot_data_Up, plot_data_Down)

plot_data_3_label <- plot_data %>%
  dplyr::select(Up_percent, Down_percent, end, Up_Counts, Down_Counts ) %>%
  dplyr::mutate(Up = Up_percent * end,
                Down = end) %>%
  dplyr::select(Up, Down, Up_Counts, Down_Counts)
 
color_assign <- colorRamp2(breaks = c(1, 2), col = c("#ef6548", "#3690c0"))

circos.genomicTrackPlotRegion(
  plot_data_3, 
  track.height = 0.08, 
  bg.border = NA, 
  stack = TRUE,  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
    ylim = get.cell.meta.data("cell.bottom.radius") - 0.5
    xlim = plot_data_3_label[get.cell.meta.data("sector.index"),1] / 2
    sector.name = plot_data_3_label[get.cell.meta.data("sector.index"),3]
    circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = FALSE) 
    xlim = (plot_data_3_label[get.cell.meta.data("sector.index"),2]+plot_data_3_label[get.cell.meta.data("sector.index"),1]) / 2
    sector.name = plot_data_3_label[get.cell.meta.data("sector.index"),4]
    circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = FALSE)  
  } )

# 第四圈 -log10(pvalue) 柱形图

plot_data_4 <- plot_data %>%
  dplyr::select(ID, start, end, total_log) %>%
  dplyr::mutate(total_log = total_log / max(total_log))

circos.genomicTrack(
  plot_data_4, 
  ylim = c(0, 1), 
  track.height = 0.5, 
  bg.col = "#f0f0f0", 
  bg.border = NA,  
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data("sector.index")  
    circos.genomicRect(region, value, col = "#3690c0", border = NA, ytop.column = 1, ybottom = 0, ...) 
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = "#d9d9d9", lwd = 0.3)  
  } )

circos.clear()


# 添加图例
updown_legend <- Legend(
  labels = c("Up", "Down"), 
  type = "points", pch = NA, background = c("#ef6548", "#3690c0"), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"))

category_legend <- Legend(
  labels = c("BP", "CC", "MF"),
  type = "points", pch = NA, background = c("#9e9ac8", "#d9f0a3", "#7bccc4"), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"))

legend_list <- lgd_list_vertical <- packLegend(updown_legend, category_legend)
pushViewport(viewport(x = 0.5, y = 0.5))
grid.draw(lgd_list_vertical)
upViewport()


dev.off()

####----sessionInfo----####
sessionInfo()
  
