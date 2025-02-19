library(ggplot2)
library(dplyr)

files <- c(
  "2025-2-17-lyq-plot/20240829_1.fastqfilt.fq.filthuman.fq.sam.depth",
# "2025-2-17-lyq-plot/20240829_2.fastqfilt.fq.filthuman.fq.sam.depth",
  "2025-2-17-lyq-plot/20240829_3.fastqfilt.fq.filthuman.fq.sam.depth",
  "2025-2-17-lyq-plot/20241127_2.fastqfilt.fq.filthuman.fq.sam.depth"
)

colors <- c("20240829_1" = "red",  "20240829_3" = "green", "20241127_2" = "orange")

chr_lengths <- data.frame(
  ID = c("CP017623.1","CP017624.1","CP017625.1","CP017626.1","CP017627.1","CP017628.1","CP017629.1","CP017630.1","CP003200.1","CP003223.1","CP003224.1","CP003225.1","CP003226.1","CP003227.1","CP003228.1"),
  Length = c(3188341,2231883,1799298,1603259,1190845,1033292,949511,2286237,5333942,122799,111195,105974,3751,3353,1308)
)

# 合并所有样本数据
df_list <- list()
sample_max_depths <- data.frame(Sample = character(), MaxDepth = numeric())

for (file in files) {
  data <- read.table(file, header = FALSE, sep = "\t", col.names = c("ID", "Position", "Depth"))
  short_name <- sub("\\.fastqfilt.*", "", basename(file))
  data <- left_join(data, chr_lengths, by = "ID") %>%
    mutate(Relative_Position = Position / Length, Sample = short_name)
  
  df_list[[short_name]] <- data
  
  # 计算每个样本的最大深度
  max_depth <- max(data$Depth, na.rm = TRUE)
  sample_max_depths <- rbind(sample_max_depths, data.frame(Sample = short_name, MaxDepth = max_depth))
}

# 按最大深度降序排列样本
sample_max_depths <- sample_max_depths %>% arrange(desc(MaxDepth))
ordered_samples <- sample_max_depths$Sample

# 绑定所有数据
all_data <- bind_rows(df_list)

# 重新设置 Sample 的因子顺序
all_data$Sample <- factor(all_data$Sample, levels = ordered_samples)

# 计算每个基因组的最大深度值，用于标注
label_data <- all_data %>% group_by(ID) %>% summarize(Depth = max(Depth, na.rm = TRUE), Length = first(Length))

# 绘制图像：每张图X轴标刻度，Y轴自适应，样本按最大深度排序
p <- ggplot(all_data, aes(x = Relative_Position, y = Depth)) +
  geom_line(aes(color = Sample)) +
  facet_grid(rows = vars(ID), cols = vars(Sample), scales = "free") +
  geom_text(data = label_data, aes(x = 0.5, y = Depth * 0.95, label = paste0(ID, " (", Length, ")")), size = 3, color = "black") +
  scale_color_manual(values = colors) +
  labs(x = "Position", y = "Depth") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray80"),
        panel.background = element_rect(fill = "white"),
        strip.text = element_text(size = 9, face = "bold"),
        plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("2025-2-17-lyq-plot/multi_sample_depth_plot_by_sample.png", plot = p, width = 16, height = 20)
