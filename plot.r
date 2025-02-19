library(ggplot2)
library(dplyr)
library(tidyr)

files <- c(
  "2025-2-17-lyq-plot/20240829_1.fastqfilt.fq.filthuman.fq.sam.depth",
  "2025-2-17-lyq-plot/20240829_3.fastqfilt.fq.filthuman.fq.sam.depth",
  "2025-2-17-lyq-plot/20241127_2.fastqfilt.fq.filthuman.fq.sam.depth"
)

colors <- c("20240829_1" = "red", "20240829_3" = "green", "20241127_2" = "orange")

chr_lengths <- data.frame(
  ID = c("CP017623.1","CP017624.1","CP017625.1","CP017626.1","CP017627.1","CP017628.1","CP017629.1","CP017630.1"),
  Length = c(3188341,2231883,1799298,1603259,1190845,1033292,949511,2286237)
)

df_list <- list()
sample_max_depths <- data.frame(Sample = character(), MaxDepth = numeric())

for (file in files) {
  data <- read.table(file, header = FALSE, sep = "\t", col.names = c("ID", "Position", "Depth"))
  short_name <- sub("\\.fastqfilt.*", "", basename(file))
  data <- left_join(data, chr_lengths, by = "ID") %>%
    mutate(Relative_Position = Position / Length, Sample = short_name)
  
  # 计算每个样本的最大深度
  max_depth <- max(data$Depth, na.rm = TRUE)
  sample_max_depths <- rbind(sample_max_depths, data.frame(Sample = short_name, MaxDepth = max_depth))
  
  df_list[[short_name]] <- data
}

sample_max_depths <- sample_max_depths %>% arrange(desc(MaxDepth))
ordered_samples <- sample_max_depths$Sample

# 绑定数据
all_data <- bind_rows(df_list)

# **确保 Length 被正确添加**
all_data <- left_join(all_data, chr_lengths, by = "ID")

# **填补缺失 ID，并保留 Length**
all_data <- all_data %>%
  complete(ID = chr_lengths$ID, Sample = ordered_samples, fill = list(Depth = 0, Relative_Position = 0))

# **确保 Length 仍然存在**
all_data <- left_join(all_data, chr_lengths, by = "ID")

# **删除不存在的染色体**
all_data <- semi_join(all_data, chr_lengths, by = "ID")

# **计算最大深度用于标注**
label_data <- all_data %>%
  group_by(ID) %>%
  summarize(Depth = max(Depth, na.rm = TRUE), Length = unique(Length)[1])

p <- ggplot(all_data, aes(x = Relative_Position, y = Depth)) +
  geom_line(aes(color = Sample)) +
  facet_grid(rows = vars(ID), cols = vars(Sample), scales = "free") +
  geom_text(data = label_data, aes(x = 0.5, y = Depth * 0.95 + 1, label = paste0(ID, " (", Length, ")")), size = 3, color = "black") +
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

