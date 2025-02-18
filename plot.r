library(ggplot2)
library(dplyr)

files <- c("2025-2-17-lyq-plot/20240829_1.fastqfilt.fq.filthuman.fq.sam.depth",
           "2025-2-17-lyq-plot/20240829_2.fastqfilt.fq.filthuman.fq.sam.depth",
           "2025-2-17-lyq-plot/20240829_3.fastqfilt.fq.filthuman.fq.sam.depth",
           "2025-2-17-lyq-plot/P-S4.fastqfilt.fq.filthuman.fq.sam.depth",
           "2025-2-17-lyq-plot/20241127_2.fastqfilt.fq.filthuman.fq.sam.depth",
           "2025-2-17-lyq-plot/P-S3.fastqfilt.fq.filthuman.fq.sam.depth",
           "2025-2-17-lyq-plot/P-S11.fastqfilt.fq.filthuman.fq.sam.depth")

# 自定义颜色 - Corrected: Use short names for color mapping
colors <- c("20240829_1" = "red",
            "20240829_2" = "blue",
            "20240829_3" = "green",
            "P-S4" = "purple",
            "20241127_2" = "orange",
            "P-S3" = "pink",
            "P-S11" = "brown")


if (!dir.exists("2025-2-17-lyq-plot")) {
  dir.create("2025-2-17-lyq-plot")
}

# 循环处理每个文件
for (file in files) {
  data <- read.table(file, header = FALSE, sep = "\t", col.names = c("ID", "Position", "Depth"))
  data$File <- basename(file)
  
  # 1. Normalize
  data <- data %>%
    group_by(ID) %>%
    mutate(Normalized_Position = (Position - min(Position)) / (max(Position) - min(Position))) %>%
    ungroup()
  
  # 获取文件名和缩短的文件名
  file_name <- basename(file)
  short_file_name <- sub("\\.fastqfilt.*", "", file_name)
  
  plot_title <- paste("Depth Plot for All Genes in", short_file_name, "(Normalized)")
  output_file <- paste0("2025-2-17-lyq-plot/", short_file_name, "_normalized_depth.png")
  
  # 2. 图像中处理Y
  p <- ggplot(data, aes(x = Normalized_Position, y = Depth)) +
    geom_line(color = colors[short_file_name]) +  
    facet_wrap(~ ID, scales = "free_y", ncol = 1) +  
    labs(
      title = plot_title,
      x = "Normalized Position (0 to 1)",
      y = "Depth"
    ) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    # background color
    theme(panel.background = element_rect(fill = "white"),  plot.background = element_rect(fill = "white"))
  
  for (gene in unique(data$ID)) {
    gene_range <- data %>% filter(ID == gene) %>% summarize(min=min(Position), max=max(Position))
    subtitle_text <- paste("Original Position Range for", gene, ":", gene_range$min, "-", gene_range$max)
    p <- p + labs(subtitle = subtitle_text)
    break
  }
  
  
  ggsave(output_file, plot = p, width = 10, height = 6 * length(unique(data$ID)) / 3)
}

