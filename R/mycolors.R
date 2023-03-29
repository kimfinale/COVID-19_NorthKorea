# these colors were derived from the stackoverflow (https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes)
mycolors25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

mycolors16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")

SB_alpha <- scales::alpha(c("steelblue"), alpha = c(0.4, 0.7, 0.9)) # symptomatic
Br_alpha <- scales::alpha(c("brown"), alpha = c(0.5)) # data
Gr_alpha <- scales::alpha(c("darkgreen"), alpha = c(0.4, 0.7, 0.9)) # infection
