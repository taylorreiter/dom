library(scales)
library(RColorBrewer)

S <- "#F8E629"
Q <- "#C0D2DD"
P <- "#8CC4E1"
I <- "#28A4D7"
H <- "#088AC9"
F <- "#0055A3" #
E <- "#2D3392" #0254A2
G <- "#5A5398"
C <- "#363189"
Y <- "#CEBB30"
O <- "#F6B344"
U <- "#F2D2C6"
W <- "#EDA17B"
N <- "#E68D63"
M <- "#E6482C"
T <- "#CE3735"
V <- "#B8525B"
D <- "#8E2E38"
B <- "#B9C771"
L <- "#98AB43"
K <- "#859065"
A <- "#6E7559"
J <- "#414C30"

unk <- c(S)
#metabolism <- c(Q, P, I , H, F, E, G, C)
cell <- c(Y, O, U, W, N, M, T, V, D)
info <- c(B, L, K, A, J)

metabolism <- brewer.pal(9, "Blues")[2:9]
# cell <- brewer.pal(9, "Oranges")[2:9]
# info <- brewer.pal(9, "YlGn")[c(2, 3, 4, 6, 8)]

cog_pal  <- c(unk, metabolism, cell, info)
show_col(cog_pal)
