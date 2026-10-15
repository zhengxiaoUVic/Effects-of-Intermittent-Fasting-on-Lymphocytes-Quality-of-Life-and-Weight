library(lattice)

df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
df$Phase <- factor(df$Phase)
png(filename="Fig.1A.png", width=750*3, height=700*4, res=300)
xyplot(LYMP ~ Time | ID, data=df, group=Phase, type=c("p","r"), layout=c(3,5),
       index.cond=list(c(11,15,5,
                         9,14,4,
                         8,13,3,
                         7,12,2,
                         6,10,1)), #change panel order
       xlab=list(label="Month", cex=1.2), 
       ylab=list(label=expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")), cex=1.2),
       sub=list(label="   Left: \"3+3\"                  Center: \"6\"                  Right: \"3+1\"", cex=1.2))
dev.off()


library(magick) #macOS needs to install XQuartz, https://www.xquartz.org/
fig1a <- image_read("~/Documents/IF/Fig.1A.png")
fig1b <- image_read("~/Documents/IF/Fig.1B.png")
fig1 <- image_append(c(fig1a, fig1b))
fig1 <- image_annotate(fig1, "A", size=100, weight=700)
fig1 <- image_annotate(fig1, "B", size=100, weight=700, location="+2300")
image_write(fig1, path="~/Documents/IF/Fig.1.png", format="png")

