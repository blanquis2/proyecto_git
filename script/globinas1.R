
                          #Alineamiento de globinas#


#cargar los paquetes
library(BSgenome)
library(BiocIO)
library(GenomicRanges)
library(rtracklayer)
library(msa)

#agregar el archivo al doc

globinas<-readDNAStringSet("datos/DivergentGlobins.fasta")
globinas

#alineamiento

clustal_globinas<-msa(globinas,method = "ClustalW")
clustal_globinas

muscle_globinas<-msa(globinas, method = "Muscle")
muscle_globinas

#árboles

#cargar el paquete

library(ape)

install.packages("seqinr")
library(seqinr)

#Cambiar el tipo de objeto#

globinas1<-msaConvert(clustal_globinas,type="seqinr::alignment")
globinas1


globinas2<-msaConvert(muscle_globinas,type="seqinr::alignment")
globinas2

#hacer las matrices de distancias
mat1<-dist.alignment(globinas1)
mat2<-dist.alignment(globinas2)
  
#árbol#
pdf("resultado/arboles.pdf")

arbol1<-nj(mat1)#con clustalw
arbol1
plot(arbol1)
arbol2<-nj(mat2)#con muscle
arbol2
plot(arbol2)
dev.off()

#Ejercicios con ggtree#
install.packages("ggtree")
library(ggtree)

class(arbol1)#verificar que la clase sea phylo
class(arbol2)

###ggtree

#Árbol random

plot(rtree(10))->arbolito

#Se puede personalizar agregando capas

ggtree(arbolito)+
  geom_tiplab(size=5, colour="bisque2")

#resaltar clados
geom_cladelab()

#señalar con color ciertos clados
geom_hilight()

#guadar el gráfico
ggsave("resultado/arbolito.png")#elegir el formato favorito

#asignar nombres
class(arbolito)
str(arbol1)

#arbolito1 
ggsave("resultado/arbolclustalw.png")
arbolito1<- ggtree(arbol1, layout = 'circular', branch.length='none')+
  geom_tiplab(colour="bisque3")+
  geom_point(aes(shape=isTip, color=isTip), size=3)
  
arbolito1
dev.off()

##arbolito2
ggsave("resultado/arbolmuscle.png")
arbolito2<- ggtree(arbol2, layout = 'circular', branch.length='none')+
  geom_tiplab(colour="bisque2")+
  geom_point(aes(shape=isTip, color=isTip), size=3)
arbolito2
dev.off()

#Prueba
prueba1<- ggtree(arbol1) + geom_tippoint(aes(size = 3), x = arbol$Nnode, shape = 1) + 
  geom_tippoint(aes(size = tarsusL), x=arbol$tip.label[1], shape = 1) + 
  geom_tippoint(aes(size = culmenL), x = arbol$tip.label[2], shape = 1) + 
  geom_tippoint(aes(size = beakD),   x = arbol$tip.label[3], shape = 1) + 
  geom_tippoint(aes(size = gonysW),  x = arbol$tip.label[4], shape = 1) +
  geom_tippoint(aes(size = gonysW),  x = arbol$tip.label[5], shape = 1)+ 
  scale_size_continuous(range = c(3,12), name="globinas") + 
  geom_text(aes(x = x, y = 0, label = lab), angle = 45) +
  geom_tiplab(offset = 1.3) + xlim(0, 3) +
  theme(legend.position = c(.1, .75)) + vexpand(.05, -1) 

prueba2<-ggtree(arbol1,layout="fan", open.angle=120)
prueba2+
  geom_rootpoint(colour="red3",size=3, alpha=.8)+
  geom_tiplab(colour="blue4")+
  geom_treescale(x=0, y=40, width=1, color='purple3')

prueba3<-ggtree(arbol1,layout_inward_circular(xlim=5))
prueba3#no sale

arbolito<-rtree(20)

ggtree(arbolito)+
  geom_tiplab(colour="bisque3")+
  geom_point(aes(shape=isTip, color=isTip), size=3)
