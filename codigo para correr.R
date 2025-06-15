library(seqinr)
getwd()

alinea = function(gen,gen2){
  m = matrix(data=0, nrow=length(gen2)+1, ncol=length(gen)+1)
  m[1, ] = seq(0,-(length(gen))*2,-2)
  m[ ,1] = seq(0,-(length(gen2))*2,-2)
  
  for (fila in seq(2,nrow(m))){
    for (col in seq(2,ncol(m))){
      if (gen[col-1]==gen2[fila-1])
        diag = m[fila-1,col-1] + 1
      else
        diag = m[fila-1,col-1] - 1
      up = m[fila-1,col] - 2
      left = m[fila,col-1] - 2
      m[fila,col] = max(diag,up,left)
    }  
  }
  
  fila = length(gen2)+1
  col = length(gen)+1
  solA = solB = c()
  while(fila>1 || col>1){
    if (col>1 && fila>1 && gen[col-1]==gen2[fila-1]){
      solA = c(gen[col-1], solA)
      solB = c(gen2[fila-1], solB)
      col = col - 1
      fila = fila -1
    }else{
      if (col>1 && fila>1 && m[fila-1, col-1]>m[fila, col-1] && m[fila-1, col-1]>m[fila-1, col]){
        solA = c(gen[col-1], solA)
        solB = c(gen2[fila-1], solB)
        col = col - 1
        fila = fila -1
      }else if (fila==1 || (col>1 && m[fila, col-1] > m[fila-1, col])){
        solA = c(gen[col-1], solA)
        solB = c("_", solB)
        col = col - 1
      }else if (col==1 || (fila>1 && m[fila, col-1] <= m[fila-1, col])){
        solA = c("_", solA)
        solB = c(gen2[fila-1], solB)
        fila = fila - 1 
      }
    }
  }
  sol = matrix(data=c(solA,solB),nrow = 2, ncol = length(solA), byrow = TRUE)
  return (sol)
}

trad = c(UUU="F", UUC="F", UUA="L", UUG="L",
         UCU="S", UCC="S", UCA="S", UCG="S",
         UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
         UGU="C", UGC="C", UGA="STOP", UGG="W",
         CUU="L", CUC="L", CUA="L", CUG="L",
         CCU="P", CCC="P", CCA="P", CCG="P",
         CAU="H", CAC="H", CAA="Q", CAG="Q",
         CGU="R", CGC="R", CGA="R", CGG="R",
         AUU="I", AUC="I", AUA="I", AUG="M",
         ACU="T", ACC="T", ACA="T", ACG="T",
         AAU="N", AAC="N", AAA="K", AAG="K",
         AGU="S", AGC="S", AGA="R", AGG="R",
         GUU="V", GUC="V", GUA="V", GUG="V",
         GCU="A", GCC="A", GCA="A", GCG="A",
         GAU="D", GAC="D", GAA="E", GAG="E",
         GGU="G", GGC="G", GGA="G", GGG="G")

datos = data.frame(
  mutacion = character(),
  cambioCodon = character(),
  cambioAmino = character(),
  pos = integer(),
  gen = character()
)

file = read.fasta("sequence.txt", forceDNAtolower = FALSE)
file2 = read.fasta("ba52 chn.fasta", forceDNAtolower = FALSE)
cat(length(file)%/%12, "vs", length(file2)/12, "secuencias \n")
vs = length(file2)/12

# Procesar solo el gen S (número 3)
i = 3
gen = file[[i]]
info = attr(gen,"Annot")
info = unlist(strsplit(info,"\\[|\\]|:|=|\\."))
gene = info[which(info=="gene")+1]
cat("Procesando Gen", i, gene, "\n")
gen[which(gen=="T")] = "U"
cat("Total de nucleótidos (Wuhan):", length(gen), "\n")

nMut = 1
secuencias_procesadas = 0  # Cambio: usar contador más preciso
secuencias_con_indels = 0
indels_mostrados = 0  # Para limitar mensajes de diagnóstico

# Diagnóstico inicial
cat("Total de secuencias en file2:", length(file2), "\n")
cat("Número estimado de muestras:", length(file2) / 12, "\n")

# Procesar todas las secuencias del gen S
# Cada muestra tiene 12 genes, el gen S está en la posición 3 de cada conjunto
indices_gen_s = seq(3, length(file2), by = 12)
cat("Índices del gen S a procesar:", head(indices_gen_s, 10), "...\n")
cat("Total de secuencias del gen S disponibles:", length(indices_gen_s), "\n")

# Limitar a 100 secuencias
max_secuencias = min(100, length(indices_gen_s))
cat("Procesando", max_secuencias, "secuencias del gen S\n")

for (idx in 1:max_secuencias) {
  j = indices_gen_s[idx]
  
  gen2 = file2[[j]]
  gen2[which(gen2=="T")] = "U"
  secuencias_procesadas = secuencias_procesadas + 1  # Cambio: incrementar siempre
  
  if (length(gen) == length(gen2)){
    # Secuencias de igual longitud - mutaciones puntuales
    diff = which(gen!=gen2)
    if (length(diff)>0){ 
      prevMut = ""
      for (pos in diff){
        ini = pos - (pos-1)%%3
        mutacion = paste(gen[pos], "to", gen2[pos], sep="")
        codOri = paste(gen[ini],gen[ini+1],gen[ini+2],sep="")
        codMut = paste(gen2[ini],gen2[ini+1],gen2[ini+2],sep="")
        codonChange = paste(codOri,"to",codMut,sep="")
        nCod = ((pos-1)%/%3) + 1
        # Verificar si los codones son válidos antes de procesar
        if (!is.na(trad[codOri]) && !is.na(trad[codMut])) {
          aminoChange = paste(trad[codOri],nCod,trad[codMut],sep="")
          
          # Filtrar mutaciones a codones STOP (opcional - comentar si quieres incluirlas)
          # if (trad[codMut] == "STOP") {
          #   cat("ADVERTENCIA: Mutación a STOP en posición", pos, ":", codOri, "->", codMut, "\n")
          # }
          
          if (trad[codOri] != trad[codMut] && prevMut != aminoChange) {
            datos[nMut, ] = list(mutacion, codonChange, aminoChange, nCod, gene)
            nMut = nMut + 1
          }
          prevMut = aminoChange
        } else {
          cat("ADVERTENCIA: Codón inválido encontrado:", codOri, "o", codMut, "en posición", pos, "\n")
        }
      }
    }
  }else{
    # Secuencias de diferente longitud - INDELs
    # CAMBIO PRINCIPAL: Remover el límite de 10 secuencias con INDELs
    secuencias_con_indels = secuencias_con_indels + 1
    
    # Solo mostrar mensaje de diagnóstico para los primeros 5 casos
    if (indels_mostrados < 5) {
      cat("Mutación de inserción o deleción \n")
      cat("Wuhan:", length(gen), "vs:", length(gen2),"\n")
      indels_mostrados = indels_mostrados + 1
    }
    
    res = alinea(gen, gen2)
    diff = which(res[1,]!=res[2,])
    if (length(diff)>0 && indels_mostrados <= 5) {
      cat("Mutaciones en las posiciones:", head(diff, 10), "...\n")  # Limitar salida
    }
    
    prevCod = ""
    align_gen = res[1, ]
    align_gen2 = res[2, ]
    for (k in seq(1, length(align_gen) - 2, by = 3)) {
      cod1 = paste(align_gen[k:(k + 2)], collapse = "")
      cod2 = paste(align_gen2[k:(k + 2)], collapse = "")
      
      if (grepl("_", cod1) || grepl("_", cod2)) next
      
      if (!is.na(trad[cod1]) && !is.na(trad[cod2]) && trad[cod1] != trad[cod2]) {
        cambio = paste(trad[cod1], k %/% 3 + 1, trad[cod2])
        if (cambio != prevCod) {
          if (indels_mostrados <= 5) {
            cat(cod1, "to", cod2, cambio, k %/% 3 + 1, gene, "\n")
          }
          datos[nMut, ] = list(paste(cod1, "to", cod2),
                               paste(cod1, "to", cod2),
                               cambio, k %/% 3 + 1, gene)
          nMut = nMut + 1
          prevCod = cambio
        }
      }
    }
  }
}

cat("Resumen para gen S\n")
cat("Total de secuencias procesadas:", secuencias_procesadas, "\n")  # Cambio: usar contador correcto
cat("Secuencias con inserciones/deleciones:", secuencias_con_indels, "\n")
cat("Porcentaje con INDELs:", round(secuencias_con_indels / secuencias_procesadas * 100, 2), "%\n")  # Cambio: usar contador correcto

str(datos)

library(dplyr)

dfgraph = filter(
  summarise(
    select(
      group_by(datos, cambioAmino),
      mutacion:gen
    ),
    mutacion = first(mutacion),
    cambioCodon = first(cambioCodon),
    pos = first(pos),
    gen = first(gen),
    cuenta = n()
  ),
  cuenta > as.integer(secuencias_procesadas*0.5)  # Cambio: usar secuencias_procesadas
)
cat("Solo frecuencias superiores al 10%:",as.integer(secuencias_procesadas*0.1),"\n")  # Cambio: usar secuencias_procesadas
dfgraph

library(ggplot2)

# Gráfico de mutaciones nucleotídicas
p = ggplot(datos)
p = p + aes(x=mutacion, fill=mutacion, label=after_stat(count))
p = p + ggtitle(paste("Cambios de nucleótidos en gen S - Wuhan vs",secuencias_procesadas,"secuencias"))  # Cambio: usar secuencias_procesadas
p = p + labs(x="Mutación", y="Frecuencia", fill="Mutación")
p = p + geom_bar(stat="count")
p = p + geom_text(stat="count", vjust=0)
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

# Gráfico de cambios de aminoácidos
p3 = ggplot(dfgraph)
p3 = p3 + aes(x=cambioAmino, y=cuenta, fill=cambioAmino, label=cuenta)
p3 = p3 + ggtitle(paste("Cambio de aminoácidos en gen S - Wuhan vs",secuencias_procesadas,"secuencias"))  # Cambio: usar secuencias_procesadas
p3 = p3 + labs(x="Cambio de aminoácido", y="Frecuencia", fill="Frecuencia")
p3 = p3 + geom_bar(stat = "identity")
p3 = p3 + geom_text(stat = "identity", vjust=1.5)
p3 = p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3