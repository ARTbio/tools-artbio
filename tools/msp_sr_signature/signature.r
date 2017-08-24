      graph_type = "${graph_type}"

      globalgraph = function () {
        ## Setup R error handling to go to stderr
        options( show.error.messages=F,
                 error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
        signature = read.delim("${output}", header=TRUE)
        signaturez=data.frame(signature[,1], (signature[,2] -mean(signature[,2]))/sd(signature[,2]))
        overlap_prob_z=data.frame(signature[,1], (signature[,3] -mean(signature[,3]))/sd(signature[,3]))
        YLIM=max(signature[,2])
        

        ## Open output2 PDF file
        pdf( "${output2}" )
        if (YLIM!=0) {
          par(mfrow=c(2,2),oma = c(0, 0, 3, 0))

          plot(signature[,1:2], type = "h", main="Numbers of pairs", cex.main=1, xlab="overlap (nt)", ylim=c(0,YLIM), ylab="Numbers of pairs", col="darkslateblue", lwd=4)

          plot(signaturez, type = "l", main="Number of pairs Z-scores", cex.main=1, xlab="overlap (nt)", ylab="z-score", pch=19, cex=0.2, col="darkslateblue", lwd=2)

          plot(signature[,1], signature[,3]*100, type = "l", main="Overlap probabilities",
             cex.main=1, xlab="overlap (nt)", ylab="Probability [%]", ylim=c(0,50),
             pch=19, col="darkslateblue", lwd=2)

          plot(overlap_prob_z, type = "l", main="Overlap Probability Z-scores", cex.main=1, xlab="overlap (nt)", ylab="z-score", pch=19, cex=0.2, col="darkslateblue", lwd=2)

          mtext("Overlap Signatures of ${minquery}-${maxquery} against ${mintarget}-${maxtarget}nt small RNAs", outer = TRUE, cex=1)
        }
        devname = dev.off()
        ## Close the PDF file
      }

      treillisgraph = function () {
        ## Open output2 PDF file
        pdf( "${output2}", paper="special", height=11.69, width=8.2677 )
        signature = read.delim("${output}", header=TRUE)
        options( show.error.messages=F,
               error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
        library(lattice)
        print (xyplot(signature[,3]*100~signature[,1]|signature[,4], type = "l", xlim=c(${minscope},${maxscope}), main="ping-pong Signature of ${minquery}-${maxquery} against ${mintarget}-${maxtarget}nt small RNAs",
             par.strip.text=list(cex=.5), strip=strip.custom(which.given=1, bg="lightblue"), scales=list(cex=0.5),
             cex.main=1, cex=.5, xlab="overlap (nt)", ylab="ping-pong signal [%]",
             pch=19, col="darkslateblue", lwd =1.5, cex.lab=1.2, cex.axis=1.2,
             layout=c(4,12), as.table=TRUE, newpage = T) )
        devnname = dev.off()
      }

      if (graph_type=="global") {
        globalgraph()

      }
      if(graph_type=="lattice") {
        treillisgraph()
      }

