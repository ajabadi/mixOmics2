exs <- list.files('examples', 'example-*.R$', full.names = TRUE)
names(exs) <- gsub('-example.R','',basename(unlist(exs)))

runtimes <- rep(-1,length(exs))
names(runtimes) <- names(exs)

for (ex in names(exs)[4:6]){
  runtimes[ex] <-
    tryCatch({
      round(system.time({
        source(exs[[ex]], ech=FALSE )
      })['elapsed'])},
      error=function(e){e}
    )

}
if(is(runtimes, 'list')) runtimes <- unlist(runtimes)

write.table(runtimes, '__exclude/examples-runetime.txt',
            sep = "\t", row.names = TRUE, col.names = 'runtime')
