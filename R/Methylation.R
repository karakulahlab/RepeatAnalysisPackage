
toMakeMethObj<-function(path, id,assmbly,treat,context){
  myobj=methRead(path,
                 sample.id=id,
                 assembly=assmbly,
                 treatment=treat,
                 context=context
  )
  return(myobj)
}



toCalcRepeats<-function(obj_meth,gr_repeats){
  repeats=regionCounts(myobj,gr_repeats)
  return(repeats)
}


