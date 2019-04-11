set(PROGS 
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
endforeach(name)
