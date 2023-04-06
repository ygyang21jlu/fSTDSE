#!/bin/bash
file="$1"
output="$2"

awk -v output_file=${output} \
'BEGIN { infty = 1e20 ; go = infty ; }
/Large amplitudes of individual field-free states/ { go = NR + 4 }
(NR>=go)&&(NF<6){go=infty;next;}
(NR>=go){
    print   $1,$2,$3,$4,$5,$6>> output_file;
  }
END {
  }
' ${file}
