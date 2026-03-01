#!/bin/bash

function clean {
  file=$1
  out=$2

  sed 's/lstlisting/algorithm/g' $file > tmp_a
  sed 's/\[language=C\]//g' tmp_a > tmp_b
  sed 's/\[language=\[90\]FORTRAN}\]//g' tmp_b > tmp_a
  sed 's/\[language=Python\]//g' tmp_a > tmp_b
  sed 's/|//g' tmp_b > tmp_a
  sed 's/\\label{/\\label{algo-/g' tmp_a > tmp_b
  #sed 's/|\\zero|/0/g' tmp_b > tmp_a
  #sed 's/|\\un|/1/g' tmp_a > tmp_b
  #sed 's/|\\deux|/2/g' tmp_b > tmp_a
  #sed 's/|\\trois|/3/g' tmp_a > tmp_b
  #sed 's/|\\mbtt{27}|/27/g' tmp_b > tmp_a
  #sed 's/|\\vbtt{x\\_pos}|/x_pos/g' tmp_a > tmp_b
  #sed 's/|\\rbtt{y\\_pos}|/y_pos/g' tmp_b > tmp_a
  #sed 's/|\\obtt{z\\_pos}|/z_pos/g' tmp_a > tmp_b
  #sed 's/|\\cdtt{aid}|/aid/g' tmp_b > tmp_a
  #sed 's/|\\bftt{nnp}|/nnp/g' tmp_a > tmp_b
  #sed 's/\\textbar\\textbar//g' tmp_b > tmp_a
  #sed 's/|\\vbtt{pix}|/pix/g' tmp_a > tmp_b
  #sed 's/|\\rbtt{pjx}|/pjx/g' tmp_b > tmp_a
  #sed 's/|\\dbtt{aid}|/aid/g' tmp_a > tmp_b
  #sed 's/|\\cbtt{bid}|/bid/g' tmp_b > tmp_a
  #sed 's/|\\cbtt{aid}|/aid/g' tmp_a > tmp_b
  #sed 's/|\\vbtt{xpos}|/xpos/g' tmp_b > tmp_a
  #sed 's/|\\rbtt{ypos}|/ypos/g' tmp_a > tmp_b
  #sed 's/|\\obtt{zpos}|/zpos/g' tmp_b > tmp_a
  #sed 's/%\/\//\/\//g' tmp_a > tmp_b
  #sed 's/|\\textquotesingle|/'"'"'/g' tmp_b > tmp_a
  #sed 's/|\\dbtt{axis}|/axis/g' tmp_a > tmp_b
  # sed 's/\\AA\\/Angstrom/g' tmp_a > tmp_b

  rm -f tmp_a
  mv tmp_b $out
}


clean 'c-code.tex' 'c-aglo.tex'
clean 'f-code.tex' 'f-algo.tex'
clean 'p-code.tex' 'p-algo.tex'

