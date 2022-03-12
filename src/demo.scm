;; demo.scm

(import (pdb parser))

(define pdb
  (parse-pdb "../data/3l0u.pdb"))

(define mol
  (pdb->molecule pdb))

(define atom-0
  (vector-ref (molecule-atoms mol) 0))
