(define-library (annotate)
  (import (gambit))
  (export pdb-parse
          pdb->molecule

          ;; Data structures
          pdb-header
          pdb-title
          pdb-compnd
          pdb-author
          pdb-jrnl
          pdb-dbref
          pdb-dbref1
          pdb-dbref2
          pdb-seqadv
          pdb-seqres
          pdb-cryst1
          pdb-origx1
          pdb-origx2
          pdb-origx3
          pdb-scale1
          pdb-scale2
          pdb-scale3
          pdb-atom
          pdb-hetatm
          pdb-conect
          pdb-natom
          pdb-nhetatm

          atom-type
          atom-serial-number
          atom-name
          atom-alternate-location-indicator
          atom-residue-name
          atom-chain-identifier
          atom-residue-sequence-number
          atom-code-for-insertion-of-residues
          atom-x-coordinate
          atom-y-coordinate
          atom-z-coordinate
          atom-occupancy
          atom-temperature-factor
          atom-segment-identifier
          atom-element-symbol
          atom-charge

          molecule-atoms
          )

  (include "annotate.scm"))
