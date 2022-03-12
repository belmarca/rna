;; pdb parser library
;; (C) 2022 Marc-André Bélanger

(define max-cols 80)
(define record-cols 6)

(define pdb-records
  '("HEADER"
    "TITLE "
    "COMPND"
    "AUTHOR"
    "JRNL  "
    "DBREF "
    "DBREF1"
    "DBREF1"
    "SEQADV"
    "SEQRES"
    "CRYST1"
    "ORIGX1"
    "ORIGX2"
    "ORIGX3"
    "SCALE1"
    "SCALE2"
    "SCALE3"
    "ATOM  "
    "HETATM"
    "CONECT"))

(define-structure pdb
  ;; PDB textual data
  header
  title
  compnd
  author
  jrnl
  dbref
  dbref1
  dbref2
  seqadv
  seqres
  cryst1
  origx1
  origx2
  origx3
  scale1
  scale2
  scale3
  atom
  hetatm
  conect
  ;; Computed
  natom
  nhetatm
  )

(define (make-pdb*)
  (make-pdb "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" 0 0))

(define (get-record line)
  (substring line 0 record-cols))

(define (valid-record record)
  (member record pdb-records))

(define (end? record)
  (string=? record "END   "))

(define (parse-pdb path)
  (call-with-input-file path parse-pdb*))

(define (parse-pdb* port)
  (define (read-line*)
    (read-line port #\newline #f (+ max-cols 1)))
  (let ((pdb (make-pdb*)))
    (let loop ((line (read-line*)))
      (if (not (eq? line #!eof))
          (let ((record (get-record line)))
            (if (not (end? record))
                (if (valid-record record)
                    (begin
                      (parse-record pdb port record line) ;; Could run sub-parser
                      (loop (read-line*)))
                    (loop (read-line*)))))))
    pdb))

(define (parse-record pdb port record line)
  (define-macro (gen-parser)
    (define records
      '(("HEADER" . s)
        ("TITLE " . m)
        ("COMPND" . m)
        ("AUTHOR" . m)
        ("JRNL  " . m)
        ("DBREF " . m)
        ("DBREF1" . m)
        ("DBREF1" . m)
        ("SEQADV" . m)
        ("SEQRES" . m)
        ("CRYST1" . s)
        ("ORIGX1" . s)
        ("ORIGX2" . s)
        ("ORIGX3" . s)
        ("SCALE1" . s)
        ("SCALE2" . s)
        ("SCALE3" . s)
        ("ATOM  " . m)
        ("HETATM" . m)
        ("CONECT" . m)))

    (define (remove-spaces r)
      (let loop ((i 0))
        (if (< i 6) ;; record-cols
            (if (eq? (string-ref r i) #\space)
                (substring r 0 i)
                (loop (+ i 1)))
            (substring r 0 i))))

    (define (gen-set! r)
      (string->symbol
       (string-downcase
        (string-concatenate
         (list "pdb-" (remove-spaces r) "-set!")))))

    (define (gen-get r)
      (string->symbol
       (string-downcase
        (string-concatenate
         (list "pdb-" (remove-spaces r))))))

    (define (gen-clauses)
      (map
       (lambda (f)
         (cond
          ((string=? (car f) "ATOM  ")
           `((string=? record "ATOM  ")
             (begin
               (pdb-atom-set! pdb (string-concatenate
                                   (list (pdb-atom pdb) line "\n")))
               (pdb-natom-set! pdb (+ (pdb-natom pdb) 1)))))
          ;; NOTE: Don't process HETATM for now.
          ;; ((string=? (car f) "HETATM")
          ;;  `((string=? record "HETATM")
          ;;    (begin
          ;;      (pdb-hetatm-set! pdb (string-concatenate
          ;;                            (list (pdb-hetatm pdb) line "\n")))
          ;;      (pdb-nhetatm-set! pdb (+ (pdb-nhetatm pdb) 1)))))
          (else
           (if (eq? (cdr f) 's)
               `((string=? record ,(car f))
                 (,(gen-set! (car f)) pdb line))
               `((string=? record ,(car f))
                 (,(gen-set! (car f)) pdb (string-concatenate
                                           (list (,(gen-get (car f)) pdb) line "\n"))))))))
       records))
    `(cond ,@(gen-clauses)))
  (gen-parser))

(define-structure atom
  type
  serial-number
  name
  alternate-location-indicator
  residue-name
  chain-identifier
  residue-sequence-number
  code-for-insertion-of-residues
  x-coordinate
  y-coordinate
  z-coordinate
  occupancy
  temperature-factor
  segment-identifier
  element-symbol
  charge)

;; PDB ATOM and HETATM record parser
(define (parse-atom line)
  (let ((type (substring line 0 6))
        (serial-number (substring line 6 11))
        (name (substring line 12 16))
        (alternate-location-indicator (substring line 16 17))
        (residue-name (substring line 17 20))
        (chain-identifier (substring line 21 22))
        (residue-sequence-number (substring line 22 26))
        (code-for-insertion-of-residues (substring line 26 27))
        (x-coordinate (substring line 30 38))
        (y-coordinate (substring line 38 46))
        (z-coordinate (substring line 46 54))
        (occupancy (substring line 54 60))
        (temperature-factor (substring line 60 66))
        (segment-identifier (substring line 72 76))
        (element-symbol (substring line 76 78))
        (charge (substring line 78 80)))
    (make-atom
     type
     serial-number
     name
     alternate-location-indicator
     residue-name
     chain-identifier
     residue-sequence-number
     code-for-insertion-of-residues
     x-coordinate
     y-coordinate
     z-coordinate
     occupancy
     temperature-factor
     segment-identifier
     element-symbol
     charge)))

(define-structure molecule
  ;; We group ATOM with HETATM records
  pdb
  atoms)

(define (pdb->molecule pdb)
  (define (read-line* port)
    (read-line port #\newline #f (+ max-cols 1)))

  (let* ((natom (pdb-natom pdb))
         (nhetatm (pdb-nhetatm pdb))
         (natoms (+ natom nhetatm))
         (mol (make-molecule pdb (make-vector natoms))))

    ;; Parse ATOM records
    (call-with-input-string
     (pdb-atom pdb)
     (lambda (port)
       (let loop ((line (read-line* port)) (i 0))
         (if (< i natom)
             (begin
               (vector-set! (molecule-atoms mol) i (parse-atom line))
               (loop (read-line* port) (+ i 1)))))))

    ;; Parse HETATM records
    (call-with-input-string
     (pdb-hetatm pdb)
     (lambda (port)
       (let loop ((line (read-line* port)) (i natom))
         (if (< i natoms)
             (begin
               (vector-set! (molecule-atoms mol) i (parse-atom line))
               (loop (read-line* port) (+ i 1)))))))
    mol))

;; Either remove leading blanks, or coalesce an all #\space string to
;; the empty string. Find the first occurring non-blank character
;; and return the suffix substring.
(define (remove-leading-blanks s)
  (if (string=? s "")
      ""
      (let ((len (string-length s)))
        (let loop ((i 0))
          (if (< i len)
              (if (eq? (string-ref s i) #\space)
                  (loop (+ i 1))
                  (substring s i len))
              ;; We received only blanks
              "")))))

;; Export molecule atomic information to a JavaScript array of objects
(define (atom->js-object atom)
  (let ((data (map remove-leading-blanks
                   (list (atom-x-coordinate atom)
                         (atom-y-coordinate atom)
                         (atom-z-coordinate atom)
                         (atom-serial-number atom)
                         (atom-residue-sequence-number atom)
                         (atom-element-symbol atom)
                         (atom-residue-name atom)
                         (atom-chain-identifier atom)
                         (atom-name atom)))))
    (string-concatenate
     (list "{"
           "x: " (list-ref data 0) ","
           "y: " (list-ref data 1) ","
           "z: " (list-ref data 2) ","
           "serialNumber: " (list-ref data 3) ","
           "residueSequenceNumber: " (list-ref data 4) ","
           "element: '" (list-ref data 5) "',"
           "residueName: '" (list-ref data 6) "',"
           "chain: '" (list-ref data 7) "',"
           "name: \"" (list-ref data 8) "\","
           "}"))))

(define (molecule->js-object mol)
  (define (atoms->js-list atoms)
    (string-concatenate
     (fold (lambda (a b)
             (cons (atom->js-object a) (cons "," b)))
           '()
           (vector->list atoms))))
  ;; Currently exports minimal information, maybe add SEQRES and others
  (let ((pdb (molecule-pdb mol)))
    (string-concatenate
     (list "{"
           "cryst1: '" (pdb-cryst1 pdb) "',"
           "origx1: '" (pdb-origx1 pdb) "',"
           "origx2: '" (pdb-origx2 pdb) "',"
           "origx3: '" (pdb-origx3 pdb) "',"
           "scale1: '" (pdb-scale1 pdb) "',"
           "scale2: '" (pdb-scale2 pdb) "',"
           "scale3: '" (pdb-scale3 pdb) "',"
           "atoms: [" (atoms->js-list (molecule-atoms mol)) "]"
           "}"
           ))))

(define (export-molecule-as-js mol dest)
  (with-output-to-file dest
    (lambda ()
      (display
       (string-concatenate
        (list "const molecule = " (molecule->js-object mol) ";"))))))
