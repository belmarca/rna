# RNA annotation in Scheme

This library requires a recently compiled version of Gambit Scheme.

Tell Gambit to search for libraries in the current directory:

``` sh
cd rna/src
gsi lib/
```

which allows you to import and use the library:

``` scheme
(import (pdb parser))
(define pdb (parse-pdb "../data/3l0u.pdb"))
(define mol (pdb->molecule pdb))
(vector-length (molecule-atoms mol))
1573
```
