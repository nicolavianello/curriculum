(TeX-add-style-hook "research-activity"
 (lambda ()
    (LaTeX-add-bibliographies
     "biblio")
    (TeX-add-symbols
     '("footlineextra" 1)
     "frame")
    (TeX-run-style-hooks
     "biblatex"
     "tikz"
     "pgf"
     "rfxcolor"
     "tangocolors"
     "multimedia"
     "amsmath"
     "listings"
     "babel"
     "english"
     "latex2e"
     "beamer10"
     "beamer"
     "t"
     "10pt")))

