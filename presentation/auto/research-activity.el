(TeX-add-style-hook "research-activity"
 (lambda ()
    (TeX-add-symbols
     '("footlineextra" 1)
     "frame")
    (TeX-run-style-hooks
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

