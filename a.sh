cd 1d
pandoc -s --filter pandoc-crossref -M "crossrefYaml=../crossref_config.yaml" --template=../github.html --mathjax --highlight-style=tango --metadata pagetitle="1d FEM" README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html

cd ../GLQ
pandoc -s --filter pandoc-crossref -M "crossrefYaml=../crossref_config.yaml" --template=../github.html --mathjax --highlight-style=tango --metadata pagetitle="Gauss-Legendre quadrature" README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html

cd ../Einstein
pandoc -s --filter pandoc-crossref -M "crossrefYaml=../crossref_config.yaml" --template=../github.html --mathjax --highlight-style=tango --metadata pagetitle="Summation convention" README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html

cd ../Deform
pandoc -s --filter pandoc-crossref -M "crossrefYaml=../crossref_config.yaml" --template=../github.html --mathjax --highlight-style=tango --metadata pagetitle="Deformation" README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html

cd ../
pandoc -s --filter pandoc-crossref -M "crossrefYaml=crossref_config.yaml" --template=github.html --mathjax --highlight-style=tango --metadata pagetitle="Computational mechanics" README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html



