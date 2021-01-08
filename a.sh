cd 1d
pandoc -s --filter pandoc-crossref -M "crossrefYaml=../crossref_config.yaml" --template=../github.html --mathjax --highlight-style=tango README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html

cd ../Einstein
pandoc -s --filter pandoc-crossref -M "crossrefYaml=../crossref_config.yaml" --template=../github.html --mathjax --highlight-style=tango README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html

cd ../
pandoc -s --filter pandoc-crossref -M "crossrefYaml=crossref_config.yaml" --template=github.html --mathjax --highlight-style=tango README.md -o index.html && sed 's/README.md/index.html/' index.html > tmp && mv tmp index.html



