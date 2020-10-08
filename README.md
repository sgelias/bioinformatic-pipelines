# Bioinformatic Fundamentals

Note: The original version is provided by Prof. Waldeyr Mendes Cordeiro da Silva, see on [Github](https://github.com/waldeyr).

This tutorial contain four main parts:

1. **01_environment_preparation**: contains all scripts to prepare a local or docker machine for bioinformatic analysis;

2. **02_de_novo_assembly**: genome assembly without reference sequence;

3. **03_with_reference_assembly**: genome assembly with reference sequence;

4. **04_rna_seq**: explore RNAseq results.

If you're interested in trying this out, just clone the repository and execute the command `pandoc *.md > markdown_book.html` in the root directory of the repo ([details for multi.md files repository preparation and compilation](https://github.com/akmassey/markdown-multiple-files-example)). Additional parameter would be supplied to produce a user friendly document:

```bash
pandoc -V geometry:margin=0.85in --highlight-style=breezedark *.md > markdown_book.html
```
