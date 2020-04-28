MOGAMUN TODO
============

## Github

- write a proper README.md

## Code

- GeneExpression file and Layers directory: add checks in load.data to nicely abort with informative error message if no such file or directory,
- NodesScores file: improve the documentation to explain that if no such file exist, the file will be created. If the file exists it will be used instead of being recomputed,
-  GeneExpression, NodesScores, Layers files: add checks and informative error message if files exist but their format is wrong,
- Layers = "123", etc. : improve documentation,
- ResultsDir: add check if directory already exist, to remove spurious Warning
- Re-add time to the results directory name (beware that ':' is forbidden in Windows file paths)

## RStudio

## Vignette

- add a proper vignette. Have to read the documentation

## Test

- add proper test. Have to read the documentation
