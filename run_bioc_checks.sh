#!/usr/bin/env bash
set -euo pipefail

rm -f NAMESPACE

echo "== ExpoRiskR: document =="
Rscript -e 'devtools::document()'

echo "== ExpoRiskR: run tests/check =="
Rscript -e 'devtools::check()'

echo "== ExpoRiskR: BiocCheck =="
Rscript -e 'BiocCheck::BiocCheck()'

echo "== ExpoRiskR: build tar.gz =="
R CMD build .

echo "Done."
