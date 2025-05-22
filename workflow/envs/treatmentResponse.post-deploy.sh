#!/usr/bin/env bash
set -o pipefail

R -e "pak::pkg_install(c('bhklab/CoreGx', 'bhklab/PharmacoGx', 'cancerrxgene/gdscIC50', 'bhklab/AnnotationGx'))"
