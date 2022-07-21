#!/bin/bash

set -ex

if test -f "./h5m_files/ARE.h5m"; then
  python ./scripts/ARE_geometry_plot.py
else
  python ./scripts/step_to_h5m.py
  python ./scripts/ARE_geometry_plot.py
fi
