#!/bin/bash

set -ex

if test -f "./h5m_files/msre.h5m"
then
  python ./scripts/msre.py
else
  python ./scripts/step_to_h5m.py &&
  python ./scripts/msre.py
fi
