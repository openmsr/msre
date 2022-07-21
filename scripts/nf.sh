#!/bin/bash

set -ex

if test -f "./h5m_files/ARE.h5m"; then
  python ./scripts/ARE_neutron_flux.py
else
  python ./scripts/step_to_h5m.py
  python ./scripts/ARE_neutron_flux.py
fi
