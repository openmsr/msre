###############################################################################
# Converting step files to h5m file to be read by openmc
###############################################################################
import numpy as np
import os
import sys
sys.path.append('../CAD_to_openMC/src/')
import CAD_to_OpenMC.assembly as ab
###############################################################################

# inputs
step_filepath = "./step_files/msre_parts.step"
h5m_out_filepath = os.getcwd() + '/h5m_files/msre.h5m'

# mesher config
ab.mesher_config['min_mesh_size'] = 0.2
ab.mesher_config['mesh_algorithm'] = 4
ab.mesher_config['threads'] = 32
#ab.mesher_config['max_mesh_size'] = 1000
#ab.mesher_config['curve_samples'] = 100
#ab.mesher_config['vetoed'] = [474,1157,1243,1341,1537]
#ab.mesher_config['angular_tolerance'] = 0.09
#ab.mesher_config['tolerance'] = 0.001
ab.mesher_config['refine']=False
# output
a=ab.Assembly()
a.remove_intermediate_files=True
a.stp_files=[step_filepath]
tags={'hold-down.*':'graphite',\
  'thimble.*':'inor', 'core\ support':'graphite',\
  'graphite':'graphite',\
  'core.*':'inor',\
  'vessel':'inconel',\
  '.*stringer.*':'graphite',\
  'inor.*':'inor',\
  '(bottom|top|basket)':'inconel',\
  '(lower|upper) lattice':'graphite',\
  'lifting stud':'inor',\
  'distributor':'inor',\
  'restrictor':'inor',\
  '(support|lower) ring':'inor',\
  '^(left|right|back|front)$':'graphite',\
  '.*grid support':'inor'}
a.import_stp_files(tags=tags, match_anywhere=True)
a.solids_to_h5m(backend='gmsh',h5m_filename=h5m_out_filepath)
