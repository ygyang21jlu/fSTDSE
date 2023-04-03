#!/bin/bash
for inp in hydrogen_1S_2P0_uniform.inp hydrogen_1S_2P0.inp hydrogen_2P0_ion.inp \
           argon_3P1_cooper.inp argon_3P1_offcooper.inp argon_3P1m_circ_l.inp \
           argon_3P1m_circ_r.inp hydrogen_1S_hhg_linear.inp hydrogen_2P0_sfi.inp \
           hydrogen_1S_hhg_elliptical.inp ; do
  out="$(echo "${inp}" | sed -e 's/\.inp/.out/')"
  if [ -r "${out}" ] ; then
    echo "${out} already exists; can't run ${inp}"
  else
    echo "Executing ${inp}"
    export HUGETLB_MORECORE=thp
    export OMP_STACKSIZE=500M
    ulimit -s 1024000
    ../spherical_tdse.x < "${inp}" > "${out}" 2>&1 
  fi
done
