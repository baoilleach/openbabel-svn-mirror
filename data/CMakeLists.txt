# 2008-03-15

set(to_install
    aromatic.txt
    atomtyp.txt
    bondtyp.txt
    element.txt
    fragments.txt
    gaff.dat
    gaff.prm
    ghemical.prm
    isotope-small.txt
    isotope.txt
    logp.txt
    MACCS.txt
    mmff94.ff
    mmff94s.ff
    mmffang.par
    mmffbndk.par
    mmffbond.par
    mmffchg.par
    mmffdef.par
    mmffdfsb.par
    mmffoop.par
    mmffpbci.par
    mmffprop.par
    mmffstbn.par
    mmfftor.par
    mmffvdw.par
    mmffs_oop.par
    mmffs_tor.par
    mr.txt
    patterns.txt
    phmodel.txt
    plugindefines.txt
    psa.txt
    qeq.txt
    resdata.txt
    ringtyp.txt
    SMARTS_InteLigand.txt
    space-groups.txt
    superatom.txt
    svgformat.script
    torlib.txt
    types.txt
    UFF.prm
)

if(NOT MSVC)
  install(FILES ${to_install} DESTINATION share/openbabel/${BABEL_VERSION})
else(NOT MSVC)
  install(FILES ${to_install} DESTINATION bin/data)
endif(NOT MSVC)
