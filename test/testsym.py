"""Test OpenBabel executables from Python

Note: Python bindings not used

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pytest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../../test/testsym.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

import os
import unittest

from testbabel import run_exec, executable, log, BaseTest

class TestSym(BaseTest):
    """Base class for a series of tests relating to symmetry"""

    def testInChItoSMI(self):
        """Verify that the InChI is read correctly"""
        output, error = run_exec(self.inchi, "babel -iinchi -ocan")
        self.assertEqual(output.rstrip(), self.cansmi)

    def testSMItoInChI(self):
        """Verify that all molecules give the same InChI"""
        output, error = run_exec("\n".join(self.smiles), "babel -ismi -oinchi")
        output = "\n".join([x.rstrip() for x in output.split("\n")])
        self.assertEqual(output.rstrip(), "\n".join([self.inchi] * len(self.smiles)))

    def testSMItoCAN(self):
        """Verify that all molecules give the same cansmi"""
        output, error = run_exec("\n".join(self.smiles), "babel -ismi -ocan")
        output = "\n".join([x.rstrip() for x in output.split("\n")])
        self.assertEqual(output.rstrip(), "\n".join([self.cansmi] * len(self.smiles)))

    def testSMIthruXML(self):
        """Verify that roundtripping through CML preserves stereo"""
        output, error = run_exec("\n".join(self.smiles), "babel -ismi -ocml tmp.cml")
        output, error = run_exec(output.rstrip(), "babel -icml tmp.cml -ocan")
        output = "\n".join([x.rstrip() for x in output.split("\n")])
        self.assertEqual(output.rstrip(), "\n".join([self.cansmi] * len(self.smiles)))
	os.remove("tmp.cml")

class TestTetSym(TestSym):
    """A series of tests relating to tetrahedral symmetry"""

    def setUp(self):
        self.canFindExecutable("babel")

        # The following all represent the same molecule
        self.cansmi = "C[C@@](F)(Cl)Br"
        self.inchi = "InChI=1S/C2H3BrClF/c1-2(3,4)5/h1H3/t2-/m0/s1"
        self.smiles = [
             'C[C@@](Cl)(Br)F',
             'C[C@](Cl)(F)Br',
             'C[C@](Br)(Cl)F',
             'C[C@@](Br)(F)Cl',
             'C[C@@](F)(Cl)Br',
             'C[C@](F)(Br)Cl',
             'Cl[C@](C)(Br)F',
             'Cl[C@@](C)(F)Br',
             'Cl[C@@](Br)(C)F',
             'Cl[C@](Br)(F)C',
             'Cl[C@](F)(C)Br',
             'Cl[C@@](F)(Br)C',
             'Br[C@@](C)(Cl)F',
             'Br[C@](C)(F)Cl',
             'Br[C@](Cl)(C)F',
             'Br[C@@](Cl)(F)C',
             'Br[C@@](F)(C)Cl',
             'Br[C@](F)(Cl)C',
             'F[C@](C)(Cl)Br',
             'F[C@@](C)(Br)Cl',
             'F[C@@](Cl)(C)Br',
             'F[C@](Cl)(Br)C',
             'F[C@](Br)(C)Cl',
             'F[C@@](Br)(Cl)C'
             ]


class TestCisTransSym(TestSym):
    """A series of tests relating to cistrans symmetry"""

    def setUp(self):
        self.canFindExecutable("babel")

        # The following all represent the same molecule
        self.cansmi = "Cl/C=C/C=C\Br"
        self.inchi = "InChI=1S/C4H4BrCl/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+"
        self.smiles = [
                "C(=C\C=C/Br)/Cl",
                "Cl/C=C/C=C\Br", 
                "Br/C=C\C=C\Cl",
                "C(=C\Cl)/C=C\Br",
                "C(=C\C=C\Cl)\Br",
                "C(=C\Br)\C=C\Cl"
                ]

class TestConversions(BaseTest):
    """A series of tests relating to file format conversions and symmetry"""
    
    def setUp(self):
        self.canFindExecutable("babel")
        
    def testSMILEStoInChI(self):
        # Tests interconversions between the SMILES on the left versus
        # the InChI on the right.
        # The canonical smiles (in the middle) were derived from the SMILES.
        data = [
('ClC=CF', 'FC=CCl',       'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H'),
('Cl/C=C/F', 'F/C=C/Cl',   'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1+'),
('Cl/C=C\\F', 'F/C=C\\Cl', 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1-'),
('Cl[C@@](Br)(F)I', 'F[C@@](Cl)(Br)I', 'InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m0/s1'),
('Cl[C@](Br)(F)I', 'F[C@](Cl)(Br)I',   'InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m1/s1'),
('ClC(Br)(F)I', 'FC(Cl)(Br)I',         'InChI=1S/CBrClFI/c2-1(3,4)5')]
        for smiles, can, inchi in data:
            output, error = run_exec(smiles, "babel -ismi -oinchi")
            self.assertEqual(output.rstrip(), inchi)
            output, error = run_exec(inchi, "babel -iinchi -ocan")
            self.assertEqual(output.rstrip(), can)
            
    def parseMDL(self, text):
        lines = text.split("\n")
        broken = lines[3].split()
        Natoms = int(broken[0])
        Nbonds = int(broken[1])
        atoms = []
        for i in range(Natoms):
            broken = lines[i+4].split()
            atoms.append({'parity':int(broken[6])})
        bonds = []
        for i in range(Nbonds):
            broken = lines[i+4+Natoms].split()
            bonds.append({'stereo':int(broken[3])})
        return atoms, bonds
        
    def testSMILESto3DMDL(self):
        """Test interconversion between SMILES and 3D MDL"""
        data = [
('ClC=CF', 'FC=CCl', [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3]),
('Cl/C=C/F', 'F/C=C/Cl', [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
('Cl/C=C\\F', 'F/C=C\\Cl', [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0]), 
('Cl[C@@](Br)(F)I', 'F[C@@](Cl)(Br)I', [0, 0, 0, 0, 1], [0, 0, 0, 1]),
('Cl[C@](Br)(F)I', 'F[C@](Cl)(Br)I', [0, 0, 0, 0, 2], [0, 0, 0, 6]),
('ClC(Br)(F)I', 'FC(Cl)(Br)I', [0, 0, 0, 0, 3], [0, 0, 0, 4])
]
        for smiles, can, atompar, bondstereo in data:
            output, error = run_exec(smiles, "babel -ismi -osdf --gen3d")
            atoms, bonds = self.parseMDL(output)
            parities = [atom['parity'] for atom in atoms]
            parities.sort()
            stereos = [bond['stereo'] for bond in bonds]
            stereos.sort()
            self.assertEqual(atompar, parities)
            self.assertEqual(bondstereo, stereos)
            output, error = run_exec(output, "babel -isdf -ocan")
            self.assertEqual(output.rstrip(), can)

    def test2DMDLto0D(self):
        """Test conversion for 2D MDL to CAN and InChI"""
        # The following file was created using RDKit starting from
        # the SMILES strings in data[x][0] below.
        filename = self.getTestFile("testsym_2Dtests.sdf")
        data = [
('ClC=CF', 'FC=CCl',       'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H'),
('Cl/C=C/F', 'F/C=C/Cl',   'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1+'),
('Cl/C=C\\F', 'F/C=C\\Cl', 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1-'),
('Cl[C@@](Br)(F)I', 'F[C@@](Cl)(Br)I', 'InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m0/s1'),
('Cl[C@](Br)(F)I', 'F[C@](Cl)(Br)I',   'InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m1/s1'),
('ClC(Br)(F)I', 'FC(Cl)(Br)I',         'InChI=1S/CBrClFI/c2-1(3,4)5')]
        
        output, error = run_exec("babel -isdf %s -ocan" % filename)
        for i, smiles in enumerate(output.rstrip().split("\n")):
            self.assertEqual(smiles.rstrip(), data[i][1])

        output, error = run_exec("babel -isdf %s -oinchi" % filename)
        for i, inchi in enumerate(output.rstrip().split("\n")):
            self.assertEqual(inchi.rstrip(), data[i][2])

class TestStereoConversion(BaseTest):
    """Random tests relating to roundtripping stereochemistry"""
    def setUp(self):
        self.canFindExecutable("babel")
    def testInChIToSMILES_Bug(self):
        """PR#2101034- InChI <-> SMILES conv misrepresents stereo"""
        test_inchi = 'InChI=1S/C10H10/c1-2-3-7-10-8-5-4-6-9-10/h2-9H,1H2/b7-3+'
        output, error = run_exec(test_inchi, "babel -iinchi -ocan")
        self.assertEqual(output.rstrip(), "C=C/C=C/c1ccccc1")
        
        test_smiles = "C=C\C=C/c1ccccc1"
        output, error = run_exec(test_smiles, "babel -ismi -oinchi")
        self.assertEqual(output.rstrip(), "InChI=1S/C10H10/c1-2-3-7-10-8-5-4-6-9-10/h2-9H,1H2/b7-3-")
        
if __name__ == "__main__":
    testsuite = []
    for myclass in [TestConversions, TestCisTransSym, TestTetSym, TestStereoConversion]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    unittest.TextTestRunner().run(unittest.TestSuite(testsuite))
