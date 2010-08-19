import os
import pybel
import time

from optparse import OptionParser

def main(input_ext, inputfile, output_ext, outputfilename,
         nconfs, rmsd_cutoff, energy_cutoff):
    
    ff = pybel._forcefields['mmff94']
    outputfile = pybel.Outputfile(output_ext, outputfilename, overwrite=True)
    for i, mol in enumerate(pybel.readfile(input_ext, inputfile)):
        t = time.time()
        
        print "**Molecule %d\n..title = %s" % (i, mol.title)
        print "..number of rotatable bonds = %d" % mol.OBMol.NumRotors()
        mol.addh()
        ff.Setup(mol.OBMol)
        ff.DiverseConfGen(rmsd_cutoff, nconfs, energy_cutoff)

        ff.GetConformers(mol.OBMol)
        confdata = pybel.ob.toConformerData(mol.OBMol.GetData(pybel.ob.ConformerData))
        energies = confdata.GetEnergies()

        N = mol.OBMol.NumConformers()
        assert N == len(energies)
        print "..generated %d conformers"
        
        u = time.time()
        data = []
        for i in range(N):
            mol.OBMol.SetConformer(i)
            outputfile.write(mol)

        print "..(overall time = %.1fs  writing results = %.1fs)" % (time.time() - t,
                                                                  time.time() -u)
        print "\n"
        
    outputfile.close()

def display_options(myd):
    for x, y in myd.iteritems():
        print "..%s = %s" % (x, y)

if __name__ == "__main__":
    # Set up the options parser
    usage = "usage: %prog [options] inputfile outputfile"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="inputformat",
                      help="input file format", metavar="FILE")
    parser.add_option("-o", dest="outputformat",
                      help="output file format", metavar="FILE")
    parser.add_option("-r", "--rmsd", dest="rmsd_cutoff", default=0.5,
                      type="float",
                      help="RMSD cutoff (default 0.5)", metavar="RMSD")
    parser.add_option("-e", "--energy", dest="energy_cutoff", default=50.0,
                      type="float", metavar="Energy",
                      help="Energy cutoff (default 50.0 kcal/mol)")
    parser.add_option("-c", "--confs", dest="conf_cutoff", default=int(1E6),
                      type="int",
                      help="max number of conformers to test (default is 1 million)", metavar="RMSD")
##    parser.add_option("--seed", dest="random_seed", default=1,
##                      help="random seed", metavar="RMSD")
    (options, args) = parser.parse_args()

    # Error handling
    if len(args) != 2:
        parser.error("You need to specify two arguments, the input file and the output file")
    inputfile, outputfile = args[0], args[1]
    if not os.path.isfile(inputfile):
        parser.error("Cannot find the input file, '%s'" % inputfile)
    
    input_ext = options.inputformat
    if not input_ext:
        input_ext = inputfile.split(".")[-1]

    if input_ext not in pybel.informats.keys():
        parser.error("The input file format is unknown.\n\n  Please specify a valid format with '-i'.")

    output_ext = options.outputformat
    if not output_ext:
        output_ext = outputfile.split(".")[-1]
    if output_ext not in pybel.outformats.keys():
        parser.error("The output file format is unknown.\n\n  Please specify a valid format with '-o'.")

    # Call main
    print "\n**Starting confab 0.1"
    display_options({'inputfile':inputfile, 'outputfile':outputfile})
    display_options(eval(str(options)))
    print "\n"
    main(input_ext, inputfile, output_ext, outputfile, options.conf_cutoff,
         options.rmsd_cutoff, options.energy_cutoff)
    print "**Finished\n"
    
