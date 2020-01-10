from pyrosetta import *
from pyrosetta.rosetta import *

def MakePyRosettaMoverFromXML(xmlfile,flags_fname,pdbname):
    '''
    Main XML-parsing function. Returns mover.
    '''

    def ProcessFlags(flags_fname):
        '''
        Helper function for properly detecting options in flag file.
        '''
        handle = open(flags_fname,'r')
        flags = [ i[:-1] for i in handle.readlines() ]
        handle.close()
        replaces = []
        # Variable replacment in the XML is handled separately:
        for line in flags:
            if "parser:script_vars" in line:
                replaces.append(line.split()[1])
        # Return all flags and variable replacements separately
        return (flags,replaces)

    # Prepare flags for input
    parsed_flags = ProcessFlags(flags_fname)
    # Initialize Rosetta
    init(extra_options=" ".join(parsed_flags[0]))
    pose = pose_from_file(filename=pdbname)
    # Now actually create mover
    parser = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()

    # Necessary magic step:
    options = basic.options.process()
    #####################

    # Manage variable replacements:
    replaces = pyrosetta.rosetta.utility.vector1_string(0)

    for i in parsed_flags[1]:
        replaces.append(i)
    ####################

    modified_pose = False
    tag = parser.create_tag_from_xml( xmlfile, replaces)
    in_mover = parser.generate_mover_for_protocol( pose, modified_pose, tag, options )

    return (in_mover,pose)
