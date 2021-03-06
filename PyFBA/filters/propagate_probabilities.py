import sys
import PyFBA

def role_probs_to_reaction_probs(role_probs, verbose=False):
    """
    Convert functional role probabilities to reaction probabilities
    using the model seed mappings.

    For a hash of roles and their probabilities, return a hash of reaction ids
    and thier associated probabilities.

    :param role_probs: a hash of roles and their associated probabilities
    :type role_probs: dict of float
    :param verbose: print error reporting
    :type verbose: bool
    :return: a hash of reaction ids and their associated probabilities
    :rtype: dict of float
    """

    # Model seed mappings:
    # key is role and value is all complexes for the role
    role_2_cmplx = PyFBA.parse.model_seed.roles()
    # key is complex and value is all reactions for the complex
    cmplx_2_rxn = PyFBA.parse.model_seed.complexes()

    # Map role probabilities to enzyme complex probabilities
    cmplx_probs = {}
    for role in role_probs:
        for cmplx in role_2_cmplx[role]:
            if cmplx not in cmplx_probs:
                cmplx_probs[cmplx] = [role_probs[role]]
            else:
                cmplx_probs[cmplx].append(role_probs[role])

    # Write out enzyme complexes and probabilites to file
    fout = open('genome_enzyme_complex_probabilities.txt','w')
    fout.write('ENZYME COMPLEX\tPROBABILITY\n')
    for cmplx in cmplx_probs:
        cmplx_probs[cmplx] = min(cmplx_probs[cmplx])
        fout.write(cmplx + '\t' + str(cmplx_probs[cmplx]) + '\n')
    fout.close()
        

    # Map enzyme complex probabilities to reaction probabilities
    rxn_probs = {}
    for cmplx in cmplx_probs:
        if cmplx not in cmplx_2_rxn:
            if verbose:
                # This can occur because there are some complexes that don't
                # have an associated reaction
                sys.stderr.write("ERROR: " + cmplx + " was not found in the complexes file.")
            continue
        for rxn in cmplx_2_rxn[cmplx]:
            if rxn not in rxn_probs:
                rxn_probs[rxn] = [cmplx_probs[cmplx]]
            else:
                rxn_probs[rxn].append(cmplx_probs[cmplx])

    # Write out
    fout = open('genome_reaction_probabilities.txt','w')
    fout.write('REACTION\tPROBABILITY\n')
    for rxn in rxn_probs:
        rxn_probs[rxn] = max(rxn_probs[rxn])
        fout.write(rxn + '\t' + str(rxn_probs[rxn]) + '\n')
    fout.close()

    return rxn_probs
        

    
                

    
