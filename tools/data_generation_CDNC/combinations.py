from itertools import product
"""
by   Atte Laakso / UEF

Forms different combinations of given parameters
"""

# for running all combinations
def make_combos(rhsel, mmrsel, nmrsel, refsel):
    
    # define options for combination template
    options = ['yes', 'no']
    
    # Determine which options are possible to be set to specific positions
    enabled_options = [rhsel, mmrsel, nmrsel, refsel]
    
    # Calculate the number of enabled options
    rep = sum(1 for option in enabled_options if option in ['rh', 'mmr', 'nmr', 'refrac'])
    
    # Generate combinations for the enabled options
    combos = list(product(options, repeat=rep))
    updated_combos = []

    # replace 'yes'-selections with corresponding characters
    for comb in combos:
        comb = list(comb)
        full_comb = ['no', 'no', 'no', 'no']
        
        enabled_index = 0
        # check all the members in a single parameter-combination
        for i, option in enumerate(enabled_options):
            # check if the option is in the right form
            if option in ['rh', 'mmr', 'nmr', 'refrac']:
                # check if the option is possible to be inserted
                if comb[enabled_index] == 'yes':
                    full_comb[i] = option
                # move on to the next member of this parameter combination
                enabled_index += 1

        # add modified combo to the list of combinations in their final forms
        updated_combos.append(full_comb)
    
    return updated_combos

# for running all combinations
def make_CDNC_combos(mmrsel, nmrsel, updraftsel, cltsel, clwsel):
    
    # define options for combination template
    options = ['yes', 'no']
    
    # Determine which options are possible to be set to specific positions
    enabled_options = [mmrsel, nmrsel, updraftsel, cltsel, clwsel]
    
    # Calculate the number of enabled options
    rep = sum(1 for option in enabled_options if option in ['mmr', 'nmr', 'updraft', 'clt', 'clw'])
    
    # Generate combinations for the enabled options
    combos = list(product(options, repeat=rep))
    updated_combos = []

    # replace 'yes'-selections with corresponding characters
    for comb in combos:
        comb = list(comb)
        full_comb = ['no', 'no', 'no', 'no', 'no']
        
        enabled_index = 0
        # check all the members in a single parameter-combination
        for i, option in enumerate(enabled_options):
            # check if the option is in the right form
            if option in ['mmr', 'nmr', 'updraft', 'clt', 'clw']:
                # check if the option is possible to be inserted
                if comb[enabled_index] == 'yes':
                    full_comb[i] = option
                # move on to the next member of this parameter combination
                enabled_index += 1

        # add modified combo to the list of combinations in their final forms
        updated_combos.append(full_comb)
    
    return updated_combos