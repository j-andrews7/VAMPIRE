import click
#import venusar.thresholds
#import venusar.motifs
#import venusar.activity
#import venusar.tf_expression
#import venusar.gene_expression
import os
import re
import math

#@click.command()
#@click.option('--as-cowboy', '-c', is_flag=True, help='Greet as a cowboy.')
#@click.argument('name', default='world', required=False)

# depends:
#     a configuration file exists; see read_config() for definition or venusar_default.config
# example use (python opened in this path same as config file):
#     import venusar
#    (method_set, var2val, ref2val) = venusar.read_config()
#    process_string = venusar.print_config(method_set, var2val)
#    process_string_t = venusar.print_config(method_set, var2val, translate=True)
#
#

def cli(name, as_cowboy):
    """
    VAMPIRE is a full-fledged bioinformatics suite geared towards easy integration of
    multiple types of big data to make reasonable and interesting biological conclusions.
    """
    greet = 'Howdy' if as_cowboy else 'Hello'
    click.echo('{0}, {1}.'.format(greet, name))

config_file = 'venusar_default.config'


def createConfigLookup():
    """
    return dictionary of lookups to convert config file variable names
    to main variable names for each program

    Caveat: for this to be right any time a new file is added or the
    main of an existing file changed both this lookup and the config
    file must be updated as needed to reflect the addition/changes

    Args:
    Returns:
        dictionary: key = (program_name.py, config_variable_name)
                    value = the abbreviated often obscure main name that matches
                        the verbose config_variable_name

    """
    config_var_lookup = {}    # dictionary

    # note: some of these definitions are not translations but defined for clarity
    #    and to avoid key missing lookup errors in rest of code
    config_var_lookup[('thresholds.py', 'motif_file')] = 'motif_file'
    config_var_lookup[('thresholds.py', 'motif_outfile')] = 'motif_outfile'
    config_var_lookup[('thresholds.py', 'baseline_file')] = 'bp'
    config_var_lookup[('thresholds.py', 'pseudocounts')] = 'pc'
    config_var_lookup[('thresholds.py', 'threshold')] = 'd_th'
    config_var_lookup[('thresholds.py', 'false_pos_rate')] = 'fpr'
    config_var_lookup[('thresholds.py', 'precision_exp')] = 'pe'
    config_var_lookup[('thresholds.py', 'overwrite')] = 'ow'

    config_var_lookup[('tf_expression.py', 'input_file')] = 'inp_file'
    config_var_lookup[('tf_expression.py', 'exp_file')] = 'exp_file'
    config_var_lookup[('tf_expression.py', 'output_file')] = 'out_file'
    config_var_lookup[('tf_expression.py', 'motif_file')] = 'motif_file'
    config_var_lookup[('tf_expression.py', 'motif_out_file')] = 'motifout_file'
    config_var_lookup[('tf_expression.py', 'threshold')] = 'th'

    config_var_lookup[('motifs.py', 'input_file')] = 'file_input'
    config_var_lookup[('motifs.py', 'file_reference')] = 'file_reference_genome'
    config_var_lookup[('motifs.py', 'file_motif')] = 'file_motif'
    config_var_lookup[('motifs.py', 'output_file')] = 'file_output'
    config_var_lookup[('motifs.py', 'baseline_file')] = 'file_baseline_prob'
    config_var_lookup[('motifs.py', 'pseudocounts')] = 'pc'
    config_var_lookup[('motifs.py', 'threshold')] = 'th'
    config_var_lookup[('motifs.py', 'wing_size')] = 'ws'
    config_var_lookup[('motifs.py', 'multi_var')] = 'multivar_distance'
    config_var_lookup[('motifs.py', 'run_homotypic')] = 'run_homotypic'
    config_var_lookup[('motifs.py', 'force_ref_match')] = 'force_ref_match'
    config_var_lookup[('motifs.py', 'chip_file')] = 'file_chip'
    config_var_lookup[('motifs.py', 'chip_out_file')] = 'file_output_chip'
    config_var_lookup[('motifs.py', 'filter_co')] = 'filter_co'
    config_var_lookup[('motifs.py', 'kary_sort')] = 'sorted_lex'
    config_var_lookup[('motifs.py', 'filter_chip')] = 'filter_chip'
    config_var_lookup[('motifs.py', 'filter_motif')] = 'filter_motif'
    config_var_lookup[('motifs.py', 'filter_novel')] = 'filter_novel'

    config_var_lookup[('activity.py', 'input_file')] = 'vcf_file'
    config_var_lookup[('activity.py', 'activity_file')] = 'act_file'
    config_var_lookup[('activity.py', 'output_vcf_file')] = 'out_vcf'
    config_var_lookup[('activity.py', 'output_bed_file')] = 'out_bed'
    config_var_lookup[('activity.py', 'threshold')] = 'thresh'
    config_var_lookup[('activity.py', 'filter_bed_num')] = 'filter_num'
    config_var_lookup[('activity.py', 'include_bed')] = 'include_bed'
    config_var_lookup[('activity.py', 'include_vcf')] = 'include_vcf'
    config_var_lookup[('activity.py', 'drop_activity_')] = 'drop_act_'

    config_var_lookup[('gene_expression.py', 'input_file')] = 'vcf_file'
    config_var_lookup[('gene_expression.py', 'expression_file')] = 'exp_file'
    config_var_lookup[('gene_expression.py', 'output_vcf_file')] = 'out_vcf'
    config_var_lookup[('gene_expression.py', 'wing_size')] = 'size'
    config_var_lookup[('gene_expression.py', 'threshold')] = 'thresh'
    config_var_lookup[('gene_expression.py', 'ethreshold')] = 'ethresh'
    config_var_lookup[('gene_expression.py', 'include_vcf')] = 'include_vcf'

    return(config_var_lookup)


def key2val_or_null(dictionary, get_key):
    """given a dictionary and key return the value or None if key is not defined"""

    try:
        return(dictionary[get_key])
    except:
        return(None)


def by_key(dict_key):
    try:
        return float(dict_key)
    except ValueError:
        return dict_key


def print_config(method_set, var2val, print_actions=True, translate=False):
    """
    Given the output of a read_config variable set,
    print what will be done in order when process_config is called
    By default does NOT translate variable names using createConfigLookup()

    Args:
        method_set: dictionary from read_config
            key = method names. value = list of var2val variables
        var2val: dictionary from read_config
            key = (method_name, variable). value = value to assign
            does not care if values are dereferenced
        print_actions: if True print action string, else just return string
            boolean, default=True
        translate: if True translate the variable names to those used by
            the main function using call to createConfigLookup()

    Returns:
        action_string: the printed content as a string
    """

    var_lookup = createConfigLookup()
    action_string = ""

    # -- determine the run order
    method_order_set = process_config_order(method_set, var2val)

    # -- process printing string by the run order
    process_number = 0
    for method_name in method_order_set:
        process_number = process_number + 1

        action_string = action_string + 'processing ' + str(process_number) + \
            ':' + method_name + '\n'

        # verify method_name is in method_set dictionary (unlikely to not be)
        if key2val_or_null(method_set, method_name) is None:
            action_string = action_string + '\tfailed config check will NOT be processed.\n'
            continue

        # for each variable name in method_set print variable and value
        for var_name in method_set[method_name]:
            if (var_name == 'run' or var_name == 'run_order'):
                continue
            if key2val_or_null(var2val, (method_name, var_name)) is not None:
                if translate:
                    print_var_name = var_lookup[(method_name, var_name)]
                else:
                    print_var_name = var_name
                action_string = action_string + '\t' + print_var_name + '=' + \
                    key2val_or_null(var2val, (method_name, var_name)) + '\n'

        action_string = action_string + '\n'

    if print_actions:
        print(action_string)

    return (action_string)


def process_config(method_set, var2val, debug=False):
    """
    Given the output of a read_config variable set,
    process what will be done in order

    Args:
        method_set: dictionary from read_config
            key = method names. value = list of var2val variables
        var2val: dictionary from read_config
            key = (method_name, variable). value = value to assign
            values MUST be dereferenced, call read_config with replace_ref=True
            because code assumes references in var2val already replaced with ref2val values
        debug: boolean, if True do not execute
            similar to print_config with translate=True but gives actual command string

    Returns:
        boolean for success; True = success; False = fail
    """

    var_lookup = createConfigLookup()
    search_pattern_file = re.compile('file')    # faster to precompile
    search_pattern_str = re.compile('str')      # faster to precompile


    # first process run & run_order variables to determine which methods get called in what order
    # then set variables and call each successive method
    #
    # how notes:
    #    ref: https://stackoverflow.com/questions/8028708/dynamically-set-local-variable
    #    ref: https://stackoverflow.com/questions/2220699/whats-the-difference-between-eval-exec-and-compile-in-python
    #        exec ignores return value
    #        eval returns the value
    # Note: both of the following seem to work if import tf_expression first
    #     exec('tf_expression.main(inp_file="../../data/FLDL_CCCB_RARE_VARIANTS.MERGED.RNA_DP10.RNA_NODUPS.CHIP_MULTIMARK.SORTED.vcf",exp_file="../../data/ALL_ARRAYS_NORMALIZED_MAXPROBE_LOG2_COORDS.sorted.txt",out_file="test_file",motif_file="../../data/HOCOMOCOv10.JASPAR_FORMAT.TF_IDS.fpr_0p001.txt",th=None,motifout_file="test_file_two.vcf")')
    #     exec("tf_expression.main(inp_file='../../data/FLDL_CCCB_RARE_VARIANTS.MERGED.RNA_DP10.RNA_NODUPS.CHIP_MULTIMARK.SORTED.vcf',exp_file='../../data/ALL_ARRAYS_NORMALIZED_MAXPROBE_LOG2_COORDS.sorted.txt',out_file='test_file',motif_file='../../data/HOCOMOCOv10.JASPAR_FORMAT.TF_IDS.fpr_0p001.txt',th=None,motifout_file='test_file_two.vcf')")
    #    also tested; changing order of the arguments; still works

    # -- determine the run order
    method_order_set = process_config_order(method_set, var2val)

    # -- set variables and process
    process_number = 0
    for method_name in method_order_set:
        process_number = process_number + 1

        print(('Start processing setup ' + method_name + ' (' + str(process_number) + ')\n'))

        # verify method_name is in method_set dictionary (unlikely to not be)
        if key2val_or_null(method_set, method_name) is None:
            print('\tfailed config check will NOT be processed.\n')
            continue

        # import the source file to run the main for
        action_string = method_name.split('.py')[0]
        if not debug:
            exec('import ' + action_string)

        # build the action string 1 variable at a time:
        action_string = action_string + '.main('
        first_variable = True
        for var_name in method_set[method_name]:
            if (var_name == 'run' or var_name == 'run_order'):
                continue
            if key2val_or_null(var2val, (method_name, var_name)) is not None:
                if first_variable:
                    first_variable = False
                else:
                    action_string = action_string + ", "

                print_var_name = var_lookup[(method_name, var_name)]

                if (re.search(search_pattern_file, var_name) is not None or
                    re.search(search_pattern_str, var_name) is not None):
                    # only quote if file or str are in variable name
                    action_string = action_string + print_var_name + '=' + \
                    '"' + key2val_or_null(var2val, (method_name, var_name)) + '"'
                else:
                    # do not quote because not file or string
                    action_string = action_string + print_var_name + '=' + \
                        key2val_or_null(var2val, (method_name, var_name))

        action_string = action_string + ')'

        if debug:
            print(action_string)
        else:
            exec(action_string)

    return (True)


def process_config_order(method_set, var2val):
    """
    Given the output of a read_config variable set, determine method call order

    Args:
        method_set: dictionary from read_config
            key = method names. value = list of var2val variables
        var2val: dictionary from read_config
            key = (method_name, variable). value = value to assign
            does not care if values are dereferenced

    Returns:
        call_order: list of method names in order
    """

    run_numbers = {}    # dictionary: keys = call order numbers, value = method name
    for method_name in list(method_set.keys()):
        if key2val_or_null(var2val, (method_name, 'run')) != 'True':
            # if run is False string or not defined (None) then do not run
            continue
        temp_value = key2val_or_null(var2val, (method_name, 'run_order'))
        if temp_value is not None:
            temp_value = int(temp_value)
            if temp_value in run_numbers:
                # run order is already defined: must increment but
                # first determine max prior collision run order value
                max_collision = temp_value
                for run_n_key in list(run_numbers.keys()):
                    if math.floor(run_n_key) == temp_value:
                        if run_n_key > max_collision:
                            max_collision = run_n_key
                # redefine run_order as half distance to next integer from current collision point
                # this allows for multiple collisions defaulting to order encountered in the file
                temp_value = max_collision + (temp_value + 1 - max_collision) / 2
            run_numbers[temp_value] = method_name

    # Could return run_numbers dictionary rather than array of names to allow same
    #    method to be called twice; though actually this is not allowed per read_config
    #    data structure definition using a dictionary for method names
    #    also: additional code would need to be used in processing to verify the correct
    #    entry is called for given method and run order.
    #    if run order collision occurs this could become very difficult.

    # convert to ordered names
    method_order_set = []
    for run_order in sorted(run_numbers, key=by_key):
        method_order_set.append(run_numbers[run_order])

    return(method_order_set)


def read_config(file_name='venusar_default.config', replace_ref=True):
    """
    Read the configuration file and define the data structures used
    to set options for methods called by this wrapper (venusar.py)

    Args:
        file_name: name of the configuration file;
            string, default = venusar_default.config
        replace_ref: replace the reference values in variable=value pairs
            with the reference value.
            QQQ: would one ever not want to do this?
            boolean, default = True

        configuration file format (used to define options passed into called methods):
            # comment line
            <method_name>   # apply successive lines to this method call
            run=True        # run must be specified for each method (MBSFEM), value = True|False
            run_order=1     # run_order MBSFEM, numeric value, smaller values processed first
            var1=value      # the variable names should match items relevant to the method options
            var2=value    REFN    # optional tab separated REFN key
                                  # REFN can be any unique character sequence, prefer: [a-zA-Z0-9]
                                  # REFN character sequence must NOT match any non-REF value
            var3=REFN             # using the previously defined REFN
            <method_name2>
            var1=value

    Returns:

        Three dictionaries:
        method_set; key = method names. value = list of var2val variables
        var2val; key = (method_name, variable). value = value to assign
        ref2val; key = reference name. value = value assigned
    """

    method_set = {}
    var2val = {}
    ref2val = {}

    if not os.path.exists(file_name):
        print(('read_config failed no file for: ' + file_name))
        return (method_set, var2val, ref2val)

    # -- read the file

    try:
        fileHan = open(file_name, "r")
    except:
        print(('read_config failed file open for: ' + file_name))
        return (method_set, var2val, ref2val)

    search_pattern_comment = re.compile('^#')      # faster to precompile
    search_pattern_method_line = re.compile('^<')  # faster to precompile
    current_program_name = ""
    reference_count = 0
    for line in fileHan:
        if len(line) < 2:
            continue
        # process comment lines
        if re.search(search_pattern_comment, line) is not None:
            continue

        # process method lines
        if re.search(search_pattern_method_line, line) is not None:
            line_set = line.strip().split('>')
            temp = line_set[0].strip('<')
            if len(temp) > 0:
                current_program_name = temp
            else:
                print(('Failed reset program name at line: ', line))
            continue

        # process variable lines
        # var = value
        # var = value    refN
        # var = refN
        # testline: line = 'help me = two five\tREFN'
        line = line.strip().split('#')[0]    # drop post comment information
        line_set = line.split('\t')          # split on tab for references
        var_set = line_set[0].split('=')     # split variable = value

        if len(var_set) > 1:
            # add variable=value pair; but first check keys for duplicate variable!
            var_name = var_set[0].rstrip().lstrip()     # variable w/o leading/trailing spaces
            var_value = var_set[1].rstrip().lstrip()    # value w/o leading/trailing spaces
            if (current_program_name, var_name) in var2val:
                # ZZZ: how to handle multiple var keys if for different methods?
                #     double up the key as a tuple; lists do NOT work
                print(('WARNING: duplicate value for ' + current_program_name + "::"
                    + var_name + " using current value."))
            var2val[(current_program_name, var_name)] = var_value
            if current_program_name in method_set:
                method_set[current_program_name].append(var_name)    # reading subsequent var
            else:
                method_set[current_program_name] = [var_name]    # reading first var for program
        else:
            continue

        if len(line_set) > 1:
            # a reference exists
            line_set[1] = line_set[1].rstrip().lstrip()
            if line_set[1] in ref2val:
                print(('WARNING: the reference is already defined! Replacing ' + ref2val))
            else:
                reference_count = reference_count + 1
            ref2val[line_set[1]] = var_set[1]

    # -- reduce cross-references
    cref = reference_count
    count_change = False
    while cref > 0:
        count_change = False    # reset
        for ref_key in list(ref2val.keys()):
            # check each key's value being a key itself
            if ref2val[ref_key] in ref2val:
                # referenced a reference! replace if not self-referential
                if not(ref2val[ref2val[ref_key]] == ref2val[ref_key]):
                    ref2val[ref_key] = ref2val[ref2val[ref_key]]
                    cref = cref - 1
                    count_change = True

        if not count_change:
            break

    # -- convert ref values in var2val to values
    if replace_ref:
        cref = reference_count
        while cref > 0:
            count_change = False    # reset
            for key_pair in list(var2val.keys()):
                # check the value of var2val being a reference key
                # then if a reference key, replace with the value
                if var2val[key_pair] in ref2val:
                    if ref2val[var2val[key_pair]] not in ref2val:
                        var2val[key_pair] = ref2val[var2val[key_pair]]
                        cref = cref - 1
                        count_change = True
            if not count_change:
                break

    return (method_set, var2val, ref2val)