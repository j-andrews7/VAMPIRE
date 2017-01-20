import click
import venusar.thresholds
import venusar.motifs
import venusar.activity
import venusar.tf_expression
import venusar.gene_expression
import os
import re


@click.command()
@click.option('--as-cowboy', '-c', is_flag=True, help='Greet as a cowboy.')
@click.argument('name', default='world', required=False)
def cli(name, as_cowboy):
    """
    Driven is a full-fledged bioinformatics suite geared towards easy integration of multiple types of big data to make
    reasonable and interesting biological conclusions.
    """
    greet = 'Howdy' if as_cowboy else 'Hello'
    click.echo('{0}, {1}.'.format(greet, name))

config_file = 'venusar_default.config'


def read_config(file_name='venusar_default.config', replace_ref=True):
    """
    Read the configuration file and define the data structures used
    to set options for methods called by this wrapper (venusar.py)

    Args:
        file_name: name of the configuration file;
            default = venusar_default.config
        replace_ref: r
            default = True

        configuration file format (used to define options passed into called methods):
            # comment line
            <method_name>    # apply successive lines to this method call
            var1=value       # the variable names should match items relevant to the method options
            var2=value    REFN    # optional tab separated REFN key
                                  # REFN can be any unique character sequence, prefer: [a-zA-Z0-9]
                                  # REFN character sequence must NOT match any non-REF value
            var3=REFN             # using the previously defined REFN
            <method_name2>
            var1=value

    Returns:

        Three dictionaries
        method_set; key = method names. value = list of var2val variables
        var2val; key = (method_name, variable). value = value to assign
        ref2val; key = reference name. value = value assigned
    """

    print('XXX: incomplete method linked to auto config using setattr to assign variables for each method; needs review and testing')

    method_set = {}
    var2val = {}
    ref2val = {}

    if not os.path.exists():
        print(('read_config failed no file for: ' + file_name))
        return (var2val, ref2val)

    # -- read the file

    try:
        fileHan = open(file_name, "r")
    except:
        return (var2val, ref2val)

    search_pattern_comment = re.compile('^#')      # faster to precompile
    search_pattern_method_line = re.compile('^<')  # faster to precompile
    current_program_name = ""
    reference_count = 0
    for line in fileHan:
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
        line_set = line.split('\t')              # split on tab for references
        var_set = line_set[0].split('=')     # split variable = value
        # variable w/o leading/trailing spaces: var_set[0].rstrip().lstrip()
        # value w/o leading/trailing spaces: var_set[1].rstrip().lstrip()
        if len(var2val) > 0:
            # variable=value pair; but first check keys for duplicate variable!
            if (current_program_name, var_set[0]) in var2val:
                # QQQ: how to handle multiple var keys if for different methods!
                #     double up the key as a tuple; lists do NOT work
                print(('WARNING: duplicate value for ' + current_program_name + "::"
                    + var_set[0] + " using current value."))
                var2val[(current_program_name, var_set[0])] = var_set[1]
            if current_program_name in method_set:
                method_set[current_program_name].append(var_set[0])
            else:
                method_set[current_program_name] = [var_set[0]]
        else:
            continue

        if len(line_set) > 1:
            # a reference exists
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
            if ref2val[ref2val[ref_key]] in ref2val:
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
                if var2val[key_pair] in ref2val:
                    if ref2val[var2val[key_pair]] not in ref2val:
                        var2val[key_pair] = ref2val[var2val[key_pair]]
                        cref = cref - 1
                        count_change = True
            if not count_change:
                break

    return (method_set, var2val, ref2val)