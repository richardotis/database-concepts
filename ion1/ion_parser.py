from pycalphad import Database, variables as v
from pycalphad.io.tdb import _sympify_string
from collections import defaultdict
import amazon.ion.simpleion as ion
import sympy

def parse_element(name: str):
    return name

def create_cond(cond: ion.IonPyList):
    op_map = {
        '<=': sympy.LessThan,
        '<': sympy.StrictLessThan,
        '>=': sympy.GreaterThan,
        '>': sympy.StrictGreaterThan
    }
    if len(cond) == 5:
        first, first_op, middle, second_op, last = cond
        first = float(first)
        first_op = op_map[first_op.text]
        middle = v.__dict__[middle.text]
        second_op = op_map[second_op.text]
        last = float(last)
        return first_op(first, middle) & second_op(middle, last)
    elif len(cond) == 3:
        first, op, last = cond
        first = v.__dict__.get(middle.text, None) or float(first)
        op = op_map[op.text]
        last = v.__dict__.get(middle.text, None) or float(last)
        return op(first, last)
    else:
        raise ValueError('Unrecognized condition: ', str(cond))

def create_expr(expr: ion.IonPySymbol):
    return _sympify_string(expr.text)

def create_piecewise(conds, exprs):
    conds = [create_cond(c) for c in conds] + [True]
    exprs = [create_expr(e) for e in exprs] + [0]
    return sympy.Piecewise(*zip(exprs, conds))

def parse_function(function_def: ion.IonPyList):
    conds = []
    exprs = []
    for cond, expr in (function_def[i:i+2] for i in range(0, len(function_def), 2)):
        conds.append(cond)
        exprs.append(expr)
    piecewise = create_piecewise(conds, exprs)
    return piecewise

def parse_parameter_gibbs(dbf, param_def):
    param_type = param_def.ion_annotations[0].text.upper()
    phase_name = param_def.ion_annotations[1].text.upper()
    interaction_def, function_def = param_def
    constituent_def, interaction_order = interaction_def
    constituent_array = []
    for subl in constituent_def:
        s = []
        for x in subl:
            s.append(v.Species(x.text))
        constituent_array.append(tuple(s))
    constituent_array = tuple(constituent_array)
    function_obj = parse_function(function_def)

    ref = None # TODO
    diffusing_species = None # TODO
    dbf.add_parameter(param_type, phase_name,
                      [[str(c) for c in sorted(lx)] for lx in constituent_array],
                      interaction_order, function_obj, ref, diffusing_species, force_insert=False)

def parse_parameter_quadruplet(dbf, param_def):
    param_type = param_def.ion_annotations[0].text.upper()
    phase_name = param_def.ion_annotations[1].text.upper()
    constituent_def, quadruplet_def = param_def
    constituent_array = []
    for subl in constituent_def:
        s = []
        for x in subl:
            s.append(x.text)
        constituent_array.append(tuple(s))
    constituent_array = tuple(constituent_array)
    quadruplet_dict = {}
    for key, value in quadruplet_def.items():
        try:
            value = value.text
        except AttributeError:
            pass
        quadruplet_dict[v.Species(key)] = sympy.sympify(value)

    ref = None # TODO
    diffusing_species = None # TODO
    # TODO
    #dbf.add_parameter(param_type, phase_name,
    #                  [[c.upper() for c in sorted(lx)] for lx in constituent_array],
    #                  interaction_order, function_obj, ref, diffusing_species, force_insert=False)

def parse_phase_cef(phase_def: ion.IonPyDict):
    constituent_array = []
    for subl in phase_def['constituents']:
        s = []
        for x in subl:
            s.append(x.text)
        constituent_array.append(tuple(s))
    constituent_array = tuple(constituent_array)
    return {'sites': list(phase_def['sites']),
            'constituent_array': constituent_array}

phase_context_map = defaultdict(lambda: parse_phase_cef)

param_context_map = {
    'G': parse_parameter_gibbs,
    'L': parse_parameter_gibbs,
    'Z': parse_parameter_quadruplet
}

context_map = {
    'ELEMENT': parse_element,
    'FUNCTION': parse_function,
    'PHASE': phase_context_map,
}

def _setitem_raise_duplicates(dictionary, key, value):
    if key in dictionary:
        raise ValueError("Database contains duplicate FUNCTION {}".format(key))
    dictionary[key] = value


def read_ion(dbf, fd):
    lines = ion.load(fd, single_value=False)
    for line in lines:
        linetype = line.ion_annotations[0].text.upper()
        if linetype == 'ELEMENT':
            element = parse_element(line.text)
            dbf.species.add(v.Species(element, {element: 1.0}, charge=0))
            dbf.elements.add(element)
        elif linetype == 'FUNCTION':
            if isinstance(line, ion.IonPyDict):
                for function_name, function_def in line.items():
                    function_obj = parse_function(function_def)
                    _setitem_raise_duplicates(dbf.symbols, function_name, function_obj)
            else:
                raise ValueError('Function must be defined inside dictionary')
        elif linetype == 'PHASE':
            if isinstance(line, ion.IonPyDict):
                for phase_name, phase_def in line.items():
                    if len(phase_def.ion_annotations) > 0:
                        phase_type = phase_def.ion_annotations[0].text.upper()
                    else:
                        phase_type = None
                    phase_obj = phase_context_map[phase_type](phase_def)
                    model_hints = {} # TODO
                    dbf.add_structure_entry(phase_name, phase_name)
                    dbf.add_phase(phase_name, model_hints, phase_obj['sites'])
                    dbf.add_phase_constituents(phase_name, phase_obj['constituent_array'])
            else:
                raise ValueError('Phase must be defined inside dictionary')
        elif linetype == 'PARAMETER':
            for param_line in line:
                param_type = param_line.ion_annotations[0].text.upper()
                param_context_map[param_type](dbf, param_line)
        elif linetype == 'REFERENCE':
            pass
        else:
            raise ValueError('Failed while parsing: ', line)
    dbf.process_parameter_queue()

Database.register_format("ion", read=read_ion, write=None)