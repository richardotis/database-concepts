from pycalphad.io.tdb import _sympify_string
from pycalphad import Database, variables as v
from sympy import Piecewise, And, Symbol
from lxml import etree


def convert_math_to_symbolic(math_nodes):
    result = 0.0
    for math_node in math_nodes:
        if isinstance(math_node, str):
            # +0 is a hack, for how the function works
            result += _sympify_string(math_node+'+0')
        elif math_node.tag == 'Expr':
            result += Symbol(math_node.xpath('./@(refid or id)')[0])
    result = result.xreplace({Symbol('T'): v.T, Symbol('P'): v.P})
    return result


def convert_intervals_to_piecewise(interval_nodes):
    exprs = []
    conds = []
    for interval_node in interval_nodes:
        if interval_node.attrib['in'] != 'T':
            raise ValueError('Unsupported interval')
        variable = interval_node.attrib['in']
        lower = float(interval_node.attrib['lower'])
        upper = float(interval_node.attrib['upper'])
        math_expr = convert_math_to_symbolic(interval_node.xpath('./Expr/@refid') + \
                                             [''.join(interval_node.itertext()).replace('\n', '').replace(' ', '').strip()])
        cond = And(lower <= getattr(v, variable, Symbol(variable)), upper > getattr(v, variable))
        conds.append(cond)
        exprs.append(math_expr)
    return Piecewise(*(list(zip(exprs, conds)) + [(0, True)]))


def parse_parameter(param_node):
    if param_node.attrib['type'] == 'G':
        int_order = 0
        constituent_array = [t.xpath('./Constituent/@refid') for t in param_node.xpath('./ConstituentArray/Site')]
    elif param_node.attrib['type'] == 'L':
        int_order = int(param_node.xpath('./Order')[0].text)
        constituent_array = [t.xpath('./Constituent/@refid') for t in param_node.xpath('./ConstituentArray/Site')]
    else:
        raise ValueError('Unknown parameter')
    return int_order, constituent_array


def parse_model(dbf, phase_name, model_node, parameters):
    site_ratios = [float(m) for m in model_node.xpath('./ConstituentArray/Site/@ratio')]
    sublattice_model = [s.xpath('./Constituent/@refid') for s in model_node.xpath('./ConstituentArray/Site')]

    model_hints = {}  # TODO
    dbf.add_structure_entry(phase_name, phase_name)
    dbf.add_phase(phase_name, model_hints, site_ratios)
    dbf.add_phase_constituents(phase_name, sublattice_model)

    for param_node in parameters:
        int_order, constituent_array = parse_parameter(param_node)
        param_nodes = param_node.xpath('./Expr/@refid') + [''.join(param_node.xpath('./text()')).strip()]
        function_obj = convert_math_to_symbolic(param_nodes)
        param_type = param_node.attrib['type']
        ref = None  # TODO
        diffusing_species = None  # TODO
        dbf.add_parameter(param_type, phase_name,
                          [[str(c) for c in sorted(lx)] for lx in constituent_array],
                          int_order, function_obj, ref, diffusing_species, force_insert=False)


def _setitem_raise_duplicates(dictionary, key, value):
    if key in dictionary:
        raise ValueError("Database contains duplicate FUNCTION {}".format(key))
    dictionary[key] = value


def read_xml(dbf, fd):
    parser = etree.XMLParser(load_dtd=False,
                             no_network=True)
    tree = etree.parse(fd, parser=parser)
    root = tree.getroot()

    for child in root:
        if child.tag == 'ChemicalElement':
            element = str(child.attrib['id'])
            dbf.species.add(v.Species(element, {element: 1.0}, charge=0))
            dbf.elements.add(element)
        elif child.tag == 'Expr':
            function_name = str(child.attrib['id'])
            function_obj = convert_intervals_to_piecewise(child)
            _setitem_raise_duplicates(dbf.symbols, function_name, function_obj)
        elif child.tag == 'Phase':
            model_node = child.xpath('./Model')[0]
            if model_node.attrib['type'] != 'CEF':
                continue
            phase_name = child.attrib['id']
            parameters = child.xpath('./Parameter')
            parse_model(dbf, phase_name, model_node, parameters)
    dbf.process_parameter_queue()


Database.register_format("xml", read=read_xml, write=None)
