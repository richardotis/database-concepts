from pycalphad.io.tdb import _sympify_string
from pycalphad import variables as v
from sympy import Piecewise, And, Symbol
from lxml import etree

parser = etree.XMLParser(load_dtd=False,
                         no_network=True)

tree = etree.parse("experiment3.xml", parser=parser)

root = tree.getroot()

elements = []
phases = []
all_parameters = []
symbols = {}


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
                                             [''.join((interval_node.text or '').replace('\n', '').strip())])
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


def parse_model(phase_name, model_node, parameters):
    site_ratios = [float(m) for m in model_node.xpath('./ConstituentArray/Site/@ratio')]
    sublattice_model = [s.xpath('./Constituent/@refid') for s in model_node.xpath('./ConstituentArray/Site')]
    print(phase_name, site_ratios, sublattice_model)
    for param_node in parameters:
        int_order, constituent_array = parse_parameter(param_node)
        param_nodes = param_node.xpath('./Expr/@refid') + [''.join(param_node.xpath('./text()')).strip()]
        param_nodes = convert_math_to_symbolic(param_nodes)
        all_parameters.append((phase_name, param_node.attrib['type'], constituent_array, int_order, param_nodes))


for child in root:
    if child.tag == 'Element':
        elements.append(child.attrib['id'])
    elif child.tag == 'Expr':
        symbols[child.attrib['id']] = convert_intervals_to_piecewise(child)
    elif child.tag == 'Phase':
        model_node = child.xpath('./Model')[0]
        if model_node.attrib['type'] != 'CEF':
            continue
        phase_name = child.attrib['id']
        parameters = child.xpath('./Parameter')
        phases.append(phase_name)
        parse_model(phase_name, model_node, parameters)

print('Elements:', elements)
print('Symbols: ', symbols)
print('Phases: ', phases)
print('Parameters: ', all_parameters)

#etree.dump(tree.getroot())
