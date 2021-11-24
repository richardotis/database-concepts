from lxml import etree

parser = etree.XMLParser(load_dtd=False,
                         no_network=True)

tree = etree.parse("experiment3.xml", parser=parser)

root = tree.getroot()

elements = []
phases = []
all_parameters = []
symbols = {}


def parse_parameter(param_node):
    if param_node.attrib['type'] == 'G':
        int_order = 0
        constituent_array = [s.xpath('./Constituent/@refid') for s in param_node.xpath('./ConstituentArray/Site')]
    elif param_node.attrib['type'] == 'L':
        int_order = int(param_node.xpath('./Order')[0].text)
        constituent_array = [s.xpath('./Constituent/@refid') for s in param_node.xpath('./ConstituentArray/Site')]
    else:
        raise ValueError('Unknown parameter')
    return int_order, constituent_array


def parse_model(phase_name, model_node, parameters):
    site_ratios = [float(m) for m in model_node.xpath('./ConstituentArray/Site/@ratio')]
    sublattice_model = [s.xpath('./Constituent/@refid') for s in model_node.xpath('./ConstituentArray/Site')]
    print(phase_name, site_ratios, sublattice_model)
    for param_node in parameters:
        int_order, constituent_array = parse_parameter(param_node)
        param_value = '+'.join(param_node.xpath('./Expr/@refid')) +\
                      ''.join(param_node.xpath('./text()')).strip()
        all_parameters.append((phase_name, param_node.attrib['type'], constituent_array, int_order, param_value))


for child in root:
    if child.tag == 'Element':
        elements.append(child.attrib['id'])
    elif child.tag == 'Expr':
        symbols[child.attrib['id']] = list(child)
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
