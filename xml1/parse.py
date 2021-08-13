from lxml import etree

parser = etree.XMLParser(load_dtd=False,
                         no_network=True)

tree = etree.parse("experiment.xml", parser=parser)

root = tree.getroot()

elements = []
phases = []
parameters = []
symbols = {}

def parse_phase(phase_node):
    for child in phase_node:
        if child.tag == 'Parameter':
            parameters.append((child.attrib['type'], child.attrib['constituents'], child.attrib['ord'], list(child)))

for child in root:
    if child.tag == 'Element':
        elements.append(child.attrib['id'])
    if child.tag == 'Expr':
        symbols[child.attrib['id']] = list(child)
    if child.tag == 'Phase':
        phases.append(child.attrib['name'])
        parse_phase(child)

print('Elements:', elements)
print('Symbols: ', symbols)
print('Phases: ', phases)
print('Parameters: ', parameters)

#etree.dump(tree.getroot())