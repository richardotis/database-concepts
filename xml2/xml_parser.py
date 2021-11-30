from pycalphad.io.tdb import _sympify_string, _process_reference_state, to_interval
from pycalphad import Database, variables as v
from pycalphad import __version__ as pycalphad_version
from sympy import Piecewise, And, Symbol
from lxml import etree, objectify


def convert_math_to_symbolic(math_nodes):
    result = 0.0
    interval_nodes = [x for x in math_nodes if (not isinstance(x, str)) and x.tag == 'Interval']
    string_nodes = [x for x in math_nodes if isinstance(x, str)]
    for math_node in string_nodes:
        # +0 is a hack, for how the function works
        result += _sympify_string(math_node+'+0')
    result += convert_intervals_to_piecewise(interval_nodes)
    result = result.xreplace({Symbol('T'): v.T, Symbol('P'): v.P})
    return result


def convert_intervals_to_piecewise(interval_nodes):
    exprs = []
    conds = []
    for interval_node in interval_nodes:
        if interval_node.attrib['in'] != 'T':
            raise ValueError('Unsupported interval')
        variable = interval_node.attrib['in']
        lower = float(interval_node.attrib.get('lower', '-inf'))
        upper = float(interval_node.attrib.get('upper', 'inf'))
        math_expr = convert_math_to_symbolic([''.join(interval_node.itertext()).replace('\n', '').replace(' ', '').strip()])
        if upper != float('inf'):
            cond = And(lower <= getattr(v, variable, Symbol(variable)), upper > getattr(v, variable))
        else:
            cond = (lower <= getattr(v, variable, Symbol(variable)))
        conds.append(cond)
        exprs.append(math_expr)
    if len(exprs) == 0:
        return 0
    return Piecewise(*(list(zip(exprs, conds)) + [(0, True)]), evaluate=False)


def convert_symbolic_to_nodes(sym):
    nodes = []
    if isinstance(sym, Piecewise):
        for expr, cond in sym.args:
            interval = to_interval(cond)
            lower = str(float(interval.start))
            upper = str(float(interval.end))
            converted_expr_nodes = [x for x in convert_symbolic_to_nodes(expr) if x != '0']
            if lower == '-inf' and upper == 'inf':
                nodes.extend(converted_expr_nodes)
                continue
            elif lower != '-inf' and upper == 'inf':
                interval_node = etree.Element("Interval", attrib={"in": "T", "lower": lower})
            else:
                interval_node = etree.Element("Interval", attrib={"in": "T", "lower": lower, "upper": upper})
            for node in converted_expr_nodes:
                if isinstance(node, str):
                    interval_node.text = node
                else:
                    interval_node.append(node)
            nodes.append(interval_node)

    else:
        str_node = str(sym).replace('log(', 'ln(')
        nodes.append(str_node)
    return nodes


def parse_cef_parameter(param_node):
    order_nodes = param_node.xpath('./Order')
    if len(order_nodes) == 0:
        int_order = 0
    else:
        int_order = int(order_nodes[0].text)
    constituent_array = [t.xpath('./Constituent/@refid') for t in param_node.xpath('./ConstituentArray/Site')]
    return int_order, constituent_array

# Symmetry options handled separately in XML
phase_options = {'ionic_liquid_2SL': 'TwoSublatticeIonicLiquid',
                 'liquid': 'Liquid',
                 'gas': 'Gas',
                 'aqueous': 'Aqueous',
                 'charged_phase': 'Charged'}
inv_phase_options = dict([reversed(i) for i in phase_options.items()])


def parse_model(dbf, phase_name, model_node, parameters):
    site_ratios = [float(m) for m in model_node.xpath('./ConstituentArray/Site/@ratio')]
    sublattice_model = [s.xpath('./Constituent/@refid') for s in model_node.xpath('./ConstituentArray/Site')]

    model_hints = {}
    magnetic_ordering_nodes = model_node.xpath('./MagneticOrdering')
    for magnetic_ordering_node in magnetic_ordering_nodes:
        if magnetic_ordering_node.attrib['type'] == 'IHJ':
            model_hints['ihj_magnetic_afm_factor'] = float(magnetic_ordering_node.attrib['afm_factor'])
            model_hints['ihj_magnetic_structure_factor'] = float(magnetic_ordering_node.attrib['structure_factor'])
        else:
            raise ValueError('Unknown magnetic ordering model')
    atomic_ordering_nodes = model_node.xpath('./AtomicOrdering')
    for atomic_ordering_node in atomic_ordering_nodes:
        model_hints['ordered_phase'] = str(atomic_ordering_node.attrib['ordered_part'])
        model_hints['disordered_phase'] = str(atomic_ordering_node.attrib['disordered_part'])
    # Simple phase options
    for ipo in inv_phase_options.keys():
        ipo_nodes = model_node.xpath('./'+str(ipo))
        for ipo_node in ipo_nodes:
            model_hints[inv_phase_options[ipo_node.tag]] = True
    # Parameter symmetry options
    symmetry_nodes = model_node.xpath('./Symmetry')
    for symmetry_node in symmetry_nodes:
        model_hints['symmetry_'+str(symmetry_node.attrib['type'])] = True
    dbf.add_structure_entry(phase_name, phase_name)
    dbf.add_phase(phase_name, model_hints, site_ratios)
    dbf.add_phase_constituents(phase_name, sublattice_model)

    for param_node in parameters:
        int_order, constituent_array = parse_cef_parameter(param_node)
        param_nodes = param_node.xpath('./Interval') + [''.join(param_node.xpath('./text()')).strip()]
        function_obj = convert_math_to_symbolic(param_nodes)
        param_type = param_node.attrib['type']
        ref = None  # TODO
        diffusing_species_nodes = model_node.xpath('./DiffusingSpecies/@refid')
        if len(diffusing_species_nodes) == 1:
            diffusing_species = str(diffusing_species_nodes[0])
        elif len(diffusing_species_nodes) == 0:
            diffusing_species = None
        else:
            raise ValueError('Multiple DiffusingSpecies nodes specified')
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
    relaxng = etree.RelaxNG(etree.parse('database.rng'))
    if not relaxng.validate(tree):
        print('Validation Error:', relaxng.error_log)
    root = tree.getroot()

    for child in root:
        if child.tag == 'ChemicalElement':
            element = str(child.attrib['id'])
            dbf.species.add(v.Species(element, {element: 1}, charge=0))
            dbf.elements.add(element)
            _process_reference_state(dbf, element, child.attrib['reference_phase'],
                                     float(child.attrib['mass']), float(child.attrib['H298']), float(child.attrib['S298']))
        elif child.tag == 'Species':
            species = str(child.attrib['id'])
            constituent_dict = {}
            species_charge = float(child.attrib.get('charge', 0))
            constituent_nodes = child.xpath('./ChemicalElement')
            for constituent_node in constituent_nodes:
                el = constituent_node.attrib['refid']
                ratio = float(constituent_node.attrib['ratio'])
                constituent_dict[el] = ratio
            dbf.species.add(v.Species(species, constituent_dict, charge=species_charge))
        elif child.tag == 'Expr':
            function_name = str(child.attrib['id'])
            function_obj = convert_intervals_to_piecewise(child)
            _setitem_raise_duplicates(dbf.symbols, function_name, function_obj)
        elif child.tag == 'Phase':
            model_nodes = child.xpath('./Model')
            if len(model_nodes) == 0:
                continue
            model_node = model_nodes[0]
            if model_node.attrib['type'] != 'CEF':
                continue
            phase_name = child.attrib['id']
            parameters = child.xpath('./Parameter')
            parse_model(dbf, phase_name, model_node, parameters)
    dbf.process_parameter_queue()


def write_xml(dbf, fd):
    root = objectify.Element("Database", version=str(0))
    metadata = objectify.SubElement(root, "metadata")
    writer = objectify.SubElement(metadata, "writer")
    writer._setText('pycalphad ' + str(pycalphad_version))
    phase_nodes = {}
    for element in sorted(dbf.elements):
        ref = dbf.refstates.get(element, {})
        refphase = ref.get('phase', 'BLANK')
        mass = ref.get('mass', 0.0)
        H298 = ref.get('H298', 0.0)
        S298 = ref.get('S298', 0.0)
        objectify.SubElement(root, "ChemicalElement", id=str(element), mass=str(mass),
                             reference_phase=refphase, H298=str(H298), S298=str(S298))
    for species in sorted(dbf.species, key=lambda s: s.name):
        if species.name not in dbf.elements:
            species_node = objectify.SubElement(root, "Species", id=str(species.name), charge=str(species.charge))
            for el_name, ratio in sorted(species.constituents.items(), key=lambda t: t[0]):
                objectify.SubElement(species_node, "ChemicalElement", refid=str(el_name), ratio=str(ratio))
    for name, expr in sorted(dbf.symbols.items()):
        expr_node = objectify.SubElement(root, "Expr", id=str(name))
        converted_nodes = convert_symbolic_to_nodes(expr)
        for node in converted_nodes:
            if isinstance(node, str):
                expr_node._setText(node)
            else:
                expr_node.append(node)
    for name, phase_obj in sorted(dbf.phases.items()):
        if phase_nodes.get(name, None) is None:
            phase_nodes[name] = objectify.SubElement(root, "Phase", id=str(name))
        # All model hints must be consumed for the writing to be considered successful
        model_hints = phase_obj.model_hints.copy()
        possible_options = set(phase_options.keys()).intersection(model_hints.keys())
        # TODO: detection for MQMQA, etc.
        if True:
            model_node = objectify.SubElement(phase_nodes[name], "Model", type="CEF")
            constit_array_node = objectify.SubElement(model_node, "ConstituentArray")
            subl_idx = 0
            for site_ratio, constituents in zip(phase_obj.sublattices, phase_obj.constituents):
                site_node = objectify.SubElement(constit_array_node, "Site", id=str(subl_idx), ratio=str(site_ratio))
                for constituent in sorted(constituents, key=str):
                    objectify.SubElement(site_node, "Constituent", refid=str(constituent))
                subl_idx += 1
            # IHJ model
            if 'ihj_magnetic_afm_factor' in model_hints.keys():
                objectify.SubElement(model_node, "MagneticOrdering",
                    type="IHJ", structure_factor=str(model_hints['ihj_magnetic_structure_factor']),
                    afm_factor=str(model_hints['ihj_magnetic_afm_factor']))
                del model_hints['ihj_magnetic_afm_factor']
                del model_hints['ihj_magnetic_structure_factor']
            # Two-part atomic ordering
            if ('ordered_phase' in model_hints.keys()):
                objectify.SubElement(model_node, "AtomicOrdering",
                    ordered_part=str(model_hints['ordered_phase']),
                    disordered_part=str(model_hints['disordered_phase']))
                del model_hints['ordered_phase']
                del model_hints['disordered_phase']
            # Symmetry options
            symmetry_node = None
            if ('symmetry_FCC_4SL' in model_hints.keys()):
                symmetry_node = objectify.SubElement(model_node, "Symmetry", type="FCC_4SL")
                del model_hints['symmetry_FCC_4SL']
            if ('symmetry_BCC_4SL' in model_hints.keys()):
                if symmetry_node is not None:
                    raise ValueError('Multiple parameter symmetry options specified')
                del model_hints['symmetry_BCC_4SL']
            # Simple phase options
            for possible_option in possible_options:
                objectify.SubElement(model_node, phase_options[possible_option])
                del model_hints[possible_option]
        if len(model_hints) > 0:
            # Some model hints were not properly consumed
            raise ValueError('Not all model hints are supported: {}'.format(model_hints))

    for param in dbf._parameters.all():
        phase_name = param['phase_name']
        # Create phase implicitly if not defined
        if phase_nodes.get(phase_name, None) is None:
            phase_nodes[phase_name] = objectify.SubElement(root, "Phase", id=str(phase_name))
        phase_node = phase_nodes[phase_name]
        param_node = objectify.SubElement(phase_node, "Parameter", type=str(param['parameter_type']))
        order_node = objectify.SubElement(param_node, "Order")
        order_node._setText(str(param['parameter_order']))
        constit_array_node = objectify.SubElement(param_node, "ConstituentArray")
        subl_idx = 0
        for constituents in param['constituent_array']:
            site_node = objectify.SubElement(constit_array_node, "Site", refid=str(subl_idx))
            for constituent in constituents:
                objectify.SubElement(site_node, "Constituent", refid=str(constituent))
            subl_idx += 1
        if param['diffusing_species'] != v.Species(None):
            objectify.SubElement(param_node, "DiffusingSpecies", refid=str(param['diffusing_species']))
        nodes = convert_symbolic_to_nodes(param['parameter'])
        for node in nodes:
            if isinstance(node, str):
                param_node._setText(node)
            else:
                param_node.append(node)
        # TODO: param['reference']
    objectify.deannotate(root, xsi_nil=True)
    etree.cleanup_namespaces(root)
    fd.write('<?xml version="1.0"?>\n')
    # XXX: href needs to be changed
    fd.write('<?xml-model href="database.rng" schematypens="http://relaxng.org/ns/structure/1.0" type="application/xml"?>\n')
    fd.write(etree.tostring(root, pretty_print=True).decode("utf-8"))


Database.register_format("xml", read=read_xml, write=write_xml)
