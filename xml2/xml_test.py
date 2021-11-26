from pycalphad import Database, Model, calculate, variables as v
import xml_parser

dbf = Database('experiment3.xml')
mod = Model(dbf, ['AL', 'ZN'], 'LIQUID')
res = calculate(dbf, ['AL', 'ZN'], 'LIQUID', T=600, P=1e5)
print(res)
